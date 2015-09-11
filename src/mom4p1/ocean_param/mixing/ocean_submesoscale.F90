!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ocean_submesoscale_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes a streamfunction within
! the upper surface boundary layer, and applies this
! streamfunction to all tracers.  
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes a streamfunction within
! the upper surface boundary layer, and applies this
! streamfunction to all tracers.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Fox-Kemper, Ferrari, and Hallberg 2008: Parameterization of 
! mixed layer eddies. Part I: theory and diagnosis
! Journal of Physical Oceanography, in press. 
! </REFERENCE>
!
! <REFERENCE>
! Fox-Kemper, Danabasoglu, Ferrari, and Hallberg 2008: 
! Parameterizing submesoscale physics in global models.
! Clivar Exchanges, vol 13, no.1,  Jan2008. pages 3-5.
! </REFERENCE>
! </INFO>
!
!<NAMELIST NAME="ocean_submesoscale_nml">
!  <DATA NAME="use_this_module=" TYPE="logical">
!  Must be .true. to use this module.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!  <DATA NAME="diag_step" TYPE="integer">
!  Number of time steps between computing max bottom value for wrho_bt_submeso.
!  Default diag_step=1200.
!  </DATA> 
!
!  <DATA NAME="submeso_skew_flux" TYPE="logical">
!  For computing the tendency as convergence of skew flux.
!  Default submeso_skew_flux=.true.
!  </DATA> 
!  <DATA NAME="submeso_advective_flux" TYPE="logical">
!  For computing the tendency as convergence of advective flux.
!  This approach uses flux limited sweby advection, which ensures
!  that the resulting tendency will not create extrema in the 
!  tracer field.  
!  This option has a bug, and so cannot be used.  
!  Default submeso_advective_flux=.false.
!  </DATA> 
!  <DATA NAME="submeso_diag_advect_transport" TYPE="logical">
!  For diagnosing the advective mass transport.
!  Doing so requires a call to a subroutine.  If not aiming 
!  to diagnose the velocity, then saves time to set 
!  submeso_diag_advect_transport=.false. 
!  Default submeso_diag_advect_transport=.false.
!  </DATA> 
!
!  <DATA NAME="use_hblt_constant" TYPE="logical">
!  For running with a constant boundary layer depth. This for the case when 
!  not using a realistic mixed layer scheme.  Default use_hblt_constant=.false.
!  </DATA> 
!  <DATA NAME="constant_hblt" UNITS="metre" TYPE="real">
!  The boundary layer depth for the case when use_hblt_constant=.true.
!  Default constant_hblt=100.0.
!  </DATA> 
!  <DATA NAME="use_hblt_equal_mld" TYPE="logical">
!  For using the diagnosed mld as the hblt for submeso.  
!  This is useful for those test models that do not have a mixed layer
!  scheme enabled, such as KPP, where the mixed layer scheme provides a
!  boundary layer depth.  In this case, it is sensible to employ the diagnosed
!  mixed layer depth for the submeso scheme. Additionally, in general it is 
!  more physical to use the mld than the KPP hblt as the depth over which 
!  the submesoscale eddies act.  Hence, default use_hblt_equal_mld=.true.
!  </DATA> 
!  <DATA NAME="min_kblt" UNITS="dimensionless" TYPE="integer">
!  The minimum number of vertical cells in the surface boundary layer 
!  that are required in order to compute the submesoscale streamfunction.
!  Default min_kblt=3.  Need at least three to fit a parabola with zero 
!  streamfunction at the top and bottom of the boundary layer.  
!  </DATA> 
!
!  <DATA NAME="minimum_hblt" TYPE="real" UNITS="metre">
!  For setting a floor to the hblt used for submesoscale scheme. 
!  Default minimum_hblt=0.0.
!  </DATA> 
!  <DATA NAME="smooth_hblt" TYPE="logical">
!  For smoothing on the hblt field. This is useful since the hblt 
!  obtained from KPP or diagnosed mld can have some grid noise. 
!  Default smooth_hblt=.false.
!  </DATA> 
!  <DATA NAME="vel_micom_smooth" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the Laplacian mixing 
!  coefficient used in the Laplacian smoothing of hblt. 
!  Default vel_micom_smooth=0.2. 
!  </DATA>
!
!  <DATA NAME="smooth_psi" TYPE="logical">
!  For doing a horizontal 1-2-1 smoothing on the psix and psiy fields. 
!  This is useful to reduce noise. Default smooth_psi=.true.
!  </DATA> 
!  <DATA NAME="limit_psi" TYPE="logical">
!  For limiting the magnitude of psi in order to reduce possibility of 
!  model crashes.  Default limit_psi=.true.
!  </DATA> 
!  <DATA NAME="limit_psi_velocity_scale" UNITS="metre/sec" TYPE="real">
!  Velocity scale used to limit the value of psi when limit_psi=.true.
!  Default limit_psi_velocity_scale=0.5.
!  </DATA> 
!  <DATA NAME="submeso_limit_flux" TYPE="logical">
!  For limiting the fluxes arising from submeso scheme, according to 
!  tmask_limit. When reach a point where tmask_limit=1.0, then set
!  the submeso flux for this cell to zero. 
!  Default submeso_limit_flux=.true.
!  </DATA> 
!
! <DATA NAME="coefficient_ce" UNITS="dimensionless" TYPE="real">
!  The dimensionless coefficient from the Fox-Kemper etal scheme. 
!  They recommend setting coefficient_ce between 0.06 and 0.08.  
!  Default coefficient_ce=0.07.
!  </DATA> 
! <DATA NAME="time_constant" UNITS="seconds" TYPE="real">
!  Timescale to mix momentum across the mixed layer.  
!  Default time_constant=86400.0 = 1day. 
!  </DATA> 
! <DATA NAME="front_length_const" UNITS="metre" TYPE="real">
!  Take constant horizontal length scale of submesoscale front. 
!  Default front_length_const=5e3.
!  </DATA> 
! <DATA NAME="front_length_deform_radius" TYPE="logical">
!  To compute the front length using the mixed layer deformation 
!  radius. Default front_length_deform_radius=.true.  Note, 
!  will have a floor on the variable front length set by the
!  nml setting for front_length_const.
!  </DATA> 
!
!  <DATA NAME="transport_units" TYPE="character">
!  The units for writing out the transport.  Either in 
!  Sv (10^9 kg/s) or mks (kg/s). Note the mks unit is requested 
!  for CMIP5 purposes.
!  Default transport_units = 'Sv'. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,     only: epsln, grav
use diag_manager_mod,  only: register_diag_field, register_static_field, need_data, send_data
use fms_mod,           only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,           only: stdout, stdlog, read_data, NOTE, FATAL, WARNING
use mpp_domains_mod,   only: mpp_update_domains, XUPDATE, YUPDATE, CGRID_NE
use mpp_mod,           only: mpp_error, mpp_chksum, mpp_max, mpp_pe

use ocean_domains_mod,    only: get_local_indices, set_ocean_domain
use ocean_operators_mod,  only: LAP_T, BDX_ET, BDY_NT
use ocean_parameters_mod, only: missing_value, DEPTH_BASED, omega_earth
use ocean_parameters_mod, only: rho0, rho0r, onehalf, onesixth, oneeigth
use ocean_tracer_diag_mod,only: calc_mixed_layer_depth
use ocean_types_mod,      only: tracer_2d_type, tracer_3d_0_nk_type, tracer_3d_1_nk_type
use ocean_types_mod,      only: ocean_time_type, ocean_domain_type, ocean_grid_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_thickness_type, ocean_density_type
use ocean_util_mod,       only: write_timestamp
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk1_2d, wrk1_v

implicit none

private

! for diagnostics 
integer :: id_kblt_submeso         =-1
integer :: id_hblt_submeso         =-1
integer :: id_mu_submeso           =-1
integer :: id_psix_submeso         =-1
integer :: id_psiy_submeso         =-1
integer :: id_tx_trans_submeso     =-1
integer :: id_ty_trans_submeso     =-1
integer :: id_neutral_rho_submeso  =-1
integer :: id_wdian_rho_submeso    =-1
integer :: id_front_length_submeso =-1
integer :: id_buoy_freq_ave_submeso=-1 
integer :: id_uhrho_et_submeso     =-1
integer :: id_vhrho_nt_submeso     =-1
integer :: id_wrho_bt_submeso      =-1
integer :: id_u_et_submeso         =-1
integer :: id_v_nt_submeso         =-1
integer :: id_w_bt_submeso         =-1

integer, dimension(:), allocatable :: id_xflux_submeso       ! i-directed flux 
integer, dimension(:), allocatable :: id_yflux_submeso       ! j-directed flux 
integer, dimension(:), allocatable :: id_zflux_submeso       ! k-directed flux 
integer, dimension(:), allocatable :: id_xflux_submeso_int_z ! vertically integrated i-flux
integer, dimension(:), allocatable :: id_yflux_submeso_int_z ! vertically integrated j-flux
integer, dimension(:), allocatable :: id_submeso             ! tendency arising from submesoscale param 

logical :: used


#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()
type(ocean_domain_type), save    :: Dom_flux_sub

character(len=128)  :: version='$$'
character (len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'


#ifdef MOM4_STATIC_ARRAYS
real, dimension(isd:ied,jsd:jed,nk)   :: uhrho_et_submeso ! i-component of transport for submeso 
real, dimension(isd:ied,jsd:jed,nk)   :: vhrho_nt_submeso ! j-component of transport for submeso 
real, dimension(isd:ied,jsd:jed,0:nk) :: wrho_bt_submeso  ! vertical component of transport for submeso 

real, dimension(isd:ied,jsd:jed,nk,0:1) :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(isd:ied,jsd:jed,0:nk)   :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(isd:ied,jsd:jed,0:1)    :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(isd:ied,jsd:jed,nk,0:1) :: psix             ! streamfunction x-component (m^2/sec) 
real, dimension(isd:ied,jsd:jed,nk,0:1) :: psiy             ! streamfunction y-component (m^2/sec) 
real, dimension(isd:ied,jsd:jed)        :: hblt             ! boundary layer depth (m) 
real, dimension(isd:ied,jsd:jed)        :: grid_length      ! grid length scale (m)
real, dimension(isd:ied,jsd:jed)        :: front_length_inv ! inverse front length (1/m)
real, dimension(isd:ied,jsd:jed)        :: buoy_freq_ave    ! buoyancy frequency averaged over mixed layer depth (1/sec)
real, dimension(isd:ied,jsd:jed)        :: time_factor      ! time factor (sec) for computing streamfunction
real, dimension(isd:ied,jsd:jed)        :: coriolis_param   ! absolute value of the Coriolis parameter (sec^-1) on T-grid 
integer, dimension(isd:ied,jsd:jed)     :: kblt             ! k-level encompassing hblt 

real, dimension(isd:ied,jsd:jed,nk)     :: flux_x      ! i-component to tracer flux
real, dimension(isd:ied,jsd:jed,nk)     :: flux_y      ! j-component to tracer flux
real, dimension(isd:ied,jsd:jed,0:nk)   :: flux_z      ! k-component to tracer flux

#else 

real, dimension(:,:,:), allocatable :: uhrho_et_submeso  ! i-component of transport for submeso 
real, dimension(:,:,:), allocatable :: vhrho_nt_submeso  ! j-component of transport for submeso 
real, dimension(:,:,:), allocatable :: wrho_bt_submeso   ! vertical component of transport for submeso 

real, dimension(:,:,:,:), allocatable :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(:,:,:),   allocatable :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(:,:,:),   allocatable :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),   allocatable :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),   allocatable :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(:,:,:),   allocatable :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(:,:,:,:), allocatable :: psix             ! streamfunction x-component (m^2/sec) 
real, dimension(:,:,:,:), allocatable :: psiy             ! streamfunction y-component (m^2/sec) 
real, dimension(:,:),     allocatable :: hblt             ! boundary layer depth (m) 
real, dimension(:,:),     allocatable :: grid_length      ! grid length scale (m)
real, dimension(:,:),     allocatable :: front_length_inv ! inverse front length (1/m)
real, dimension(:,:),     allocatable :: buoy_freq_ave    ! buoyancy frequency averaged over mixed layer depth (1/sec)
real, dimension(:,:),     allocatable :: time_factor      ! time factor (sec) for computing streamfunction
real, dimension(:,:),     allocatable :: coriolis_param   ! absolute value of the Coriolis parameter (sec^-1) on T-grid 
integer, dimension(:,:),  allocatable :: kblt             ! k-level encompassing hblt 

real, dimension(:,:,:),   allocatable :: flux_x      ! i-component to tracer flux
real, dimension(:,:,:),   allocatable :: flux_y      ! j-component to tracer flux
real, dimension(:,:,:),   allocatable :: flux_z      ! k-component to tracer flux


#endif 

! for advecting tracers with sweby advection 
real, dimension(:,:,:), allocatable :: tmask_mdfl
real, dimension(:,:,:), allocatable :: tracer_mdfl
type  :: tracer_mdfl_type
  real, dimension(:,:,:), pointer :: field => NULL()
end type tracer_mdfl_type
type(tracer_mdfl_type),  dimension(:), allocatable  :: tracer_mdfl_all ! tracer array for sweby advection   
type(ocean_domain_type), save :: Dom_mdfl

type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdx ! tracer partial derivative (tracer/m)
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdy ! tracer partial derivative (tracer/m)
type(tracer_3d_0_nk_type), dimension(:), allocatable  :: dTdz ! tracer partial derivative (tracer/m)

public ocean_submesoscale_init 
public submeso_restrat
private tracer_derivs 
private compute_flux_x
private compute_flux_y
private compute_flux_z
private compute_psi
private compute_transport 
private compute_submeso_skewsion
private compute_submeso_advection
private maximum_bottom_w_submeso


integer :: index_temp=-1
integer :: index_salt=-1 
integer :: num_prog_tracers=0

! for diagnosing fluxes
real    :: flux_sign

! vertical coordinate 
integer :: vert_coordinate_class=1

! for output
integer :: unit=6

real    :: fiveover21
real    :: eightover21
real    :: dtime 
real    :: ce_grav_rho0r 
real    :: time_constant2_r
real    :: grav_rho0r
real    :: front_length_const_inv
real    :: front_length_max=1e10

logical :: module_is_initialized=.false.

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_units='Sv' 
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 


! nml parameters 
logical :: use_this_module               = .false.
logical :: debug_this_module             = .false. 
logical :: submeso_skew_flux             = .true.
logical :: submeso_advective_flux        = .false.
logical :: submeso_diag_advect_transport = .false.
logical :: use_hblt_constant             = .false.    
logical :: use_hblt_equal_mld            = .true. 
logical :: smooth_hblt                   = .false.    
logical :: smooth_psi                    = .true.    
logical :: front_length_deform_radius    = .true.
logical :: limit_psi                     = .true.
logical :: submeso_limit_flux            = .true.
real    :: limit_psi_velocity_scale      = 0.5
real    :: vel_micom_smooth              = 0.2 
real    :: coefficient_ce                = 0.07
real    :: time_constant                 = 86400.0 
real    :: front_length_const            = 5e3 
real    :: constant_hblt                 = 100.0
real    :: minimum_hblt                  = 0.0
integer :: min_kblt                      = 3
integer :: diag_step                     = 1200

namelist /ocean_submesoscale_nml/ use_this_module, debug_this_module, diag_step,     &
                                  use_hblt_constant, use_hblt_equal_mld,             &
                                  smooth_hblt, vel_micom_smooth, constant_hblt,      &
                                  coefficient_ce, time_constant, front_length_const, &
                                  min_kblt, minimum_hblt, smooth_psi,                &
                                  front_length_deform_radius,                        &
                                  limit_psi, limit_psi_velocity_scale,               &
                                  submeso_limit_flux,                                &
                                  submeso_skew_flux, submeso_advective_flux,         &
                                  submeso_diag_advect_transport, transport_units

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_submesoscale_init">
!
! <DESCRIPTION>
! Initialization for the ocean_submesoscale module.
! </DESCRIPTION>
  subroutine ocean_submesoscale_init(Grid, Domain, Time, T_prog, Ocean_options, dtime_t, &
                                     ver_coordinate_class, debug)
  
    type(ocean_grid_type),        intent(in), target   :: Grid
    type(ocean_domain_type),      intent(in), target   :: Domain
    type(ocean_time_type),        intent(in)           :: Time
    type(ocean_prog_tracer_type), intent(in)           :: T_prog(:)
    type(ocean_options_type),     intent(inout)        :: Ocean_options 
    real,                         intent(in)           :: dtime_t
    integer,                      intent(in)           :: ver_coordinate_class 
    logical,                      intent(in), optional :: debug

    integer :: unit, io_status, ierr
    integer :: i,j,k,n
    integer :: num_methods=0

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
      call mpp_error(FATAL,&
      '==>Error from ocean_submesoscale_init_mod: module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    unit = open_namelist_file()
    read(unit, ocean_submesoscale_nml,iostat=io_status)
    write (stdoutunit,'(/)')
    write(stdoutunit,ocean_submesoscale_nml)    
    write(stdlogunit,ocean_submesoscale_nml)
    ierr = check_nml_error(io_status, 'ocean_submesoscale_nml')
    call close_file(unit)

    Dom => Domain
    Grd => Grid

#ifndef MOM4_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
#endif

    if(transport_units=='Sv') then
        transport_convert=1.0e-9 
        transport_dims   = 'Sv (10^9 kg/s)'
    else
        transport_convert=1.0
        transport_dims   = 'kg/s'
    endif

    if (PRESENT(debug) .and. .not. debug_this_module) then
        debug_this_module = debug
    endif
    if(debug_this_module) then 
        write(stdoutunit,'(a)') '==>Note: running ocean_submesoscale_mod with debug_this_module=.true.'  
    endif

    if(use_this_module) then 
      call mpp_error(NOTE, '==>Note: USING ocean_submesoscale_mod')
      Ocean_options%submesoscale = 'Used submesoscale closure for surface restratification.'
    else 
      call mpp_error(NOTE, '==>Note: NOT using ocean_submesoscale_mod')
      Ocean_options%submesoscale = 'Did NOT use submesoscale closure for surface restratification.'
      return 
    endif 

    if(use_hblt_equal_mld) then 
      write(stdoutunit,'(a)') &
      '==>Note: For ocean_submesoscale, setting bldepth equal to diagnosed mld.'
    endif 
    if(use_hblt_constant) then 
      write(stdoutunit,'(a)') &
      '==>Note: For ocean_submesoscale, setting bldepth equal to prescribed constant.'
    endif 
    if(use_hblt_equal_mld .and. use_hblt_constant) then 
      call mpp_error(FATAL, &
      '==>Error: in ocean_submesoscale, use_hblt_equal_mld & use_hblt_constant cannot both be true.')
    endif 

    if(submeso_advective_flux) then 
      write(stdoutunit,'(a)') &
      '==>Error: For ocean_submesoscale, advective flux calculation remains incomplete.  Use skew flux instead.'
      call mpp_error(FATAL, &
      '==>Error: in ocean_submesoscale, advective flux calculation remains incomplete.  Use skew flux instead.')
      num_methods=num_methods+1
      flux_sign = 1.0
    endif 
    if(submeso_skew_flux) then 
      write(stdoutunit,'(a)') &
      '==>Note: For ocean_submesoscale, computing tendency as convergence of skew flux.'
      num_methods=num_methods+1
      flux_sign = -1.0
    endif 
    if(num_methods > 1) then 
      call mpp_error(FATAL, &
      '==>Error: in ocean_submesoscale, can choose only one method for computing tendency.')
    endif 
    if(num_methods == 0) then 
      call mpp_error(FATAL, &
      '==>Error: in ocean_submesoscale, must choose one method for computing tendency.')
    endif 
   

    fiveover21             = 5.0/21.0
    eightover21            = 8.0/21.0
    grav_rho0r             = grav*rho0r
    ce_grav_rho0r          = coefficient_ce*grav*rho0r
    time_constant2_r       = 1.0/(time_constant**2) 
    dtime                  = dtime_t 
    vert_coordinate_class  = ver_coordinate_class
    front_length_const_inv = 1.0/front_length_const

    call set_ocean_domain(Dom_flux_sub, Grid,xhalo=Dom%xhalo,yhalo=Dom%yhalo,name='flux dom submeso',&
                          maskmap=Dom%maskmap)

    num_prog_tracers = size(T_prog(:))
    do n=1,num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo

    allocate( dTdx(num_prog_tracers) )
    allocate( dTdy(num_prog_tracers) )
    allocate( dTdz(num_prog_tracers) )

#ifndef MOM4_STATIC_ARRAYS
    allocate (uhrho_et_submeso(isd:ied,jsd:jed,nk))
    allocate (vhrho_nt_submeso(isd:ied,jsd:jed,nk))
    allocate (wrho_bt_submeso(isd:ied,jsd:jed,0:nk))

    allocate (delqc(isd:ied,jsd:jed,nk,0:1))
    allocate (dzwtr(isd:ied,jsd:jed,0:nk))
    allocate (dtew(isd:ied,jsd:jed,0:1))
    allocate (dtns(isd:ied,jsd:jed,0:1))
    allocate (dtwedyt(isd:ied,jsd:jed,0:1))
    allocate (dxtdtsn(isd:ied,jsd:jed,0:1))

    allocate (psix(isd:ied,jsd:jed,nk,0:1))
    allocate (psiy(isd:ied,jsd:jed,nk,0:1))
    allocate (hblt(isd:ied,jsd:jed))
    allocate (grid_length(isd:ied,jsd:jed))
    allocate (front_length_inv(isd:ied,jsd:jed))
    allocate (buoy_freq_ave(isd:ied,jsd:jed))
    allocate (time_factor(isd:ied,jsd:jed))
    allocate (coriolis_param(isd:ied,jsd:jed))
    allocate (kblt(isd:ied,jsd:jed))

    allocate (flux_x(isd:ied,jsd:jed,nk) )
    allocate (flux_y(isd:ied,jsd:jed,nk) )
    allocate (flux_z(isd:ied,jsd:jed,0:nk) )

    do n=1,num_prog_tracers
       allocate ( dTdx(n)%field(isd:ied,jsd:jed,nk) )
       allocate ( dTdy(n)%field(isd:ied,jsd:jed,nk) )
       allocate ( dTdz(n)%field(isd:ied,jsd:jed,0:nk) )
    enddo

#endif 

    uhrho_et_submeso = 0.0
    vhrho_nt_submeso = 0.0
    wrho_bt_submeso  = 0.0

    do n=1,num_prog_tracers
       dTdx(n)%field(:,:,:) = 0.0
       dTdy(n)%field(:,:,:) = 0.0
       dTdz(n)%field(:,:,:) = 0.0
    enddo

    psix = 0.0
    psiy = 0.0
    kblt = 0
    hblt = 0.0

    flux_x = 0.0
    flux_y = 0.0
    flux_z = 0.0

    dzwtr = 0.0
    delqc = 0.0

    grid_length(:,:)      = 2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))
    front_length_inv(:,:) = front_length_const_inv
    buoy_freq_ave(:,:)    = 0.0
    coriolis_param(:,:)   = 2.0*omega_earth*abs(sin(Grd%phit(:,:)))

    dtew(:,:,0) = Grd%dtw(:,:)
    dtew(:,:,1) = Grd%dte(:,:)
    dtns(:,:,0) = Grd%dts(:,:)
    dtns(:,:,1) = Grd%dtn(:,:)

    dtwedyt(:,:,:) = 0.0
    dtwedyt(:,:,0) = Grd%dte(:,:)*Grd%dyt(:,:)
    do i=isc-1,iec
       dtwedyt(i,:,1) = Grd%dtw(i+1,:)*Grd%dyt(i+1,:)
    enddo

    dxtdtsn(:,:,:) = 0.0
    dxtdtsn(:,:,0) = Grd%dxt(:,:)*Grd%dtn(:,:)
    do j=jsc-1,jec
       dxtdtsn(:,j,1) = Grd%dxt(:,j+1)*Grd%dts(:,j+1)
    enddo

    do j=jsd,jed
       do i=isd,ied
          time_factor(i,j) = 1.0/(sqrt( time_constant2_r + coriolis_param(i,j)**2 ))
       enddo
    enddo


    ! for the sweby advection scheme 
    if(submeso_advective_flux) then 
        allocate(tracer_mdfl_all(num_prog_tracers))
        allocate(tmask_mdfl (isc-2:iec+2,jsc-2:jec+2,nk))    
        allocate(tracer_mdfl(isc-2:iec+2,jsc-2:jec+2,nk))    
        do n=1,num_prog_tracers
           allocate (tracer_mdfl_all(n)%field(isc-2:iec+2,jsc-2:jec+2,nk))
        enddo
        tmask_mdfl  = 0.0
        tracer_mdfl = 0.0 
        do n=1,num_prog_tracers
           tracer_mdfl_all(n)%field(:,:,:) = 0.0
        enddo
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 tmask_mdfl(i,j,k) = Grd%tmask(i,j,k)
              enddo
           enddo
        enddo
        call mpp_update_domains(tmask_mdfl,Dom_mdfl%domain2d)  
    endif


    ! diagnostics 
    id_kblt_submeso = register_diag_field ('ocean_model', 'kblt_submeso', &
         Grid%tracer_axes(1:2), Time%model_time,                          &
         'Number of k-levels in boundary layer for submesoscale closure', &
         'dimensionless', missing_value=missing_value, range=(/-1.0,1e6/))
    id_hblt_submeso = register_diag_field ('ocean_model', 'hblt_submeso', &
         Grid%tracer_axes(1:2), Time%model_time,                          &
         'Boundary layer depth used for submesoscale closure',            &
         'metre', missing_value=missing_value, range=(/-1.0,1e6/))
    id_front_length_submeso = register_diag_field ('ocean_model', 'front_length_submeso', &
         Grid%tracer_axes(1:2), Time%model_time,                                          &
         'Front length used for submesoscale closure',                                    &
         'metre', missing_value=missing_value, range=(/-1.0,1e6/))
    id_buoy_freq_ave_submeso = register_diag_field ('ocean_model', 'buoy_freq_ave_submeso',&
         Grid%tracer_axes(1:2), Time%model_time,                                           &
         'Buoyancy frequency averaged over depth of mixed layer for submesoscale closure', &
         '1/sec', missing_value=missing_value, range=(/-1.0,1e6/))
    id_mu_submeso = register_diag_field ('ocean_model', 'mu_submeso',   &
         Grd%tracer_axes(1:3), Time%model_time,                         &
         'vertical structure function for submesoscale streamfunction', &
         'dimensionless', missing_value=missing_value, range=(/-1.e2,1e2/))
    id_psix_submeso = register_diag_field ('ocean_model', 'psix_submeso', &
         Grd%tracer_axes_flux_y(1:3), Time%model_time,                    &
         'i-comp of submesoscale streamfunction',                         &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))
    id_psiy_submeso = register_diag_field ('ocean_model', 'psiy_submeso', &
         Grd%tracer_axes_flux_x(1:3), Time%model_time,                    &
         'j-comp of submesoscale streamfunction',                         &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))
    id_tx_trans_submeso = register_diag_field ('ocean_model', 'tx_trans_submeso', &
         Grd%tracer_axes_flux_x(1:3), Time%model_time,                            &
         'T-cell mass i-transport from submesoscale param',                       &
         trim(transport_dims), missing_value=missing_value, range=(/-1.e10,1e10/))
    id_ty_trans_submeso = register_diag_field ('ocean_model', 'ty_trans_submeso', &
         Grd%tracer_axes_flux_y(1:3), Time%model_time,                            &
         'T-cell mass j-transport from submesoscale param',                       &
         trim(transport_dims), missing_value=missing_value, range=(/-1.e10,1e10/))

    id_u_et_submeso = register_diag_field ('ocean_model', 'u_et_submeso', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                     &
       'i-component of submesoscale transport velocity',                  &
       'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

    id_v_nt_submeso = register_diag_field ('ocean_model', 'v_nt_submeso', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                     &
       'j-component of submesoscale transport velocity',                  &
       'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

    id_w_bt_submeso = register_diag_field ('ocean_model', 'w_bt_submeso', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                     &
       'vertical component of submesoscale transport velocity',           &
       'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

    id_uhrho_et_submeso = register_diag_field ('ocean_model', 'uhrho_et_submeso', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                             &
       'i-component of submesoscale transport',                                   &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

    id_vhrho_nt_submeso = register_diag_field ('ocean_model', 'vhrho_nt_submeso', &
        Grd%tracer_axes_flux_y(1:3), Time%model_time,                             &
       'j-component of submesoscale transport',                                   &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

    id_wrho_bt_submeso = register_diag_field ('ocean_model', 'wrho_bt_submeso', &
        Grd%vel_axes_wt(1:3), Time%model_time,                                  &
       'k-component of submesoscale transport',                                 &
       '(kg/m^3)*m/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))



  id_neutral_rho_submeso = register_diag_field ('ocean_model', 'neutral_rho_submeso', &
        Grd%tracer_axes(1:3), Time%model_time,                      &
        'update of neutral density from submeso param',             &
        'rho*rho_dz/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_wdian_rho_submeso = register_diag_field ('ocean_model', 'wdian_rho_submeso', &
       Grd%tracer_axes(1:3), Time%model_time,                                     &
       'dianeutral velocity component due to submeso closure',                    &
       'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

    allocate (id_xflux_submeso(num_prog_tracers))
    allocate (id_yflux_submeso(num_prog_tracers))
    allocate (id_zflux_submeso(num_prog_tracers))
    allocate (id_xflux_submeso_int_z(num_prog_tracers))
    allocate (id_yflux_submeso_int_z(num_prog_tracers))
    allocate (id_submeso(num_prog_tracers))

    do n=1,num_prog_tracers

       if(n == index_temp) then 

           id_xflux_submeso(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_xflux_submeso',              &
                Grd%tracer_axes_flux_x(1:3), Time%model_time,        &
                'cp*submeso_xflux*dyt*rho_dzt*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_yflux_submeso(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_yflux_submeso',              &
                Grd%tracer_axes_flux_y(1:3), Time%model_time,        &
                'cp*submeso_yflux*dxt*rho_dzt*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_zflux_submeso(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_zflux_submeso',              &
                Grd%tracer_axes_wt(1:3), Time%model_time,            &
                'cp*submeso_zflux*dxt*dyt*rho*temp',                 &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_xflux_submeso_int_z(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_xflux_submeso_int_z',              &
                Grd%tracer_axes_flux_y(1:2), Time%model_time,              &
                'z-integral cp*submeso_xflux*dyt*rho_dzt*temp',            &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_yflux_submeso_int_z(n) = register_diag_field ('ocean_model', &
                trim(T_prog(n)%name)//'_yflux_submeso_int_z',              &
                Grd%tracer_axes_flux_y(1:2), Time%model_time,              &
                'z-integral cp*submeso_yflux*dxt*rho_dzt*temp',            &
                'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
           id_submeso(n) = register_diag_field ('ocean_model',           &
                trim(T_prog(n)%name)//'_submeso',                        &              
                Grd%tracer_axes(1:3), Time%model_time,                   &
                'rho*dzt*cp*submesoscale tendency (heating)',            &
                trim(T_prog(n)%flux_units), missing_value=missing_value, &
                range=(/-1.e10,1.e10/))

       else

           id_xflux_submeso(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_xflux_submeso',                       &
                Grd%tracer_axes_flux_x(1:3), Time%model_time,                 &
                'submeso_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
                trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,&
                range=(/-1.e18,1.e18/))
           id_yflux_submeso(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_yflux_submeso',                       &
                Grd%tracer_axes_flux_y(1:3), Time%model_time,                 &
                'submeso_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
                trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,&
                range=(/-1.e18,1.e18/))
           id_zflux_submeso(n) = register_diag_field ('ocean_model',          &
                trim(T_prog(n)%name)//'_zflux_submeso',                       &
                Grd%tracer_axes_wt(1:3), Time%model_time,                     &
                'submeso_yflux*dxt*dyt*rho*tracer for'//trim(T_prog(n)%name), &
                trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,&
                range=(/-1.e18,1.e18/))
           id_xflux_submeso_int_z(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_xflux_submeso_int_z',                           &
                Grd%tracer_axes_flux_x(1:2), Time%model_time,                           &
                'z-integral submeso_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name),&
                trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,          &
                range=(/-1.e18,1.e18/))
           id_yflux_submeso_int_z(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_yflux_submeso_int_z',                           &
                Grd%tracer_axes_flux_y(1:2), Time%model_time,                           &
                'z-integral submeso_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),&
                trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,          &
                range=(/-1.e18,1.e18/))
           id_submeso(n) = register_diag_field ('ocean_model',              &
                trim(T_prog(n)%name)//'_submeso',                           &
                Grd%tracer_axes(1:3), Time%model_time,                      &
                'rho*dzt*submesoscale tendency for '//trim(T_prog(n)%name), &
                trim(T_prog(n)%flux_units), missing_value=missing_value,    &
                range=(/-1.e10,1.e10/))

       endif

    enddo


end subroutine ocean_submesoscale_init
! </SUBROUTINE> NAME="ocean_submesoscale_init"


!#######################################################################
! <SUBROUTINE NAME="submeso_restrat">
!
! <DESCRIPTION>
! This routine computes a thickness and density weighted time tendency
! for each tracer, arising from the effects of parameterized 
! submesoscale eddies in the surface mixed layer.  
! </DESCRIPTION>
!
  subroutine submeso_restrat(Time, Thickness, Dens, T_prog, surf_blthick)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(in)    :: surf_blthick

  integer :: i,j,k,kp1,n
  integer :: ip,jq
  integer :: tau, taum1
  real    :: temporary

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_submesoscale_mode (ocean_submeso): needs initialization')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1 

  ! time dependent delqc geometric factor 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           delqc(i,j,k,0) = Grd%fracdz(k,0)*Thickness%rho_dzt(i,j,k,tau)
           delqc(i,j,k,1) = Grd%fracdz(k,1)*Thickness%rho_dzt(i,j,k,tau)
        enddo
     enddo
  enddo

  ! time dependent inverse dzwt 
  do k=0,nk
     do j=jsd,jed
        do i=isd,ied
           dzwtr(i,j,k) = 1.0/Thickness%dzwt(i,j,k) 
        enddo
     enddo
  enddo

  call compute_bldepth(Time, Thickness, Dens, T_prog, surf_blthick)
  call tracer_derivs(taum1, T_prog) 
  call compute_psi(Time, Dens, Thickness)

  if(submeso_diag_advect_transport .or. submeso_advective_flux) then 
     call compute_transport(Time, Dens, Thickness)
  endif 

  ! compute tracer flux components and their convergence
  if(submeso_skew_flux) then 
      call compute_submeso_skewsion(Thickness, Dens, Time, T_prog)
  else 
      call compute_submeso_advection(Thickness, Dens, Time, T_prog)
  endif


  ! diagnostics 
  if(id_neutral_rho_submeso > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*                          &
                    (Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k)  &
                    +Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k))
            enddo
         enddo
      enddo
      used = send_data (id_neutral_rho_submeso, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),        &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_wdian_rho_submeso > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               temporary   = Dens%drhodz_zt(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
               wrk1(i,j,k) = Grd%tmask(i,j,k)*                              &
                        (Dens%drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k)  &
                        +Dens%drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)) &
                        /(epsln+temporary)
            enddo
         enddo
      enddo
      used = send_data (id_wdian_rho_submeso, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),      &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


end subroutine submeso_restrat
! </SUBROUTINE> NAME="submeso_restrat"



!#######################################################################
! <SUBROUTINE NAME="compute_bldepth">
!
! <DESCRIPTION>
! Compute the boundary layer depth and kblt.
! </DESCRIPTION>
!
subroutine compute_bldepth(Time, Thickness, Dens, T_prog, surf_blthick) 

  type(ocean_time_type),       intent(in) :: Time
  type(ocean_thickness_type),  intent(in) :: Thickness
  type(ocean_density_type),    intent(in) :: Dens
  type(ocean_prog_tracer_type),intent(in) :: T_prog(:)
  real, dimension(isd:,jsd:),  intent(in) :: surf_blthick

  integer :: i, j, k, tau
  real    :: active_cells 
  real    :: mld_thickness

  tau  = Time%tau 
  hblt = 0.0

  if(use_hblt_constant) then 

      do j=jsc,jec
         do i=isc,iec
            hblt(i,j) = Grd%tmask(i,j,1)*min(constant_hblt, Grd%ht(i,j))
         enddo
      enddo

  elseif(use_hblt_equal_mld) then 

      call calc_mixed_layer_depth(Thickness,                &
           T_prog(index_salt)%field(isd:ied,jsd:jed,:,tau), &
           T_prog(index_temp)%field(isd:ied,jsd:jed,:,tau), &
           Dens%rho(isd:ied,jsd:jed,:,tau),                 &
           Dens%pressure_at_depth(isd:ied,jsd:jed,:),       &
           Time%model_time, hblt, smooth_hblt)

      do j=jsc,jec
         do i=isc,iec
            hblt(i,j) = Grd%tmask(i,j,1)*min(hblt(i,j), Grd%ht(i,j))
         enddo
      enddo

  else 

      do j=jsc,jec
         do i=isc,iec
            hblt(i,j) = Grd%tmask(i,j,1)*min(surf_blthick(i,j), Grd%ht(i,j))
         enddo
      enddo

  endif


  ! set floor to hblt 
  do j=jsc,jec
     do i=isc,iec
        hblt(i,j) = Grd%tmask(i,j,1)*max(minimum_hblt, hblt(i,j))
     enddo
  enddo

  ! halo values needed 
  call mpp_update_domains(hblt(:,:), Dom%domain2d) 


  ! k-index at bottom of hblt
  ! also move the hblt to equal depth_zwt(kblt); this is 
  ! necessary to get a full extent of the streamfunction 
  ! contained in the boundary layer. 
  kblt=0
  do j=jsd,jed
     do i=isd,ied
        if(Grd%kmt(i,j) > 1) then  
            kloop: do k=1,nk
               if(Thickness%depth_zwt(i,j,k) >= hblt(i,j)) then
                   kblt(i,j) = k
                   hblt(i,j) = Thickness%depth_zwt(i,j,k)
                   exit kloop 
               endif
            enddo kloop
        endif
     enddo
  enddo


  ! inverse front length = f/(<N> H), with H=hblt and <N> ave buoyancy freq over hblt
  if(front_length_deform_radius) then 
      do j=jsd,jed
         do i=isd,ied
            buoy_freq_ave(i,j)    = 0.0
            mld_thickness         = epsln
            front_length_inv(i,j) = front_length_const_inv 
            if(kblt(i,j) >= min_kblt) then  
                buoy_freq_ave(i,j)= epsln
                mld_thickness     = epsln
                do k=1,kblt(i,j)
                   mld_thickness      = mld_thickness      + Thickness%dzt(i,j,k)  
                   buoy_freq_ave(i,j) = buoy_freq_ave(i,j) - grav_rho0r*Thickness%dzt(i,j,k)*Dens%drhodz_zt(i,j,k)  
                enddo
                buoy_freq_ave(i,j)    = sqrt(abs(buoy_freq_ave(i,j))/mld_thickness) 
                front_length_inv(i,j) = min(front_length_const_inv, coriolis_param(i,j)/(epsln+mld_thickness*buoy_freq_ave(i,j))) 
            endif 
         enddo
      enddo
  endif


  ! diagnostics  
  if (id_front_length_submeso > 0) then 
      do j=jsd,jed
         do i=isd,ied
            wrk1_2d(i,j) = min(front_length_max, 1.0/(front_length_inv(i,j)+epsln))
         enddo
      enddo
      used = send_data (id_front_length_submeso, wrk1_2d(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),         &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if (id_buoy_freq_ave_submeso > 0) then 
      used = send_data (id_buoy_freq_ave_submeso, buoy_freq_ave(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),                &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if (id_hblt_submeso > 0) then 
      used = send_data (id_hblt_submeso, hblt(:,:),    &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if (id_kblt_submeso > 0) then 
      do j=jsd,jed
         do i=isd,ied
            wrk1_2d(i,j) = kblt(i,j)
         enddo
      enddo
      used = send_data (id_kblt_submeso, wrk1_2d(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


end subroutine compute_bldepth
! </SUBROUTINE> NAME="compute_bldepth"


!#######################################################################
! <SUBROUTINE NAME="tracer_derivs">
!
! <DESCRIPTION>
! Compute the tracer derivatives: 
! horizontal derivative (constant k-level)
! and vertical derivative. 
! </DESCRIPTION>
!
subroutine tracer_derivs(taum1, T_prog)

  integer,                      intent(in) :: taum1
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  integer :: i, j, k, n
  integer :: kp1, kbot
  real    :: tmaski, tmaskj

  do n=1,num_prog_tracers

     dTdx(n)%field(:,:,:) = 0.0
     dTdy(n)%field(:,:,:) = 0.0
     dTdz(n)%field(:,:,:) = 0.0

     do j=jsc-1,jec
        do i=isc-1,iec

           if(kblt(i,j) >= min_kblt) then  
               do k=1,kblt(i,j)
                  tmaski = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
                  tmaskj = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                  dTdx(n)%field(i,j,k) = (T_prog(n)%field(i+1,j,k,taum1)-T_prog(n)%field(i,j,k,taum1)) &
                                          *Grd%dxter(i,j)*tmaski
                  dTdy(n)%field(i,j,k) = (T_prog(n)%field(i,j+1,k,taum1)-T_prog(n)%field(i,j,k,taum1)) &
                                          *Grd%dytnr(i,j)*tmaskj
               enddo
           endif

        enddo
     enddo

     do j=jsd,jed
        do i=isd,ied
           if(kblt(i,j) >= min_kblt) then  
               do k=1,kblt(i,j)
                  kp1 = min(k+1,nk)
                  dTdz(n)%field(i,j,k) = (T_prog(n)%field(i,j,k,taum1)-T_prog(n)%field(i,j,kp1,taum1)) &
                                          *Grd%tmask(i,j,kp1)*dzwtr(i,j,k)
               enddo
           endif
        enddo
     enddo

  enddo

end subroutine tracer_derivs
! </SUBROUTINE> NAME="tracer_derivs"


!#######################################################################
! <SUBROUTINE NAME="compute_psi">
!
! <DESCRIPTION>
! Compute the vector streamfunction 
!
! Units of psi are m^2/sec
!
! If computing skewsion tendency, then need psi at depth_zt. 
! If computing advection tendency, then need psi at depth_zwt. 
!
! </DESCRIPTION>
!
subroutine compute_psi(Time, Dens, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness
  real    :: gradx, grady
  real    :: tmaski, tmaskj
  real    :: gradxrho(0:1), gradyrho(0:1)
  real    :: factor, coefficient
  real    :: active_cells 
  real    :: max_psi
  real    :: hblt_r
  integer :: i, j, k
  integer :: ip, jq
  integer :: tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau = Time%tau

  wrk1 = 0.0   ! depth 
  wrk2 = 0.0
  wrk3 = 0.0   ! mu_submeso

  ! vector streamfunction components
  psix = 0.0
  psiy = 0.0


  ! for holding the depths 
  if(submeso_skew_flux) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = Thickness%depth_zt(i,j,k) 
            enddo
         enddo
      enddo
  else
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = Thickness%depth_zwt(i,j,k) 
            enddo
         enddo
      enddo
  endif

  do j=jsd,jec
     do i=isd,iec

        if(kblt(i,j) >= min_kblt) then  

            gradxrho(:) = 0.0
            gradyrho(:) = 0.0

            do k=1,kblt(i,j)
               do ip=0,1 
                  jq=ip
                  gradx        =   Dens%drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k) &
                                  +Dens%drhodS(i+ip,j,k)*dTdx(index_salt)%field(i,j,k) 
                  grady        =   Dens%drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k) &
                                  +Dens%drhodS(i,j+jq,k)*dTdy(index_salt)%field(i,j,k)
                  gradxrho(ip) = gradxrho(ip) + gradx*Thickness%dzt(i,j,k)  
                  gradyrho(jq) = gradyrho(jq) + grady*Thickness%dzt(i,j,k)  
               enddo

            enddo   ! enddo for k=1,kblt(i,j)

            ! normalization by depth of boundary layer 
            do ip=0,1 
               jq=ip
               gradxrho(ip) = gradxrho(ip)/hblt(i,j)
               gradyrho(jq) = gradyrho(jq)/hblt(i,j)
            enddo

            ! coefficient has units m^6/(sec*kg)
            coefficient = time_factor(i,j)*ce_grav_rho0r*grid_length(i,j)*front_length_inv(i,j)*hblt(i,j)**2
            hblt_r      = 1.0/hblt(i,j)

            ! compute the vector streamfunction (m^2/sec)
            do k=1,kblt(i,j) 

               tmaski = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
               tmaskj = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)

               factor      = (1.0 - 2.0*wrk1(i,j,k)*hblt_r)**2
               wrk3(i,j,k) = (1.0 - factor)*(1.0 + fiveover21*factor) 

               do ip=0,1 
                  jq=ip
                  psix(i,j,k,jq) = -coefficient*wrk3(i,j,k)*gradyrho(jq)*tmaskj
                  psiy(i,j,k,ip) =  coefficient*wrk3(i,j,k)*gradxrho(ip)*tmaski
               enddo

            enddo

        endif ! endif for kblt(i,j) >= min_kblt

     enddo
  enddo


  ! limit magnitude of psi to reduce potential for crashes 
  if(limit_psi) then 
      do j=jsd,jec
         do i=isd,iec
            do k=1,kblt(i,j) 
               max_psi = limit_psi_velocity_scale*Thickness%dzt(i,j,k)
               do ip=0,1 
                  jq=ip
                    psix(i,j,k,jq) = sign(1.0,psix(i,j,k,jq))*min(max_psi,abs(psix(i,j,k,jq))) 
                    psiy(i,j,k,ip) = sign(1.0,psiy(i,j,k,ip))*min(max_psi,abs(psiy(i,j,k,ip))) 
               enddo
            enddo
         enddo
      enddo
  endif


  ! spatially smooth psi
  if(smooth_psi) then 

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

      do ip=0,1
         jq=ip
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec

                  if(Grd%tmask(i,j,k)==1.0) then 

                      active_cells = 4.0       +&
                           Grd%tmask(i-1,j,k)  +&
                           Grd%tmask(i+1,j,k)  +&
                           Grd%tmask(i,j-1,k)  +&
                           Grd%tmask(i,j+1,k)

                      if (active_cells > 4.0) then
                          wrk1(i,j,k) = & 
                               (4.0*psix(i,j,k,jq) +&
                               psix(i-1,j,k,jq)    +&
                               psix(i+1,j,k,jq)    +&
                               psix(i,j-1,k,jq)    +&
                               psix(i,j+1,k,jq)) / active_cells
                          wrk2(i,j,k) =  &
                               (4.0*psiy(i,j,k,ip) +&
                               psiy(i-1,j,k,ip)    +&
                               psiy(i+1,j,k,ip)    +&
                               psiy(i,j-1,k,ip)    +&
                               psiy(i,j+1,k,ip)) / active_cells

                      else
                          wrk1(i,j,k) = psix(i,j,k,jq)
                          wrk2(i,j,k) = psiy(i,j,k,ip)
                      endif

                  endif

               enddo
            enddo

            do j=jsc,jec
               do i=isc,iec
                  psix(i,j,k,jq) = wrk1(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                  psiy(i,j,k,ip) = wrk2(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
               enddo
            enddo

         enddo
      enddo

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

  endif


  ! send diagnostics 
  if (id_psix_submeso > 0) then 
      used = send_data (id_psix_submeso, onehalf*(psix(:,:,:,0)+psix(:,:,:,1)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                           &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_psiy_submeso > 0) then 
      used = send_data (id_psiy_submeso, onehalf*(psiy(:,:,:,0)+psiy(:,:,:,1)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                           &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_mu_submeso > 0) then 
      used = send_data (id_mu_submeso, wrk3(:,:,:),   &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 


  if(debug_this_module) then 
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_submeso_mod: chksums'
      call write_timestamp(Time%model_time)

      wrk1(isc:iec,jsc:jec,:) = psix(isc:iec,jsc:jec,:,0)*Grd%tmask(isc:iec,jsc:jec,:)
      write(stdoutunit,*) 'chksum for psix(:,:,:,0)  = ',  mpp_chksum(wrk1(isc:iec,jsc:jec,:))
      wrk1(isc:iec,jsc:jec,:) = psix(isc:iec,jsc:jec,:,1)*Grd%tmask(isc:iec,jsc:jec,:)
      write(stdoutunit,*) 'chksum for psix(:,:,:,1)  = ',  mpp_chksum(wrk1(isc:iec,jsc:jec,:))
      wrk1(isc:iec,jsc:jec,:) = psiy(isc:iec,jsc:jec,:,0)*Grd%tmask(isc:iec,jsc:jec,:)
      write(stdoutunit,*) 'chksum for psiy(:,:,:,0)  = ',  mpp_chksum(wrk1(isc:iec,jsc:jec,:))
      wrk1(isc:iec,jsc:jec,:) = psiy(isc:iec,jsc:jec,:,1)*Grd%tmask(isc:iec,jsc:jec,:)
      write(stdoutunit,*) 'chksum for psiy(:,:,:,1)  = ',  mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  endif


end subroutine compute_psi
! </SUBROUTINE> NAME="compute_psi"


!#######################################################################
! <SUBROUTINE NAME="compute_transport">
!
! <DESCRIPTION>
! Compute the mass transport from submeso
! </DESCRIPTION>
!
subroutine compute_transport(Time, Dens, Thickness)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness

  real    :: gradx, grady
  real    :: tmaski, tmaskj
  real    :: gradxrho(0:1), gradyrho(0:1)
  real    :: factor, coefficient
  real    :: active_cells 
  real    :: max_psi
  real    :: hblt_r
  integer :: i, j, k
  integer :: ip, jq
  integer :: tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau = Time%tau

  ! store density according to vertical coordinate class 
  wrk1(:,:,:) = 0.0
  if(vert_coordinate_class==DEPTH_BASED) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = rho0
            enddo
         enddo
      enddo
  else 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk1(i,j,k) = Dens%rho(i,j,k,tau)
            enddo
         enddo
      enddo
  endif

  ! compute horiz transport components by taking derivatives of the streamfunction 
  uhrho_et_submeso = 0.0
  vhrho_nt_submeso = 0.0

  ! assume psix and psiy are centred on bottom face of tracer cells.
  ! take thickness to be that of tracer cell, so rho_dzt*/dzt = rho.
  ! psix = psiy = 0 at k=0.
  do k=1,1
     do j=jsd,jec
        do i=isd,iec
           uhrho_et_submeso(i,j,k) =  0.5*(psiy(i,j,k,0)+psiy(i,j,k,1))*wrk1(i,j,k)
           vhrho_nt_submeso(i,j,k) = -0.5*(psix(i,j,k,0)+psix(i,j,k,1))*wrk1(i,j,k)
        enddo
     enddo
  enddo
  do k=2,nk
     do j=jsd,jec
        do i=isd,iec
           uhrho_et_submeso(i,j,k) =                     &
                -0.5*( (psiy(i,j,k-1,0)+psiy(i,j,k-1,1)) &
                      -(psiy(i,j,k,0)  +psiy(i,j,k,1)) ) &
                *wrk1(i,j,k)
           vhrho_nt_submeso(i,j,k) =                    &
                0.5*( (psix(i,j,k-1,0)+psix(i,j,k-1,1)) &
                     -(psix(i,j,k,0)  +psix(i,j,k,1)) ) &
                *wrk1(i,j,k)
        enddo
     enddo
  enddo

  ! compute vertical component by taking divergence of horizontal 
  wrho_bt_submeso = 0.0  
  do k=1,nk
     wrho_bt_submeso(:,:,k) = Grd%tmask(:,:,k)*(wrho_bt_submeso(:,:,k-1) &
       + BDX_ET(uhrho_et_submeso(:,:,k)) + BDY_NT(vhrho_nt_submeso(:,:,k)))    
  enddo

  ! send diagnostics 
  if (id_u_et_submeso > 0) then
      wrk2 = 0.0 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk2(i,j,k) = Grd%tmask(i,j,k)*uhrho_et_submeso(i,j,k) &
                             /(epsln+Thickness%rho_dzt(i,j,k,tau))
            enddo
         enddo
      enddo
      used = send_data (id_u_et_submeso, wrk2(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),   &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_v_nt_submeso > 0) then
      wrk2 = 0.0 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk2(i,j,k) = Grd%tmask(i,j,k)*vhrho_nt_submeso(i,j,k) &
                             /(epsln+Thickness%rho_dzt(i,j,k,tau))
            enddo
         enddo
      enddo
      used = send_data (id_v_nt_submeso, wrk2(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),   &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_w_bt_submeso > 0) then
      wrk2 = 0.0 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk2(i,j,k) = Grd%tmask(i,j,k)*wrho_bt_submeso(i,j,k) &
                             /(epsln+wrk1(i,j,k))
            enddo
         enddo
      enddo
      used = send_data (id_w_bt_submeso, wrk2(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),   &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_uhrho_et_submeso > 0) then 
      used = send_data (id_uhrho_et_submeso, uhrho_et_submeso(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                 &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_vhrho_nt_submeso > 0) then 
      used = send_data (id_vhrho_nt_submeso, vhrho_nt_submeso(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                 &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_wrho_bt_submeso > 0) then 
      used = send_data (id_wrho_bt_submeso, wrho_bt_submeso(:,:,1:nk), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                  &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 

  if (diag_step > 0) then
      if (mod(Time%itt,diag_step) == 0) then
          call maximum_bottom_w_submeso(wrho_bt_submeso(:,:,1:nk),Grd%xt,Grd%yt,Thickness%depth_zwt,rho0r,Grd%kmt)
      endif
  endif

  if(debug_this_module) then 
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_submeso_mod: chksums for transport'
      call write_timestamp(Time%model_time)

      wrk1(isc:iec,jsc:jec,:) = uhrho_et_submeso(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
      write(stdoutunit,*) 'uhrho_et_submeso = ',  mpp_chksum(wrk1(isc:iec,jsc:jec,:))
      wrk1(isc:iec,jsc:jec,:) = vhrho_nt_submeso(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
      write(stdoutunit,*) 'vhrho_nt_submeso = ',  mpp_chksum(wrk1(isc:iec,jsc:jec,:))
      wrk1(isc:iec,jsc:jec,:) = wrho_bt_submeso(isc:iec,jsc:jec,1:nk)*Grd%tmask(isc:iec,jsc:jec,:)
      write(stdoutunit,*) 'wrho_bt_submeso  = ',  mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  endif

end subroutine compute_transport
! </SUBROUTINE> NAME="compute_transport"


!#######################################################################
! <SUBROUTINE NAME="compute_submeso_skewsion">
!
! <DESCRIPTION>
!
! Compute tendency from skewsion for submeso. 
!
! </DESCRIPTION>
!
subroutine compute_submeso_skewsion(Thickness, Dens, Time, T_prog)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer :: i,j,k,n, tau

  tau = Time%tau 
  
  do n=1,num_prog_tracers

     call compute_flux_x(Time,n,T_prog(n))
     call compute_flux_y(Time,n,T_prog(n))
     call compute_flux_z(Time,n,T_prog(n))

     if (Grd%tripolar) then 
         call mpp_update_domains(flux_x(:,:,:), flux_y(:,:,:), Dom_flux_sub%domain2d, &
              gridtype=CGRID_NE) 
     endif

     ! tracer tendency (units rho*dzt * tracer concentration/sec)
     T_prog(n)%wrk1 = 0.0
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%wrk1(i,j,k) = Grd%tmask(i,j,k)*(flux_z(i,j,k-1)-flux_z(i,j,k)           &
                   +(flux_x(i,j,k)-flux_x(i-1,j,k)+flux_y(i,j,k)-flux_y(i,j-1,k))*Grd%datr(i,j) &
                   )
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
           enddo
        enddo
     enddo

     if(id_submeso(n) > 0) then 
         used = send_data (id_submeso(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion, &
              Time%model_time, rmask=Grd%tmask(:,:,:),                              &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

  enddo ! enddo for n=1,num_prog_tracers


  ! diagnose mass transports for sending to diagnostics.
  if (id_tx_trans_submeso > 0 .or. id_ty_trans_submeso > 0) then 

      wrk1_v = 0.0

      if(vert_coordinate_class==DEPTH_BASED) then 
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1_v(i,j,k,1) = &
                    -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*rho0
                   wrk1_v(i,j,k,2) = &
                     transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*rho0
                enddo
             enddo
          enddo
      else
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1_v(i,j,k,1) = &
                    -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*Dens%rho(i,j,k,tau)
                   wrk1_v(i,j,k,2) = &
                     transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*Dens%rho(i,j,k,tau)
                enddo
             enddo
          enddo
      endif

      ! change signs to agree with convention used for ty_trans_gm 
      if (id_tx_trans_submeso > 0) then 
          used = send_data (id_tx_trans_submeso, -1.0*wrk1_v(:,:,:,1), &
               Time%model_time, rmask=Grd%tmask(:,:,:),                &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_ty_trans_submeso > 0) then 
          used = send_data (id_ty_trans_submeso, -1.0*wrk1_v(:,:,:,2), &
               Time%model_time, rmask=Grd%tmask(:,:,:),                &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif

  endif


end subroutine compute_submeso_skewsion
! </SUBROUTINE> NAME="compute_submeso_skewsion"



!#######################################################################
! <SUBROUTINE NAME="compute_flux_x">
!
! <DESCRIPTION>
! Subroutine computes the zonal submesoscale tracer flux component.
!
! fx has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine compute_flux_x(Time,n,Tracer)

  type(ocean_time_type),        intent(in) :: Time
  integer,                      intent(in) :: n
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  integer :: i, j, k
  integer :: ip, kr, kpkr
  real    :: tensor_13(isd:ied,jsd:jed,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)

  flux_x = 0.0

  do k=1,nk

     ! tracer-independent part of the calculation 
     tensor_13(:,:,:) = 0.0
     do ip=0,1 
        do j=jsc,jec
           do i=isc-1,iec
              tensor_13(i,j,ip) = -psiy(i,j,k,ip)
           enddo
        enddo
     enddo

     ! tracer-dependent part of the calculation
     sumz(:,:,:) = 0.0
     do kr=0,1
        do ip=0,1
           do j=jsc,jec
              do i=isc-1,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) &
                      *tensor_13(i,j,ip)*dTdz(n)%field(i+ip,j,k-1+kr)                              &
                      *min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc,jec
        do i=isc-1,iec
           flux_x(i,j,k) = Grd%dxter(i,j)*(sumz(i,j,0)+sumz(i,j,1))
        enddo
     enddo

     if(submeso_limit_flux) then 
         do j=jsc,jec
            do i=isc-1,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then 
                   flux_x(i,j,k) = 0.0
               endif
            enddo
         enddo
     endif

  enddo   ! k-loop


  ! send fluxes to diag_manager 
  if(id_xflux_submeso(n) > 0) then 
      used = send_data (id_xflux_submeso(n), flux_sign*Tracer%conversion*flux_x(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                                     &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  if(id_xflux_submeso_int_z(n) > 0) then 
      wrk1_2d = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + flux_x(i,j,k)
            enddo
         enddo
      enddo
      used = send_data (id_xflux_submeso_int_z(n), flux_sign*Tracer%conversion*wrk1_2d(:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,1),                                          & 
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif



end subroutine compute_flux_x
! </SUBROUTINE> NAME="compute_flux_x"



!#######################################################################
! <SUBROUTINE NAME="compute_flux_y">
!
! <DESCRIPTION>
! Subroutine computes the meridional submesoscale tracer flux component.
!
! fy has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine compute_flux_y(Time,n,Tracer)

  type(ocean_time_type),        intent(in) :: Time
  integer,                      intent(in) :: n
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  integer :: i, j, k
  integer :: jq, kr, kpkr
  real    :: tensor_23(isd:ied,jsd:jed,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)

  flux_y = 0.0

  do k=1,nk

     ! tracer-independent part of the calculation 
     tensor_23(:,:,:) = 0.0
     do jq=0,1  
        do j=jsc-1,jec
           do i=isc,iec
              tensor_23(i,j,jq) = psix(i,j,k,jq)
           enddo
        enddo
     enddo

     ! tracer-dependent part of the calculation
     sumz(:,:,:) = 0.0
     do kr=0,1
        do jq=0,1
           do j=jsc-1,jec
              do i=isc,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) &
                      *tensor_23(i,j,jq)*dTdz(n)%field(i,j+jq,k-1+kr)                              &
                      *min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc-1,jec
        do i=isc,iec
           flux_y(i,j,k) = Grd%dytnr(i,j)*(sumz(i,j,0)+sumz(i,j,1))
        enddo
     enddo

     if(submeso_limit_flux) then 
         do j=jsc-1,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then 
                   flux_y(i,j,k) = 0.0
               endif
            enddo
         enddo
     endif

  enddo   ! k-loop


  ! diagnostics 
  if(id_yflux_submeso(n) > 0) then 
      used = send_data (id_yflux_submeso(n), flux_sign*Tracer%conversion*flux_y(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                                     &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  if(id_yflux_submeso_int_z(n) > 0) then 
      wrk1_2d = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + flux_y(i,j,k)
            enddo
         enddo
      enddo
      used = send_data (id_yflux_submeso_int_z(n), flux_sign*Tracer%conversion*wrk1_2d(:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,1),                                          &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


end subroutine compute_flux_y
! </SUBROUTINE> NAME="compute_flux_y"


!#######################################################################
! <SUBROUTINE NAME="compute_flux_z">
!
! <DESCRIPTION>
! Subroutine computes the vertical submeso tracer flux component.
!
! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
!
! fz has physical dimensions (density*diffusivity*tracer gradient)
!
! </DESCRIPTION>
!
subroutine compute_flux_z(Time,n,Tracer)

  type(ocean_time_type),        intent(in) :: Time
  integer,                      intent(in) :: n
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  integer :: i, j, k, ip, jq, kr
  real :: sumx_0, sumx_1, sumy_0, sumy_1
  real :: temparray31(isc:iec,jsc:jec,0:1,0:1)
  real :: temparray32(isc:iec,jsc:jec,0:1,0:1)
  real :: tensor_31(isd:ied,jsd:jed,0:1)
  real :: tensor_32(isd:ied,jsd:jed,0:1)

  flux_z      = 0.0
  tensor_31   = 0.0
  tensor_32   = 0.0
  temparray31 = 0.0
  temparray32 = 0.0

  do k=1,nk-1

     do ip=0,1
        jq=ip  
        do j=jsc-1,jec
           do i=isc-1,iec
              tensor_31(i,j,ip) =  psiy(i,j,k,ip)
              tensor_32(i,j,jq) = -psix(i,j,k,jq)
           enddo
        enddo
     enddo

     do kr=0,1

        do ip=0,1
           do j=jsc,jec
              do i=isc,iec
                 temparray31(i,j,ip,kr) = tensor_31(i,j,ip)*dtew(i,j,ip) &
                      *min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr))
              enddo
           enddo
        enddo
        do jq=0,1
           do j=jsc,jec
              do i=isc,iec
                 temparray32(i,j,jq,kr) = tensor_32(i,j,jq)*dtns(i,j,jq) &
                      *min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr))  
              enddo
           enddo
        enddo

     enddo

     do j=jsc,jec
        do i=isc,iec
           sumx_0 =  temparray31(i,j,0,0)*dTdx(n)%field(i-1,j,k) &
                  +  temparray31(i,j,0,1)*dTdx(n)%field(i-1,j,k+1)
           sumx_1 =  temparray31(i,j,1,0)*dTdx(n)%field(i,j,k)   &
                  +  temparray31(i,j,1,1)*dTdx(n)%field(i,j,k+1)
           sumy_0 =  temparray32(i,j,0,0)*dTdy(n)%field(i,j-1,k) &
                  +  temparray32(i,j,0,1)*dTdy(n)%field(i,j-1,k+1)
           sumy_1 =  temparray32(i,j,1,0)*dTdy(n)%field(i,j,k)   &
                  +  temparray32(i,j,1,1)*dTdy(n)%field(i,j,k+1)

           flux_z(i,j,k) = Grd%tmask(i,j,k+1)                                      &
                *( Grd%dxtr(i,j)*(sumx_0+sumx_1) + Grd%dytr(i,j)*(sumy_0+sumy_1) ) &
                *dzwtr(i,j,k)

        enddo
     enddo

     if(submeso_limit_flux) then 
         do j=jsc,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then 
                   flux_z(i,j,k) = 0.0
               endif
            enddo
         enddo
     endif


  enddo   ! end of k=1,nk-1 loop


  ! diagnostics 
  if(id_zflux_submeso(n) > 0) then 
      wrk2 = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk2(i,j,k) = Grd%dat(i,j)*flux_z(i,j,k)
            enddo
         enddo
      enddo
      used = send_data (id_zflux_submeso(n), flux_sign*Tracer%conversion*wrk2(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                                   &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


end subroutine compute_flux_z
! </SUBROUTINE> NAME="compute_flux_z"


!#######################################################################
! <SUBROUTINE NAME="compute_submeso_advection">
!
! <DESCRIPTION>
!
! Sweby scheme to compute the tendency from advection for submeso. 
! Presently this scheme has problems, and is not usable.  
!
! </DESCRIPTION>
!
subroutine compute_submeso_advection(Thickness, Dens, Time, T_prog)

  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  real,dimension(isc:iec,jsc:jec) :: ftp
  real,dimension(isc:iec,jsc:jec) :: fbt
  real,dimension(isc:iec,jsc:jec) :: wkm1
  real,dimension(isd:ied,jsd:jed) :: tmp_flux

  integer  :: i, j, k, n
  integer  :: kp1, kp2, km1
  integer  :: tau, taum1
  real     :: Rjm, Rj, Rjp, cfl, massflux
  real     :: d0, d1, thetaP, psiP 
  real     :: thetaM, psiM

  tau     = Time%tau
  taum1   = Time%taum1
  flux_x  = 0.0
  flux_y  = 0.0
  flux_z  = 0.0
  wrk1    = 0.0

  ! calculate flux at bottom face of the T-cells
  do n=1,num_prog_tracers

     ftp  = 0.0
     fbt  = 0.0
     wkm1 = 0.0

     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              T_prog(n)%wrk1(i,j,k) = 0.0
           enddo
        enddo
     enddo

     do k=1,nk

        do j=jsc,jec
           do i=isc,iec
              tracer_mdfl_all(n)%field(i,j,k) = T_prog(n)%field(i,j,k,taum1)
           enddo
        enddo

        kp1 = min(k+1,nk)
        kp2 = min(k+2,nk)
        km1 = max(k-1,1)

        do j=jsc,jec
           do i=isc,iec

              Rjp = (T_prog(n)%field(i,j,km1,taum1) - T_prog(n)%field(i,j,k,taum1))    &
                   *tmask_mdfl(i,j,km1)*tmask_mdfl(i,j,k)
              Rj  = (T_prog(n)%field(i,j,k,taum1) - T_prog(n)%field(i,j,kp1,taum1))    &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j,kp1)
              Rjm = (T_prog(n)%field(i,j,kp1,taum1) - T_prog(n)%field(i,j,kp2,taum1))  &
                   *tmask_mdfl(i,j,kp1)*tmask_mdfl(i,j,kp2)

              massflux = Grd%dat(i,j) * wrho_bt_submeso(i,j,k)
              cfl = abs(wrho_bt_submeso(i,j,k) * dtime  / Thickness%rho_dzt(i,j,k,tau))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              fbt(i,j) =  0.5 * ( ( massflux + abs(massflux) )         &
                   * ( T_prog(n)%field(i,j,kp1,taum1) + psiP * Rj )    &
                   + ( massflux - abs(massflux) )                      &
                   * ( T_prog(n)%field(i,j, k ,taum1) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)

              wrk1(i,j,k) = Grd%datr(i,j)*( fbt(i,j) - ftp(i,j))     &
                   + T_prog(n)%field(i,j,k,taum1)*(wkm1(i,j) - wrho_bt_submeso(i,j,k)) 

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)  &
                   + wrk1(i,j,k) * dtime/Thickness%rho_dzt(i,j,k,tau) 

              flux_z(i,j,k) = fbt(i,j)
              ftp(i,j)      = fbt(i,j)

           enddo
        enddo

        ! update vertical velocity for next k-level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = wrho_bt_submeso(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                              flags=XUPDATE, complete=T_prog(n)%complete)

  enddo ! end of n-loop for tracers  


  ! calculate flux at the eastern face of the T-cells
  do n=1,num_prog_tracers

     do k=1,nk

        do j=jsc,jec
           do i=isc-1,iec

              Rjp = (tracer_mdfl_all(n)%field(i+2,j,k) - tracer_mdfl_all(n)%field(i+1,j,k))    &
                   *tmask_mdfl(i+2,j,k)*tmask_mdfl(i+1,j,k)
              Rj  = (tracer_mdfl_all(n)%field(i+1,j,k) - tracer_mdfl_all(n)%field( i ,j,k) )   &
                   *tmask_mdfl(i+1,j,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k) - tracer_mdfl_all(n)%field(i-1,j,k))      &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i-1,j,k)

              massflux = Grd%dyte(i,j) * uhrho_et_submeso(i,j,k)
              cfl = abs(uhrho_et_submeso(i,j,k) * dtime * 2.0    &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i+1,j,k,tau)) * Grd%dxte(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_x(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )      &
                   * ( tracer_mdfl_all(n)%field( i ,j,k) + psiP * Rj )   &
                   + ( massflux - abs(massflux) )                        &
                   * ( tracer_mdfl_all(n)%field(i+1,j,k) - psiM * Rj ) ) &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
           enddo
        enddo

        ! update the tracer
        do j=jsc,jec
           do i=isc,iec
              wrk1(i,j,k) = tmask_mdfl(i,j,k) * Grd%datr(i,j)       &
                   * (flux_x(i-1,j,k) - flux_x(i,j,k)               &
                   + T_prog(n)%field(i,j,k,taum1)* (                &
                     Grd%dyte(i,j)   * uhrho_et_submeso( i ,j,k)    &
                   - Grd%dyte(i-1,j) * uhrho_et_submeso(i-1,j,k) ) )

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k) &
                   + wrk1(i,j,k)*dtime/Thickness%rho_dzt(i,j,k,tau)
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, &
                              flags=YUPDATE, complete=T_prog(n)%complete)

  enddo ! end of n-loop for tracers  


  ! calculate flux at the northern face of the T-cells 
  do n=1,num_prog_tracers 

     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = 0.0
        enddo
     enddo

     do k=1,nk

        do j=jsc-1,jec
           do i=isc,iec

              Rjp = (tracer_mdfl_all(n)%field(i,j+2,k) - tracer_mdfl_all(n)%field(i,j+1,k))   &
                   *tmask_mdfl(i,j+2,k)*tmask_mdfl(i,j+1,k)
              Rj  = (tracer_mdfl_all(n)%field(i,j+1,k) - tracer_mdfl_all(n)%field(i,j,k))     &
                   *tmask_mdfl(i,j+1,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k)   - tracer_mdfl_all(n)%field(i,j-1,k))   &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j-1,k)

              massflux = Grd%dxtn(i,j) * vhrho_nt_submeso(i,j,k)
              cfl = abs(vhrho_nt_submeso(i,j,k) * dtime * 2.0              &
                   / ((Thickness%rho_dzt(i,j,k,tau) + Thickness%rho_dzt(i,j+1,k,tau)) * Grd%dytn(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_y(i,j,k) =  0.5 * ( ( massflux + abs(massflux) )       &
                   * ( tracer_mdfl_all(n)%field(i,j,k) + psiP * Rj )      &
                   + ( massflux - abs(massflux) )                         &
                   * ( tracer_mdfl_all(n)%field(i,j+1,k) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)

           enddo
        enddo

        ! calculate the overall tendency and update tracer 
        do j=jsc,jec
           do i=isc,iec

              wrk1(i,j,k) = tmask_mdfl(i,j,k)*Grd%datr(i,j)*(flux_y(i,j-1,k)-flux_y(i,j,k))  &
                            + T_prog(n)%field(i,j,k,taum1) *                                 &
                   (wrho_bt_submeso(i,j,k) - wkm1(i,j)                                       &
                   + Grd%datr(i,j)*                                                          &
                   ( Grd%dyte(i-1,j) * uhrho_et_submeso(i-1,j,k)                             &
                   - Grd%dyte( i ,j) * uhrho_et_submeso( i ,j,k)))

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)          &
                   + wrk1(i,j,k)*dtime/Thickness%rho_dzt(i,j,k,tau)

              T_prog(n)%wrk1(i,j,k) = &
                Thickness%rho_dzt(i,j,k,tau)*(tracer_mdfl_all(n)%field(i,j,k)-T_prog(n)%field(i,j,k,taum1))/dtime  &
                   *tmask_mdfl(i,j,k) 

           enddo
        enddo

        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
          enddo
        enddo

        ! update vertical velocity for next level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = wrho_bt_submeso(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop



     ! send fluxes to diag_manager 
     if(id_xflux_submeso(n) > 0) then 
         used = send_data (id_xflux_submeso(n), flux_sign*T_prog(n)%conversion*flux_x(:,:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,:),                                        &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

     if(id_xflux_submeso_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_x(i,j,k)
               enddo
            enddo
         enddo
         used = send_data (id_xflux_submeso_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1),                                             & 
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif

     if(id_yflux_submeso(n) > 0) then 
         used = send_data (id_yflux_submeso(n), flux_sign*T_prog(n)%conversion*flux_y(:,:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,:),                                        &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

     if(id_yflux_submeso_int_z(n) > 0) then 
         wrk1_2d = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_2d(i,j) = wrk1_2d(i,j) + flux_y(i,j,k)
               enddo
            enddo
         enddo
         used = send_data (id_yflux_submeso_int_z(n), flux_sign*T_prog(n)%conversion*wrk1_2d(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1),                                             &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif

     if(id_zflux_submeso(n) > 0) then 
         wrk2 = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk2(i,j,k) = Grd%dat(i,j)*flux_z(i,j,k)
               enddo
            enddo
         enddo
         used = send_data (id_zflux_submeso(n), flux_sign*T_prog(n)%conversion*wrk2(:,:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,:),                                      &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     
     if(id_submeso(n) > 0) then 
         used = send_data (id_submeso(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion, &
              Time%model_time, rmask=Grd%tmask(:,:,:),                                &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif


  enddo ! end of n-loop
  

  ! diagnose mass transports for sending to diagnostics.
  if (id_tx_trans_submeso > 0 .or. id_ty_trans_submeso > 0) then 

      wrk1_v = 0.0

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = transport_convert*uhrho_et_submeso(i,j,k)*Grd%dyte(i,j)
               wrk1_v(i,j,k,2) = transport_convert*vhrho_nt_submeso(i,j,k)*Grd%dxtn(i,j)
            enddo
         enddo
      enddo

      if (id_tx_trans_submeso > 0) then 
          used = send_data (id_tx_trans_submeso, wrk1_v(:,:,:,1), &
               Time%model_time, rmask=Grd%tmask(:,:,:),           &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
      if (id_ty_trans_submeso > 0) then 
          used = send_data (id_ty_trans_submeso, wrk1_v(:,:,:,2), &
               Time%model_time, rmask=Grd%tmask(:,:,:),           &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif

  endif



end subroutine compute_submeso_advection
! </SUBROUTINE> NAME="compute_submeso_advection"


!#######################################################################
! <SUBROUTINE NAME="maximum_bottom_w_general">
!
! <DESCRIPTION>
! Compute maximum vertical velocity from submeso.
! </DESCRIPTION>
!
subroutine maximum_bottom_w_submeso(wb,x_array,y_array,depth_array,scale,km_array)

  real,    dimension(isd:,jsd:,:) , intent(in) :: wb
  real,    dimension(isd:,jsd:)   , intent(in) :: x_array
  real,    dimension(isd:,jsd:)   , intent(in) :: y_array
  real,    dimension(isd:,jsd:,:) , intent(in) :: depth_array
  real,                             intent(in) :: scale
  integer, dimension(isd:,jsd:)   , intent(in) :: km_array

  real    :: wbot, wbot0, fudge
  integer :: i, j, k, iwbot, jwbot, kwbot

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_submesoscale_mod (maximum_bottom_w_submeso): module needs initialization ')
  endif 

  wbot=0.0; iwbot=isc; jwbot=jsc; kwbot=1
  do j=jsc,jec
    do i=isc,iec
      k = km_array(i,j)
      if (k /= 0 .and. (abs(wb(i,j,k)) > abs(wbot))) then
        wbot  = wb(i,j,k)
        iwbot = i
        jwbot = j
        kwbot = k
      endif
    enddo
  enddo
  wbot = scale*abs(wbot)


  write (stdoutunit,'(//60x,a/)') ' Summary of bottom vertical velocity from submeso:'
  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when tracer is independent of processor
  wbot  = wbot*fudge
  wbot0 = wbot
  wbot  = abs(wbot)
  call mpp_max(wbot)

  if (abs(wbot0) == wbot) then
    wbot = wbot0/fudge
    write (unit,9912) wbot, iwbot+Dom%ioff, jwbot+Dom%joff, kwbot, &
                      x_array(iwbot,jwbot), y_array(iwbot,jwbot), depth_array(iwbot,jwbot,kwbot)
  endif
  
9912  format(/' Maximum submeso bottom velocity (',es10.3,' m/s){error}  at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m)')

end subroutine maximum_bottom_w_submeso
! </SUBROUTINE> NAME="maximum_bottom_wrho_submeso"




end module ocean_submesoscale_mod
