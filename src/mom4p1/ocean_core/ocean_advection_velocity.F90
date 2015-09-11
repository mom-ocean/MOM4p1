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
module ocean_advection_velocity_mod
!  
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
!</CONTACT>
!
!<OVERVIEW>
! Advection velocity components for tracer and momenta transport    
!</OVERVIEW>
!
!<DESCRIPTION>
! The module computes the horizontal and vertical components to the 
! advection velocity on the face of tracer and velocity cells.
!
! The three components are related by continuity. 
!
! All components are density weighted, and the 
! horizontal components are thickness weighted. 
!
! Some diagnostics related to fluid mass transport
! classified according to both depth and density classes 
! are also computed.
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies
! The MOM3 Manual (1999)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies
! Elements of MOM4p1 (2006)
! </REFERENCE>
!
! <NOTE>
! The expressions for the horizontal components for tracer advection
! allow for a proper conversion between pressure work and buoyancy.
! </NOTE>
!
! <NOTE>
! The remapping operators are derived from considerations of linear
! interpolation and volume conservation.  A "remapping error" is
! computed to determine consistency between the tracer and velocity
! grid advection velocities.  This error is roundoff only for cases
! where the horizontal tracer and velocity grids are linearly related,
! as is the case for the spherical coordinate version of mom4.  The
! tripolar version of mom4 does not have tracer and velocity grids
! related linearly, and so the "remapping error" is nontrivial.  The 
! significance of this error is unclear.  No adverse effects have been 
! identified. 
! </NOTE>
!
! <NOTE>
! The vertical velocity components for both the tracer and velocity cells
! are diagnosed via continuity (either volume or mass conservation depending
! on the use of the Boussinesq or non-Boussinesq versions of mom4).
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_advection_velocity_nml">
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.   
!  </DATA> 
!
!  <DATA NAME="inflow_nboundary" TYPE="logical">
!  For adding an inflow velocity from the northern boundary.
!  Default is inflow_nboundary=.false.
!  </DATA> 
!
!  <DATA NAME="read_advection_velocity" TYPE="logical">
!  For reading in a file with specified zonal, meridional,
!  and vertical components to the advection velocity.
!  The file should have velocity at the east face of T-cell,
!  north face, and bottom, just as on a C-grid.  The units 
!  should be m/s for each component.  MOM then multiplies
!  but the appropriate thickness and density factors to 
!  generate transport for use in the model.  
!  Default read_zonal_advection_velocity=.false.
!  </DATA> 
!
!  <DATA NAME="constant_advection_velocity" TYPE="logical">
!  When reading in the advection velocity components, we 
!  may choose to keep them constant in time.  This facilitates
!  idealized tests of tracer advection.
!  Default constant_advection_velocity=.false.
!  </DATA> 
!
!  <DATA NAME="max_advection_velocity" UNITS="meter/sec" TYPE="real">
!  This is a check value used to determine if the time steps will result in 
!  linearly stable advection.  If set to a number < 0, then model will estimate the
!  value as a function of maximum grid size.
!  Note that this time step check is not rigorous, and it depends on the details of
!  the advection scheme.  Nonetheless, it provides some useful warning for setting the 
!  time steps in the model.  
!  </DATA> 
!
!</NAMELIST>

use axis_utils_mod,      only: nearest_index 
use constants_mod,       only: epsln
use diag_manager_mod,    only: send_data, register_diag_field, register_static_field
use fms_mod,             only: write_version_number, FATAL, open_namelist_file, close_file, check_nml_error
use fms_mod,             only: read_data, file_exist
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: mpp_error, mpp_chksum, stdout, stdlog
use mpp_mod,             only: mpp_min, mpp_max, mpp_error, mpp_pe

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: BAX, BAY, BDX_ET, BDY_NT, BDX_EU, BDY_NU
use ocean_operators_mod, only: REMAP_ET_TO_EU, REMAP_NT_TO_NU, REMAP_BT_TO_BU
use ocean_parameters_mod,only: missing_value, rho0, rho0r, DEPTH_BASED
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,     only: ocean_grid_type, ocean_thickness_type
use ocean_types_mod,     only: ocean_adv_vel_type, ocean_velocity_type
use ocean_types_mod,     only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,     only: ocean_density_type
use ocean_workspace_mod, only: wrk1
use ocean_obc_mod,       only: ocean_obc_adjust_advel
use ocean_util_mod,      only: write_timestamp

implicit none

private

integer :: unit

! for vertical coordinate class
integer :: vert_coordinate_class

! for diagnostics 
integer :: id_wt               =-1
integer :: id_wrho_bt          =-1
integer :: id_uhrho_et         =-1
integer :: id_vhrho_nt         =-1
integer :: id_wrho_bu          =-1
integer :: id_uhrho_eu         =-1
integer :: id_vhrho_nu         =-1
integer :: id_rhodz_vorticity_z=-1
integer :: id_mass_source      =-1
integer :: id_horz_diverge_t   =-1
integer :: id_horz_diverge_u   =-1
integer :: id_courant_uet      =-1
integer :: id_courant_vnt      =-1
integer :: id_courant_wbt      =-1
integer :: id_courant_ueu      =-1
integer :: id_courant_vnu      =-1
integer :: id_courant_wbu      =-1
integer :: id_ue               =-1
integer :: id_vn               =-1
integer :: id_wb               =-1
integer :: id_ue_rhodzt        =-1
integer :: id_vn_rhodzt        =-1
integer :: id_wb_rho           =-1

real    :: dtts
real    :: dtuv

logical :: used

#include <ocean_memory.h>

#ifdef MOM4_STATIC_ARRAYS
 real, dimension(isd:ied,jsd:jed)    :: tmp 
#else
 real, dimension(:,:),   allocatable :: tmp 
#endif

real, dimension(:,:,:), allocatable :: vt_inflow  ! m/s
real, dimension(:,:,:), allocatable :: mask_inflow 
real, dimension(:,:,:), allocatable :: ue !m/s from file 
real, dimension(:,:,:), allocatable :: vn !m/s from file 
real, dimension(:,:,:), allocatable :: wb !m/s from file

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_advection_velocity_init
public ocean_advection_velocity

character(len=128) :: version=&
     '$Id: ocean_advection_velocity.F90,v 16.0.2.1.84.1 2009/10/10 00:41:49 nnz Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: have_obc              = .false.
logical :: module_is_initialized = .FALSE.

real    :: max_advection_velocity       = -1.0 
logical :: debug_this_module            = .false.
logical :: inflow_nboundary             = .false. 
logical :: read_advection_velocity      = .false.
logical :: constant_advection_velocity  = .false. 

namelist /ocean_advection_velocity_nml/ debug_this_module, max_advection_velocity, inflow_nboundary,  &
                                        read_advection_velocity, constant_advection_velocity

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_advection_velocity_init">
!
! <DESCRIPTION>
! Initialize the advection velocity module
! </DESCRIPTION>
!
subroutine ocean_advection_velocity_init(Grid, Domain, Time, Time_steps, Thickness, Adv_vel, &
                                         ver_coordinate_class, obc, debug)

  type(ocean_grid_type),       intent(in), target   :: Grid
  type(ocean_domain_type),     intent(in), target   :: Domain  
  type(ocean_time_type),       intent(in), target   :: Time
  type(ocean_time_steps_type), intent(in)           :: Time_steps
  type(ocean_thickness_type),  intent(in)           :: Thickness
  type(ocean_adv_vel_type),    intent(inout)        :: Adv_vel
  integer,                     intent(in)           :: ver_coordinate_class
  logical,                     intent(in)           :: obc
  logical,                     intent(in), optional :: debug

  integer :: ioun, io_status, ierr
  real    :: gridmin, gridmax, dtadv, a, b
  integer :: i, j, icg, jcg
  real    :: vertical_factor = 2.0
  real    :: max_dt_for_advection
  real    :: max_dt_for_advection0

  real, dimension(7) :: res = (/0.0625, 0.125, 0.25, 0.50, 1.00, 2.0, 4.0/) ! max resolution (deg)
  real, dimension(7) :: vel = (/2.0,    1.8,   1.6,  1.45, 1.25, 0.9, 0.4/) ! estimated max advection velocity (m/sec)

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
     call mpp_error(FATAL,&
    '==>Error from ocean_advection_velocity_mod: module already initialized.')
  endif 

  call write_version_number( version, tagname )

  module_is_initialized = .TRUE.

  have_obc = obc
  vert_coordinate_class = ver_coordinate_class

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file ()
  read  (ioun, ocean_advection_velocity_nml,iostat=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_advection_velocity_nml)
  write (stdlogunit, ocean_advection_velocity_nml)
  ierr = check_nml_error(io_status,'ocean_advection_velocity_nml')
  call close_file(ioun)

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_advection_velocity_mod with debug_this_module=.true.'  
  endif 

  Dom => Domain
  Grd => Grid

  dtts = Time_steps%dtts
  dtuv = Time_steps%dtuv
  
#ifndef MOM4_STATIC_ARRAYS
  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grd%nk
 
  ! allocate variables
  allocate ( Adv_vel%uhrho_et(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%vhrho_nt(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%uhrho_eu(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%vhrho_nu(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%wrho_bt(isd:ied,jsd:jed,0:nk) )
  allocate ( Adv_vel%wrho_bu(isd:ied,jsd:jed,0:nk) )
  allocate ( Adv_vel%diverge_t(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%diverge_u(isd:ied,jsd:jed,nk) )
  allocate ( tmp(isd:ied,jsd:jed) )

#endif

  Adv_vel%uhrho_et(:,:,:)  = 0.0
  Adv_vel%vhrho_nt(:,:,:)  = 0.0
  Adv_vel%wrho_bt(:,:,:)   = 0.0
  Adv_vel%uhrho_eu(:,:,:)  = 0.0
  Adv_vel%vhrho_nu(:,:,:)  = 0.0
  Adv_vel%wrho_bu(:,:,:)   = 0.0
  Adv_vel%diverge_t(:,:,:) = 0.0
  Adv_vel%diverge_u(:,:,:) = 0.0
  tmp(:,:)                 = 0.0

  if(inflow_nboundary) then 
    call inflow_nboundary_init
  endif 

  if(read_advection_velocity) then 
    call read_advect_velocity(Time, Thickness, Adv_vel)
  endif 
  if(constant_advection_velocity) then 
    write(stdoutunit,'(a)') &
    '==>Note: Using ocean_advection_velocity_mod with constant_advection_velocity=.true.'  
    write(stdoutunit,'(a)') &
    '         Will hold the advection velocity components constant in time. '
  endif 


  ! estimate the max advective speed for the given resolution
  if (max_advection_velocity < 0.0) then
    gridmax = 0.0
    do j=jsc,jec
      do i=isc,iec
        if (Grd%kmt(i,j) /= 0) then
          gridmax = max(gridmax,Grd%dyt(i,j),Grd%dxt(i,j))
        endif
      enddo
    enddo
    call mpp_max (gridmax)
    gridmax = gridmax/111.324e3 ! convert to degrees from metres
    i = nearest_index(gridmax, res)
    if (gridmax <= res(1)) then
      max_advection_velocity =  vel(1)
    elseif (gridmax >= res(7)) then
      max_advection_velocity =  vel(7)
    elseif (gridmax < res(i)) then
      a = res(i)-gridmax
      b = gridmax-res(i-1)
      max_advection_velocity = (a*vel(i-1) + b*vel(i))/(a+b)      
    elseif (gridmax > res(i)) then
      a = gridmax-res(i)
      b = res(i+1)-gridmax
      max_advection_velocity = (a*vel(i+1) + b*vel(i))/(a+b)      
    endif
    write (stdoutunit,'(/a,f5.2,a/a/)') &
    ' Note: The max_advection_velocity is estimated to be ',max_advection_velocity,&
    ' m/s based on grid resolution.','       A more appropriate value can be specified via namelist.'
  endif

  icg = isc; jcg = jsc; gridmin = 1.0e20; max_dt_for_advection = 1.e6; gridmax = 0.0
  do j=jsc,jec
    do i=isc,iec
      if (Grd%kmt(i,j) /= 0) then
        gridmin = min(gridmin,Grd%dyt(i,j),Grd%dxt(i,j))
        gridmax = max(gridmax,Grd%dyt(i,j),Grd%dxt(i,j))
        dtadv = 0.5*gridmin/max_advection_velocity
        if (dtadv < max_dt_for_advection) then
          max_dt_for_advection = dtadv; icg  = i; jcg  = j
        endif
      endif
    enddo
  enddo
  
  max_dt_for_advection = nint(max_dt_for_advection/vertical_factor) ! Questions? -> Ronald.Pacanowski@noaa.gov
  max_dt_for_advection =  max_dt_for_advection + 0.001*mpp_pe()     ! to separate redundancies
  max_dt_for_advection0 = max_dt_for_advection
  call mpp_min (max_dt_for_advection)
  call mpp_max (gridmax)
  
  if (max_dt_for_advection == max_dt_for_advection0) then
    if (max(dtts,dtuv) <= max_dt_for_advection) then
      write (unit,'(/a,i4,a,i4,a,f8.3,a,f8.3,a)') &
       ' Note: Advection stability most nearly violated at T-cell (i,j) = (',&
        icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xt(icg,jcg),',',Grd%yt(icg,jcg),').'
    else
      write (unit,'(/a,i4,a,i4,a,f8.3,a,f8.3,a)') &
      '=>Error: Advection stability violated at T-cell (i,j) = (',&
      icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xt(icg,jcg),',',Grd%yt(icg,jcg),').'
      call mpp_error(FATAL,&
      '==>Error from ocean_advection_velocity_mod: tracer advection stability is violated.')
    endif
    write (unit,'(a,f7.2,a,/a,f12.2,a,/a,f12.2,a,f12.2,a)')                                  &
     ' Assuming a maximum advection velocity of ',max_advection_velocity,' m/s,',            &
     ' Linear stability requires max(dtuv,dtts) be less than ',max_dt_for_advection,' sec.', &
     ' Model is now using (dtuv,dtts) = (',dtuv,',',dtts,') sec.'
  endif
  max_dt_for_advection = nint(max_dt_for_advection)

  ! register advective velocity components for diagnostic output

  id_horz_diverge_t = register_diag_field ('ocean_model', 'horz_diverge_t',                 &
     Grd%vel_axes_wt(1:3), Time%model_time,                                                 &
    'horizontal (on k=constant) divergence of (uhrho,vhrho) at T-points', '(kg/m^3)*m/sec', &
    missing_value=missing_value, range=(/-10.e4,10.e4/))

  id_wt = register_diag_field ('ocean_model', 'wt', Grd%vel_axes_wt(1:3), Time%model_time, &
    'dia-surface velocity T-points', 'm/sec', missing_value=missing_value, range=(/-10.e4,10.e4/))

  id_wrho_bt = register_diag_field ('ocean_model', 'wrhot', Grd%vel_axes_wt(1:3), Time%model_time, &
    'rho*dia-surface velocity T-points', '(kg/m^3)*m/sec', missing_value=missing_value, range=(/-10.e4,10.e4/))
  id_uhrho_et = register_diag_field ('ocean_model', 'uhrho_et', Grd%tracer_axes_flux_x(1:3), Time%model_time, &
    'uhrho_et on horz face of T-cells', '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-100.e4,100.e4/))
  id_vhrho_nt = register_diag_field ('ocean_model', 'vhrho_nt', Grd%tracer_axes_flux_y(1:3), Time%model_time, &
    'vhrho_nt on horz face of T-cells', '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-100.e4,100.e4/))

  id_horz_diverge_u = register_diag_field ('ocean_model', 'horz_diverge_u',                 &
     Grd%vel_axes_wu(1:3), Time%model_time,                                                 &
    'horizontal (on k=constant) divergence of (uhrho,vhrho) at U-points', '(kg/m^3)*m/sec', &
    missing_value=missing_value, range=(/-10.e4,10.e4/))

  id_wrho_bu = register_diag_field ('ocean_model', 'wrhou', Grd%vel_axes_wu(1:3), Time%model_time, &
    'rho*dia-surface velocity U-points', '(kg/m^3)*m/sec', missing_value=missing_value, range=(/-10.e4,10.e4/))
  id_uhrho_eu = register_diag_field ('ocean_model', 'uhrho_eu', Grd%tracer_axes(1:3), Time%model_time, &
    'uhrho_eu on horz face of U-cells', '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-100.e4,100.e4/))
  id_vhrho_nu = register_diag_field ('ocean_model', 'vhrho_nu', Grd%tracer_axes(1:3), Time%model_time, &
    'vhrho_nu on horz face of U-cells', '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-100.e4,100.e4/))

  ! mass-based vorticity 
  id_rhodz_vorticity_z = register_diag_field ('ocean_model', 'rhodz_vorticity_z', Grd%vel_axes_uv(1:3), &
     Time%model_time, 'vertical rhodz vorticity: (rho*dz*v)_x-(rho*dz*u)_y', &
     '(kg*m^2)/sec', missing_value=missing_value, range=(/-1e8,1e8/))

  ! register Courant numbers for diagnostic output
  id_courant_uet = register_diag_field ('ocean_model', 'courant_uet', Grd%tracer_axes(1:3), Time%model_time, &
    'Courant number [uet*dtts/dxt]', 'none', missing_value=missing_value, range=(/-1e9,1e9/))
  id_courant_vnt = register_diag_field ('ocean_model', 'courant_vnt', Grd%tracer_axes(1:3), Time%model_time, &
    'Courant number [vnt*dtts/dyt]', 'none', missing_value=missing_value, range=(/-1e9,1e9/))
  id_courant_wbt = register_diag_field ('ocean_model', 'courant_wbt', Grd%tracer_axes(1:3), Time%model_time, &
    'Courant number [wbt*dtts/dzt]', 'none', missing_value=missing_value, range=(/-1e9,1e9/))
  id_courant_ueu = register_diag_field ('ocean_model', 'courant_ueu', Grd%vel_axes_uv(1:3), Time%model_time, &
    'Courant number [ueu*dtuv/dxu]', 'none', missing_value=missing_value, range=(/-1e9,1e9/))
  id_courant_vnu = register_diag_field ('ocean_model', 'courant_vnu', Grd%vel_axes_uv(1:3), Time%model_time, &
    'Courant number [vnu*dtuv/dyu]', 'none', missing_value=missing_value, range=(/-1e9,1e9/))
  id_courant_wbu = register_diag_field ('ocean_model', 'courant_wbu', Grd%vel_axes_uv(1:3), Time%model_time, &
    'Courant number [wnu*dtuv/dzu]', 'none', missing_value=missing_value, range=(/-1e9,1e9/))

  ! register mass source 
  id_mass_source = register_diag_field ('ocean_model', 'mass_source', Grd%tracer_axes(1:3), Time%model_time, &
    'Mass source on T-grid', '(kg/m^3)*(m/s)', missing_value=missing_value, range=(/-1e9,1e9/))


end subroutine ocean_advection_velocity_init
! </SUBROUTINE> NAME="ocean_advection_velocity_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_advection_velocity">
!
! <DESCRIPTION>
! Compute thickness weighted and density weighted advection velocity 
! components for the B-grid on the T-cells and U-cells. 
! </DESCRIPTION>
!
subroutine ocean_advection_velocity (Velocity, Time, Thickness, Dens, pme, river, Adv_vel)
  
  type(ocean_velocity_type),  intent(in)    :: Velocity
  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_density_type),   intent(in)    :: Dens
  real, dimension(isd:,jsd:), intent(in)    :: pme
  real, dimension(isd:,jsd:), intent(in)    :: river
  type(ocean_adv_vel_type),   intent(inout) :: Adv_vel

  integer ::  i, j, k, tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) &
   call mpp_error(FATAL, &
    '==>Error from ocean_advection_velocity_mod (ocean_advection_velocity): module must be initialized')

  tau  = Time%tau
  wrk1 = 0.0
  tmp  = 0.0

  ! compute thickness weighted and density weighted C-grid advection 
  ! velocity components (uhrho_et, vhrho_nt) to advect tracer fields. 
  ! Remap these to get the corresponding advection velocity components 
  ! (uhrho_eu,vhrho_nu) to advect B-grid horizontal momentum.

  if(.not. constant_advection_velocity) then 

      do k=1,nk
         Adv_vel%uhrho_et(:,:,k) = &
         BAY(Velocity%u(:,:,k,1,tau)*Grd%dyu(:,:)*Thickness%rho_dzu(:,:,k,tau))/Grd%dyte(:,:) 
         Adv_vel%vhrho_nt(:,:,k) = &
         BAX(Velocity%u(:,:,k,2,tau)*Grd%dxu(:,:)*Thickness%rho_dzu(:,:,k,tau))/Grd%dxtn(:,:) 
         Adv_vel%uhrho_eu(:,:,k) = REMAP_ET_TO_EU(Adv_vel%uhrho_et(:,:,k))
         Adv_vel%vhrho_nu(:,:,k) = REMAP_NT_TO_NU(Adv_vel%vhrho_nt(:,:,k))
         Adv_vel%diverge_t(:,:,k)= Grd%tmask(:,:,k) &
           *(BDX_ET(Adv_vel%uhrho_et(:,:,k)) + BDY_NT(Adv_vel%vhrho_nt(:,:,k)))    
         Adv_vel%diverge_u(:,:,k)= Grd%umask(:,:,k) &
           *(BDX_EU(Adv_vel%uhrho_eu(:,:,k)) + BDY_NU(Adv_vel%vhrho_nu(:,:,k)))    
      enddo

      ! dia-surface mass flux through ocean surface arises from water crossing ocean surface. 
      ! minus sign arises from the model sign convention.
      Adv_vel%wrho_bt(:,:,0) = -(pme(:,:) + river(:,:))
      Adv_vel%wrho_bu(:,:,0) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Adv_vel%wrho_bt(:,:,0))         

      do k=1,nk
         tmp(:,:) = Thickness%rho_dzt_tendency(:,:,k) - Thickness%mass_source(:,:,k)
         Adv_vel%wrho_bt(:,:,k) = (tmp(:,:) + Adv_vel%diverge_t(:,:,k) + Adv_vel%wrho_bt(:,:,k-1)) &
                                  *Grd%tmask(:,:,k)
         tmp(:,:) = Grd%umask(:,:,k)*REMAP_BT_TO_BU(tmp(:,:))
         Adv_vel%wrho_bu(:,:,k) = tmp(:,:) + Adv_vel%diverge_u(:,:,k) + Adv_vel%wrho_bu(:,:,k-1) 
      enddo

  endif

  if (have_obc) then
     call ocean_obc_adjust_advel(Adv_vel) 
  endif


  ! specify a value for the inflow meridional transport 
  if(inflow_nboundary) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Adv_vel%vhrho_nt(i,j,k) = Adv_vel%vhrho_nt(i,j,k)  &
                + Grd%tmask(i,j,k)*mask_inflow(i,j,k)*rho0*Thickness%dzt(i,j,k)*vt_inflow(i,j,k)
            enddo
         enddo
      enddo
  endif

  
  if(debug_this_module) then 
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_advection_velocity_mod: chksums'
      call write_timestamp(Time%model_time)
      write(stdoutunit,*) 'Adv_vel%uhrho_et   ==> ', mpp_chksum(Adv_vel%uhrho_et(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'Adv_vel%vhrho_nt   ==> ', mpp_chksum(Adv_vel%vhrho_nt(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'Adv_vel%uhrho_eu   ==> ', mpp_chksum(Adv_vel%uhrho_eu(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'Adv_vel%vhrho_nu   ==> ', mpp_chksum(Adv_vel%vhrho_nu(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'Adv_vel%wrho_bt    ==> ', mpp_chksum(Adv_vel%wrho_bt(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'Adv_vel%wrho_bu    ==> ', mpp_chksum(Adv_vel%wrho_bu(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'Adv_vel%wrho_bt(0) ==> ', mpp_chksum(Adv_vel%wrho_bt(isc:iec,jsc:jec,0))
      write(stdoutunit,*) 'Adv_vel%wrho_bu(0) ==> ', mpp_chksum(Adv_vel%wrho_bu(isc:iec,jsc:jec,0))
  endif

  ! send advective velocity components to diagnostic manager 

  if (id_wt > 0)  then
      if(vert_coordinate_class==DEPTH_BASED) then 
          used = send_data (id_wt, rho0r*Adv_vel%wrho_bt(:,:,1:nk), &
                 Time%model_time, rmask=Grd%tmask(:,:,:),           &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      else 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   if(Grd%tmask(i,j,k) > 0.0) then 
                       wrk1(i,j,k) = Adv_vel%wrho_bt(i,j,k)/Dens%rho(i,j,k,tau)
                   endif
                enddo
             enddo
          enddo
          used = send_data (id_wt, wrk1(:,:,:),           &
                 Time%model_time, rmask=Grd%tmask(:,:,:), &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif
  endif

  if (id_uhrho_et > 0) used = send_data (id_uhrho_et, Adv_vel%uhrho_et(:,:,:), &
                           Time%model_time, rmask=Grd%tmask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_vhrho_nt > 0) used = send_data (id_vhrho_nt, Adv_vel%vhrho_nt(:,:,:), &
                           Time%model_time, rmask=Grd%tmask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_wrho_bt > 0) used =  send_data (id_wrho_bt, Adv_vel%wrho_bt(:,:,1:nk), &
                           Time%model_time, rmask=Grd%tmask(:,:,:),&
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_uhrho_eu > 0) used = send_data (id_uhrho_eu, Adv_vel%uhrho_eu(:,:,:), &
                           Time%model_time, rmask=Grd%umask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_vhrho_nu > 0) used = send_data (id_vhrho_nu, Adv_vel%vhrho_nu(:,:,:), &
                           Time%model_time, rmask=Grd%umask(:,:,:),&
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_wrho_bu > 0) used =  send_data (id_wrho_bu, Adv_vel%wrho_bu(:,:,1:nk), &
                           Time%model_time, rmask=Grd%umask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if(id_horz_diverge_t > 0) then 
      used = send_data (id_horz_diverge_t, Adv_vel%diverge_t(:,:,:),  &
                 Time%model_time, rmask=Grd%tmask(:,:,:),             &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_horz_diverge_u > 0) then 
      used = send_data (id_horz_diverge_u, Adv_vel%diverge_u(:,:,:), &
                 Time%model_time, rmask=Grd%umask(:,:,:),            &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  ! send mass and thickness weighted vertical vorticity 
  if(id_rhodz_vorticity_z > 0) then 
      wrk1=0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = (Adv_vel%vhrho_nt(i+1,j,k)-Adv_vel%vhrho_nt(i,j,k))*Grd%dxur(i,j) &
                            -(Adv_vel%uhrho_et(i,j+1,k)-Adv_vel%uhrho_et(i,j,k))*Grd%dyur(i,j) 
            enddo
         enddo
      enddo
      used = send_data (id_rhodz_vorticity_z, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%umask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 


  ! send Courant numbers to diagnostic manager 
  if (id_courant_uet > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtts*Adv_vel%uhrho_et(i,j,k)/(Thickness%rho_dzt(i,j,k,tau)*Grd%dxt(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_uet, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_vnt > 0) then 
      do k=1,nk  
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtts*Adv_vel%vhrho_nt(i,j,k)/(Thickness%rho_dzt(i,j,k,tau)*Grd%dyt(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_vnt, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_wbt > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtts*Adv_vel%wrho_bt(i,j,k)/Thickness%rho_dzt(i,j,k,tau)
            enddo
         enddo
      enddo
      used = send_data (id_courant_wbt, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_ueu > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtuv*Adv_vel%uhrho_eu(i,j,k)/(Thickness%rho_dzu(i,j,k,tau)*Grd%dxu(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_ueu, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%umask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_courant_vnu > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtuv*Adv_vel%vhrho_nu(i,j,k)/(Thickness%rho_dzu(i,j,k,tau)*Grd%dyu(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_vnu, wrk1(:,:,:),  &
             Time%model_time, rmask=Grd%umask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_wbu > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtuv*Adv_vel%wrho_bu(i,j,k)/Thickness%rho_dzu(i,j,k,tau)
            enddo
         enddo
      enddo
      used = send_data (id_courant_wbu, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%umask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


end subroutine ocean_advection_velocity
! </SUBROUTINE> NAME="advection_velocity"



!#######################################################################
! <SUBROUTINE NAME="read_advect_velocity">
!
! <DESCRIPTION>
! For reading in the advection velocity components.  Assume that 
! the advection velocity components read in from a file are in units 
! of meter/sec and placed on the T-cell faces, as in a C-grid ocean model.  
! 
! This routine assumes that the read-in velocity components already 
! have the proper masking.  
!
! The main application of this routine is for developing idealized
! test cases for tracer advection.   
! </DESCRIPTION>
!
subroutine read_advect_velocity(Time, Thickness, Adv_vel)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_adv_vel_type),   intent(inout) :: Adv_vel

  integer            :: k
  character(len=128) :: filename

  allocate (ue(isd:ied,jsd:jed,nk))
  ue(:,:,:) = 0.0
  allocate (vn(isd:ied,jsd:jed,nk))
  vn(:,:,:) = 0.0
  allocate (wb(isd:ied,jsd:jed,nk))
  wb(:,:,:) = 0.0

  filename = 'INPUT/ocean_advect_velocity.nc'
  if (file_exist(trim(filename))) then 
      call read_data(filename,'ue',ue,Dom%domain2d,timelevel=1)
      call mpp_update_domains(ue(:,:,:), Dom%domain2d)
      call read_data(filename,'vn',vn,Dom%domain2d,timelevel=1)
      call mpp_update_domains(vn(:,:,:), Dom%domain2d)
      call read_data(filename,'wb',wb,Dom%domain2d,timelevel=1)
      call mpp_update_domains(wb(:,:,:), Dom%domain2d)
  else 
      call mpp_error(FATAL,&
      '==>Error from ocean_advection_velocity_mod: read_advect_velocity cannot find ocean_advect_velocity.nc.')
  endif

  do k=1,nk
     Adv_vel%uhrho_et(:,:,k) = rho0*ue(:,:,k)*BAY(Thickness%dzu(:,:,k)) 
     Adv_vel%vhrho_nt(:,:,k) = rho0*vn(:,:,k)*BAX(Thickness%dzu(:,:,k)) 
     Adv_vel%wrho_bt(:,:,k)  = rho0*wb(:,:,k)
     Adv_vel%uhrho_eu(:,:,k) = REMAP_ET_TO_EU(Adv_vel%uhrho_et(:,:,k))                                       
     Adv_vel%vhrho_nu(:,:,k) = REMAP_NT_TO_NU(Adv_vel%vhrho_nt(:,:,k))                                       
     Adv_vel%wrho_bu(:,:,k)  = REMAP_BT_TO_BU(Adv_vel%wrho_bt(:,:,k))
  enddo

  id_ue = register_static_field('ocean_model','ue', Grd%tracer_axes_flux_x(1:3), &
          'zonal C-grid velocity from file','m/s', range=(/-1e3,1e3/))
  if (id_ue > 0) used = send_data(id_ue, ue(isc:iec,jsc:jec,:), Time%model_time)

  id_ue_rhodzt = register_static_field('ocean_model','ue_rhodzt', Grd%tracer_axes_flux_x(1:3), &
                 'zonal C-grid velocity*rhodzt from file','(kg/m3)*m^2/s', range=(/-1e6,1e6/))
  if (id_ue_rhodzt > 0) used = send_data(id_ue_rhodzt, Adv_vel%uhrho_et(isc:iec,jsc:jec,:), Time%model_time)


  id_vn = register_static_field('ocean_model','vn', Grd%tracer_axes_flux_y(1:3), &
          'meridional C-grid velocity from file','m/s', range=(/-1e3,1e3/))
  if (id_vn > 0) used = send_data(id_vn, vn(isc:iec,jsc:jec,:), Time%model_time)

  id_vn_rhodzt = register_static_field('ocean_model','vn_rhodzt', Grd%tracer_axes_flux_y(1:3), &
                 'merid C-grid velocity*rhodzt from file','(kg/m3)*m^2/s', range=(/-1e6,1e6/))
  if (id_vn_rhodzt > 0) used = send_data(id_vn_rhodzt, Adv_vel%vhrho_nt(isc:iec,jsc:jec,:), Time%model_time)


  id_wb = register_static_field('ocean_model','wb', Grd%tracer_axes_wt(1:3), &
          'vertical C-grid velocity from file','m/s', range=(/-1e3,1e3/))
  if (id_wb > 0) used = send_data(id_wb, wb(isc:iec,jsc:jec,:), Time%model_time)

  id_wb_rho = register_static_field('ocean_model','wb_rho', Grd%tracer_axes_wt(1:3), &
              'vertical C-grid velocity*rho from file','(kg/m3)*m/s', range=(/-1e6,1e6/))
  if (id_wb_rho > 0) used = send_data(id_wb_rho, Adv_vel%wrho_bt(isc:iec,jsc:jec,1:nk), Time%model_time)



end subroutine read_advect_velocity
! </SUBROUTINE> NAME="read_advect_velocity"




!#######################################################################
! <SUBROUTINE NAME="inflow_nboundary_init">
!
! <DESCRIPTION>
! Initialize the advection velocity used for specifying a nonzero 
! southward inflow introduced to the domain from the northern boundary.  
!
! </DESCRIPTION>
!
subroutine inflow_nboundary_init

  character(len=128) :: filename

  allocate ( vt_inflow(isd:ied,jsd:jed,nk) )
  allocate ( mask_inflow(isd:ied,jsd:jed,nk) )

  vt_inflow(:,:,:)   = 0.0
  mask_inflow(:,:,:) = 0.0


  filename = 'INPUT/mask_inflow.nc'
  if (file_exist(trim(filename))) then 
      call read_data(filename,'mask_inflow',mask_inflow,Dom%domain2d,timelevel=1)
      call mpp_update_domains(mask_inflow(:,:,:), Dom%domain2d)
  else 
     call mpp_error(FATAL,&
    '==>Error from ocean_advection_velocity_mod: inflow_nboundary_init cannot find mask_inflow file.')
  endif 

  filename = 'INPUT/vt_inflow.nc'
  if (file_exist(trim(filename))) then 
      call read_data(filename,'vt_inflow',vt_inflow,Dom%domain2d,timelevel=1)
      call mpp_update_domains(vt_inflow(:,:,:), Dom%domain2d)
  else 
     call mpp_error(FATAL,&
    '==>Error from ocean_advection_velocity_mod: inflow_nboundary_init cannot find vt_inflow.')
  endif 

end subroutine inflow_nboundary_init
! </SUBROUTINE> NAME="inflow_nboundary_init"




end module ocean_advection_velocity_mod
