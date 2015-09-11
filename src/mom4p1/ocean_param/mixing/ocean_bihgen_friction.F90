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
module ocean_bihgen_friction_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module computes thickness weighted and density weighted
! time tendency for horizontal velocity arising from 
! biharmonic friction. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes thickness weighted and density weighted
! time tendency for horizontal velocity arising from biharmonic
! friction. 
!
! The viscosity used to determine the strength of the tendency 
! can be a general function of space and time as specified by 
! the Smagorinsky approach as well as a grid-scale dependent
! background viscosity.  The form of the friction operator 
! can be isotropic or anisotropic in the lateral plane. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies and R.W. Hallberg, 
! Biharmonic friction with a Smagorinsky viscosity for use
! in large-scale eddy-permitting ocean models
! Monthly Weather Review vol 128 (2000) pages 2935--2946
! </REFERENCE>
!
! <REFERENCE>
! R.D. Smith and J.C. McWilliams, 
! Anisotropic viscosity for ocean models
! Ocean Modelling 
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! Fundamentals of Ocean Climate Models 
! Princeton University Press 2004
! </REFERENCE>
!
! <REFERENCE>
!  Elements of mom4p1 (2005), S.M. Griffies 
! </REFERENCE>
!
! <NOTE>
! The ocean model can generally run with both Laplacian and biharmonic
! friction enabled at the same time.  Such has been found useful 
! for some eddying ocean simulations. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_bihgen_friction_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging by printing checksums.  
!  </DATA> 
!
!  <DATA NAME="k_smag_iso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="k_smag_aniso" UNITS="dimensionless" TYPE="real">
!  This is the dimensionless Smagorinsky coefficient used to set the scale 
!  of the Smagorinsky anisotropic viscosity. 
!  </DATA> 
!  <DATA NAME="vel_micom_iso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity. 
!  </DATA> 
!  <DATA NAME="vel_micom_aniso" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic viscosity. 
!  </DATA> 
!
!  <DATA NAME="visc_crit_scale" UNITS="dimensionless" TYPE="real">
!  Scaling factor used to determine the critical viscosity, above which 
!  the viscosity is not allowed to reach. 
!  Use visc_crit_scale < 1.0 for cases where the visc_crit from linear stability 
!  allows for still too large of a viscosity.  Use visc_crit_scale>1.0 when wish 
!  to allow for larger viscosity. Default is visc_crit_scale=1.0.
!  </DATA> 
!
!  <DATA NAME="equatorial_zonal" TYPE="real">
!  Orient the anisotropic friction within a latitudinal band according 
!  to zonal direction. 
!  </DATA> 
!  <DATA NAME="equatorial_zonal_lat" TYPE="real">
!  Latitudinal band to use the zonal friction orientation. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_iso" TYPE="real">
!  Velocity scale that is used for computing the MICOM isotropic viscosity 
!  within a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_vel_micom_aniso" TYPE="real">
!  Velocity scale that is used for computing the MICOM anisotropic 
!  viscosity within a user specified equatorial band. 
!  </DATA> 
!  <DATA NAME="eq_lat_micom" TYPE="real">
!  Equatorial latitude band (degrees) within which the MICOM viscosity
!  is set according to eq_vel_micom_iso and eq_vel_micom_aniso.
!  </DATA> 
!  <DATA NAME="bottom_5point" TYPE="logical">
!  To alleviate problems with small partial cells, it is often necessary
!  to reduce the operator to the traditional 5-point Laplacian at the 
!  ocean bottom.  This logical implements this mixing. 
!  Default bottom_5point=.false. 
!  </DATA> 
!  <DATA NAME="vel_micom_bottom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM viscosity for 
!  5point Laplacian at the bottom. 
!  </DATA> 
!
!  <DATA NAME="ncar_boundary_scaling" TYPE="logical">
!  To enhance the velocity scale used in western boundaries 
!  for the isotropic and anisotropic  background viscosities, 
!  we compute a scaling using the algorithm from the laplacian
!  NCAR anisotropic scheme.
!  Default ncar_boundary_scaling=.false.
!  </DATA> 
!  <DATA NAME="ncar_rescale_power" UNITS="dimensionless" TYPE="integer">
!  For determining rescaling of the viscosity so to enhance the 
!  friction near the western boundaries. Default ncar_rescale_power=1.
!  </DATA> 
!  <DATA NAME="ncar_vconst_4" UNITS="1/cm" TYPE="real">
!  Inverse damping length for exponential falloff of the velocity scale 
!  as move eastward away from western boundary. Default ncar_vconst_4=2.e-8.
!  </DATA> 
!  <DATA NAME="ncar_vconst_5" UNITS="dimensionless" TYPE="integer">
!  For determining number of grid points in boundary calculation.
!  Default ncar_vconst_5=3.
!  </DATA> 
!
!</NAMELIST>
!
use constants_mod,       only: pi, radius, epsln
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use fms_mod,             only: open_namelist_file, check_nml_error
use fms_mod,             only: write_version_number, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains, mpp_global_field, mpp_global_min, mpp_global_max
use mpp_mod,             only: mpp_error, mpp_max, mpp_sum, mpp_pe, mpp_chksum 

use ocean_domains_mod,    only: get_local_indices
use ocean_obc_mod,        only: ocean_obc_enhance_visc_back
use ocean_operators_mod,  only: BAY, BAX, BDX_EU, BDY_NU, FDX_U, FDY_U, FMX, FMY
use ocean_parameters_mod, only: missing_value, oneeigth
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_thickness_type, ocean_velocity_type, ocean_options_type
use ocean_util_mod,       only: write_timestamp
use ocean_workspace_mod,  only: wrk1, wrk2, wrk1_v

implicit none

private

public ocean_bihgen_friction_init
public bihgen_friction
public bihgen_viscosity_check
public bihgen_reynolds_check

private ncar_boundary_scale 
private BDX_EU_smag
private BDY_NU_smag

real :: k_smag_iso          = 0.0    ! smag scaling coeff for isotropic viscosity (dimensionless) 
real :: k_smag_aniso        = 0.0    ! smag scaling coeff for anisotripic viscosity (dimensionless) 
real :: vel_micom_iso       = 0.0    ! background scaling velocity for isotropic viscosity (m/sec) 
real :: vel_micom_aniso     = 0.0    ! background scaling velocity for anisotropic viscosity (m/sec) 
real :: eq_vel_micom_iso    = 0.2    ! background scaling velocity (m/sec) within equatorial band for isotropic visc
real :: eq_vel_micom_aniso  = 0.0    ! background scaling velocity (m/sec) within equatorial band for anisotropic visc
real :: eq_lat_micom        = 0.0    ! equatorial latitude band for micom (degrees)
real :: vel_micom_bottom    = 0.01   ! velocity scale for determining viscosity at the bottom 
real :: equatorial_zonal_lat= 0.0    ! latitudinal band for orienting the friction along zonal direction
real :: visc_crit_scale     = 1.0    ! constant for setting the critical viscosity on the instability constraint 
logical :: equatorial_zonal = .false.! zonally orient anisotropic friction w/i equatorial_zonal_lat
logical :: bottom_5point    =.false.  ! for bottom Laplacian 5point mixing to avoid problems with thin partial cells. 

! for the western boundary scaling of the background viscosity
! nx and ny are number of global points in generalized x and y 
! directions
logical :: ncar_boundary_scaling = .false.
real    :: ncar_vconst_4         = 2.e-8
integer :: ncar_vconst_5         = 3  
integer :: ncar_rescale_power    = 1
integer :: nx, ny 

! for diagnostics 
integer :: id_aiso          =-1
integer :: id_aaniso        =-1
integer :: id_along         =-1
integer :: id_across        =-1
integer :: id_aiso_back     =-1
integer :: id_aaniso_back   =-1
integer :: id_along_back    =-1
integer :: id_across_back   =-1
integer :: id_visc_crit_bih =-1
integer :: id_bih_fric_u    =-1
integer :: id_bih_fric_v    =-1
integer :: id_horz_bih_diss =-1
integer :: id_ncar_rescale  =-1
logical :: used

! time step 
real ::  dtime = 0.0  ! time step used for friction tendency (2*dtuv if threelevel, dtuv if twolevel) 

#include <ocean_memory.h>


#ifdef MOM4_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed)    :: fsmag_iso       ! (m^4) combination of terms for computing isotropic smag visc
real, dimension(isd:ied,jsd:jed)    :: fsmag_aniso     ! (m^4) combination of terms for computing anisotropic smag visc
real, dimension(isd:ied,jsd:jed,nk) :: aiso_back       ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed,nk) :: aaniso_back     ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(isd:ied,jsd:jed,nk) :: aiso            ! isotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(isd:ied,jsd:jed,nk) :: ncar_rescale       ! viscosity scale (cm^2/s) based on NCAR western boundary scaling 

real, dimension(isd:ied,jsd:jed)    :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed)    :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(isd:ied,jsd:jed,nk) :: massqc      ! (1/4)*dxu*dyu*rho_dzu (kg) time-independent quarter-cell mass
real, dimension(isd:ied,jsd:jed)    :: visc_crit   ! critical value of the viscosity (m^4/sec) for linear stability 
real, dimension(isd:ied,jsd:jed)    :: visc_bottom ! Laplacian viscosity at the bottom (m^2/sec)

real, dimension(isd:ied,jsd:jed,0:1,0:1,2) :: stress   !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(isd:ied,jsd:jed,2)         :: tmplap
real, dimension(isd:ied,jsd:jed,2)         :: tmpfdelx
real, dimension(isd:ied,jsd:jed,2)         :: tmpfdely
real, dimension(isd:ied,jsd:jed)           :: tmp
real, dimension(isc:iec,jsc:jec)           :: cos2theta
real, dimension(isc:iec,jsc:jec)           :: sin2theta

#else

real, dimension(:,:), allocatable   :: fsmag_iso   ! (m^4) combination of terms for computing isotropic smag visc
real, dimension(:,:), allocatable   :: fsmag_aniso ! (m^4) combination of terms for computing anisotropic smag visc
real, dimension(:,:,:), allocatable :: aiso_back   ! minimum allowable isotropic viscosity (m^4/s)
real, dimension(:,:,:), allocatable :: aaniso_back ! minimum allowable anisotropic viscosity (m^4/s)
real, dimension(:,:,:), allocatable :: aiso        ! isotropic viscosity (m^4/s) averaged to U cell for diagnostics
real, dimension(:,:,:), allocatable :: ncar_rescale   ! viscosity scale (cm^2/s) based on NCAR western boundary scaling 

real, dimension(:,:), allocatable   :: daur_dxur   ! 1/(dxu*dxu*dyu) (m^-3)
real, dimension(:,:), allocatable   :: daur_dyur   ! 1/(dxu*dyu*dyu) (m^-3)
real, dimension(:,:,:), allocatable :: massqc      ! (1/4)*dxu*dyu*rho_dzu (kg) time-independent quarter-cell mass
real, dimension(:,:), allocatable   :: visc_crit   ! critical value of the viscosity (m^4/sec) for linear stability 
real, dimension(:,:), allocatable   :: visc_bottom ! Laplacian viscosity at the bottom (m^2/sec)

real, dimension(:,:,:,:,:), allocatable :: stress  !stress tensor: last index 1=stress_xx, 2=stress_xy
real, dimension(:,:,:), allocatable     :: tmplap
real, dimension(:,:,:), allocatable     :: tmpfdelx
real, dimension(:,:,:), allocatable     :: tmpfdely
real, dimension(:,:), allocatable       :: tmp
real, dimension(:,:), allocatable       :: cos2theta
real, dimension(:,:), allocatable       :: sin2theta

#endif

real, dimension(:,:), allocatable   :: aiso_back_obc   ! for enhancing visc next to OBC boundaries 
real, dimension(:,:), allocatable   :: aaniso_back_obc ! for enhancing visc next to OBC boundaries 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

character(len=256) :: version=&
     '=>Using: ocean_bihgen_friction.F90 ($Id: ocean_bihgen_friction.F90,v 16.0.2.1.84.1 2009/10/10 00:42:16 nnz Exp $)'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc              = .FALSE.

namelist /ocean_bihgen_friction_nml/ use_this_module, debug_this_module,           &
          k_smag_iso, k_smag_aniso, vel_micom_iso,                                 &
          vel_micom_aniso, eq_vel_micom_iso, eq_vel_micom_aniso, eq_lat_micom,     &
          vel_micom_bottom, bottom_5point, equatorial_zonal, equatorial_zonal_lat, &
          visc_crit_scale,                                                         &
          ncar_boundary_scaling, ncar_rescale_power, ncar_vconst_4, ncar_vconst_5

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bihgen_friction_init">

! <DESCRIPTION>
! Initialize the horizontal biharmonic friction module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that viscosity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_bihgen_friction_init(Grid, Domain, Time, Ocean_options, d_time, &
                                      obc, use_bihgen_friction, debug)

  type(ocean_grid_type),    target, intent(in)    :: Grid
  type(ocean_domain_type),  target, intent(in)    :: Domain
  type(ocean_time_type),            intent(in)    :: Time
  type(ocean_options_type),         intent(inout) :: Ocean_options
  real,                             intent(in)    :: d_time
  logical,                          intent(in)    :: obc
  logical,                          intent(inout) :: use_bihgen_friction
  logical,                optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: i, j, k
  real    :: coeff_iso, coeff_aniso, dxdy

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_bihgen_friction_mod(ocean_bihgen_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_bihgen_friction_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_bihgen_friction_nml)  
  write (stdlogunit,ocean_bihgen_friction_nml)
  ierr = check_nml_error(io_status,'ocean_bihgen_friction_nml')
  call close_file(ioun)

#ifndef MOM4_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
#endif

  if(use_this_module) then 
      call mpp_error(NOTE, '==> NOTE: USING ocean_bihgen_friction_mod.')
      Ocean_options%horz_bih_friction = 'Used general horizontal biharmonic friction.'
      use_bihgen_friction=.true.
  else 
      call mpp_error(NOTE, '==> NOTE: NOT using ocean_bihgen_friction_mod.')
      Ocean_options%horz_bih_friction = 'Did NOT use horizontal biharmonic friction.'
      use_bihgen_friction=.false.
      return 
  endif

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running with debug_this_module=.true. so will print lots of checksums.'  
  endif 

  Grd => Grid
  Dom => Domain
  nx  =  Grd%ni
  ny  =  Grd%nj

  dtime = d_time

  write(stdoutunit,'(a)') ' ' 

  write(stdoutunit,'(/a,f10.2)') &
   '==> Note from ocean_bihgen_friction_mod: using forward time step of (secs)', dtime 

  have_obc       = obc
  if (have_obc) then 
     write(stdoutunit,'(/a,f10.2)') '==> Note from ocean_bihgen_friction_mod: considering obc' 
  endif 

  if (Dom%xhalo > 1 .or. Dom%yhalo > 1) then
    write (stdoutunit,'(/1x,a)') ' ==> NOTE: General biharmonic friction allows xhalo=yhalo=1.'
  endif

  if(bottom_5point) then
    write(stdoutunit,'(/)')  
    write(stdoutunit,'(/1x,a)') ' ==> NOTE: Will make horizontal friction to a 5point Laplacian on the bottom'
    write(stdoutunit,'(a)')     '     This helps alleviate numerical problems with thin bottom partial cells.'
  endif 

  if(k_smag_iso > 0.0) then 
     write( stdoutunit,'(a)')' Computing horziontal isotropic biharmonic viscosity via Smagorinsky.'
  endif
  if(k_smag_iso==0.0)  then 
    write( stdoutunit,*)' Setting horzizontal isotropic biharmonic Smagorinsky viscosity to zero.'
  endif 
  if(k_smag_aniso > 0.0) then 
     write( stdoutunit,'(a)')' Computing horzizontal anisotropic biharmonic viscosity via Smagorinsky.'
     write( stdoutunit,'(1x,a)')' ==>WARNING: anisotropic biharmonic Smagorinsky not well tested.'
  endif
  if(k_smag_aniso==0.0) then 
    write( stdoutunit,*)' Setting horzizontal anisotropic biharmonic Smagorinsky viscosity to zero.'
  endif 
  if(vel_micom_iso > 0.0) then 
     write( stdoutunit,'(a)')' Computing background horzizontal biharmonic isotropic viscosity via MICOM.'
  endif
  if(vel_micom_iso==0.0) then 
    write( stdoutunit,*)' Setting background horzizontal biharmonic isotropic viscosity to zero.'
  endif 
  if(vel_micom_aniso > 0.0) then
     write( stdoutunit,'(a)')' Computing background horzizontal anisotropic biharmonic viscosity via MICOM.'
     write( stdoutunit,'(1x,a)')' ==> WARNING: anisotropic biharmonic background viscosity not well tested.'
  endif 
  if(vel_micom_aniso==0.0) then
     write(stdoutunit,'(a)')' Setting background horziontal biharmonic anisotropic viscosity to zero.' 
  endif 
  if(k_smag_aniso > 2.0*k_smag_iso) then
     call mpp_error(FATAL, &
     '==>Error: Smag horizontal bih anisotropic visc too large. Must be < 2*isotropic-visc.')
  endif  
  if(vel_micom_aniso > 2.0*vel_micom_iso) then
    call mpp_error(FATAL, &
     '==>Error: Background bih horz aniso visc too large. Must be < 2.0*background iso visc.')
  endif  
  if(eq_lat_micom > 0) then 
    write( stdoutunit,'(a)')'Setting different background horz bih viscosity w/i equatorial zone.'
  endif 
  if(equatorial_zonal) then
    write( stdoutunit,'(a)') &
    'If using anisotropic biharmonic friction, zonally orient friction w/i latitudinal'
    write( stdoutunit,'(a)') &
    'band of width ',equatorial_zonal_lat,' degrees.'
  endif 


#ifndef MOM4_STATIC_ARRAYS
  
  allocate (fsmag_iso(isd:ied,jsd:jed))
  allocate (fsmag_aniso(isd:ied,jsd:jed))
  allocate (ncar_rescale(isd:ied,jsd:jed,nk))
  allocate (aiso_back(isd:ied,jsd:jed,nk))
  allocate (aaniso_back(isd:ied,jsd:jed,nk))
  allocate (aiso(isd:ied,jsd:jed,nk)) 
  allocate (daur_dxur(isd:ied,jsd:jed))
  allocate (daur_dyur(isd:ied,jsd:jed))
  allocate (massqc(isd:ied,jsd:jed,nk))
  allocate (visc_crit(isd:ied,jsd:jed))
  allocate (visc_bottom(isd:ied,jsd:jed))
  allocate (stress(isd:ied,jsd:jed,0:1,0:1,2))  !stress tensor: last index 1=stress_xx, 2=stress_xy
  allocate (tmplap(isd:ied,jsd:jed,2))    
  allocate (tmpfdelx(isd:ied,jsd:jed,2))    
  allocate (tmpfdely(isd:ied,jsd:jed,2))   
  allocate (tmp(isd:ied,jsd:jed))
  allocate (cos2theta(isc:iec,jsc:jec))
  allocate (sin2theta(isc:iec,jsc:jec))
  
#endif

  allocate (aiso_back_obc(isd:ied,jsd:jed))
  allocate (aaniso_back_obc(isd:ied,jsd:jed))
 
  tmp(:,:) = 0.0
  tmplap(:,:,:) = 0.0
  
  ! Some commonly used grid arrays 
  daur_dxur(:,:) = 0.0
  daur_dyur(:,:) = 0.0
  daur_dxur(:,:) = Grd%daur(:,:)*Grd%dxur(:,:)
  daur_dyur(:,:) = Grd%daur(:,:)*Grd%dyur(:,:)

  fsmag_iso(:,:)   = 0.0
  fsmag_aniso(:,:) = 0.0
  coeff_iso        = 0.125*(k_smag_iso/pi)**2
  coeff_aniso      = 0.125*(k_smag_aniso/pi)**2
  fsmag_iso(:,:)   = Grd%umask(:,:,1)*coeff_iso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4
  fsmag_aniso(:,:) = Grd%umask(:,:,1)*coeff_aniso* &
                     ((2.0*Grd%dxu(:,:)*Grd%dyu(:,:))/(epsln + Grd%dxu(:,:) + Grd%dyu(:,:)))**4


  aiso_back_obc(:,:)   = 0.0
  aaniso_back_obc(:,:) = 0.0
  visc_bottom(:,:)     = 0.0
  do j=jsd,jed
    do i=isd,ied 
      visc_bottom(i,j)   = Grd%umask(i,j,1)*vel_micom_bottom* &
                           ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))
      aiso_back_obc(i,j)     = Grd%umask(i,j,1)*vel_micom_iso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      aaniso_back_obc(i,j)   = Grd%umask(i,j,1)*vel_micom_aniso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      if(abs(Grd%yu(i,j)) < eq_lat_micom) then
        aiso_back_obc(i,j)   = Grd%umask(i,j,1)*eq_vel_micom_iso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
        aaniso_back_obc(i,j) = Grd%umask(i,j,1)*eq_vel_micom_aniso* &
                               ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
      endif
    enddo  
  enddo

  aiso(:,:,:)        = 0.0
  aiso_back(:,:,:)   = 0.0
  aaniso_back(:,:,:) = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied 

           aiso_back(i,j,k)     = Grd%umask(i,j,k)*vel_micom_iso* &
                ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           aaniso_back(i,j,k)   = Grd%umask(i,j,k)*vel_micom_aniso* &
                ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3

           if(abs(Grd%yu(i,j)) < eq_lat_micom) then
               aiso_back(i,j,k)   = Grd%umask(i,j,k)*eq_vel_micom_iso* &
                    ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
               aaniso_back(i,j,k) = Grd%umask(i,j,k)*eq_vel_micom_aniso* &
                    ((2.0*Grd%dxu(i,j)*Grd%dyu(i,j))/(epsln + Grd%dxu(i,j) + Grd%dyu(i,j)))**3
           endif

        enddo
     enddo
  enddo

  ! rescale the background viscosities according to NCAR boundary scaling  
  if(ncar_boundary_scaling) then 
     write(stdoutunit,'(a)') ' Note: rescaling background viscosities so they are larger in western boundaries.'
     call ncar_boundary_scale(Time)
  endif 

  ! enhance background viscosity near open boundaries if needed
  if (have_obc) then 
      call ocean_obc_enhance_visc_back(aiso_back_obc, 'bih')
      call ocean_obc_enhance_visc_back(aaniso_back_obc, 'bih')
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               if(aiso_back_obc(i,j)   > aiso_back(i,j,k))   aiso_back(i,j,k)   = aiso_back_obc(i,j)
               if(aaniso_back_obc(i,j) > aaniso_back(i,j,k)) aaniso_back(i,j,k) = aaniso_back_obc(i,j)
            enddo
         enddo
      enddo
  endif


  ! critical value of viscosity, above which 2-dim linear stability is not satisfied.
  ! visc_crit taken from equation (18.26) of Griffies book. 
  visc_crit(:,:) = 0.0
  do j=jsc,jec
    do i=isc,iec
      dxdy = 1.0/( 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j)) )
      visc_crit(i,j) =  visc_crit_scale*oneeigth*dxdy**2/(dtime+epsln)
    enddo
  enddo

  id_visc_crit_bih = register_static_field ('ocean_model', 'visc_crit_bih', &
                     Grd%vel_axes_uv(1:2), 'critical viscosity', 'm^4/sec', &
                     missing_value=missing_value, range=(/0.0,1.e20/))
  if (id_visc_crit_bih > 0) used = send_data (id_visc_crit_bih, visc_crit(isc:iec,jsc:jec), &
                            Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,1))


  ! ensure that background viscosities are not too large 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(aiso_back(i,j,k)   > visc_crit(i,j))  aiso_back(i,j,k)   = visc_crit(i,j)
           if(aaniso_back(i,j,k) > visc_crit(i,j))  aaniso_back(i,j,k) = visc_crit(i,j)
        enddo
     enddo
  enddo

  ! ensure viscosities maintain proper relative values so that friction dissipates 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(aaniso_back(i,j,k) >= 2.0*aiso_back(i,j,k) .and. aaniso_back(i,j,k) > 0.0)  then 
               write(stdoutunit,'(a,i4,a,i4,a,i3,a)')'Violating bih iso/aniso constraint at (',i,',',j,')'
               aaniso_back(i,j,k) = 1.9*aiso_back(i,j,k)
           endif
        enddo
     enddo
  enddo

  call mpp_update_domains (aiso_back(:,:,:),    Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:),  Dom%domain2d)    

  ! viscosity for diagnostic output 
  id_aiso    = register_diag_field ('ocean_model', 'aiso_bih', Grd%vel_axes_uv(1:3), &
               Time%model_time,'U-cell isotropic bih visc', 'm^4/sec',               &
               missing_value=-10.0, range=(/-10.0,1.e20/),                           &
               standard_name='ocean_momentum_xy_biharmonic_diffusivity')
  id_aaniso  = register_diag_field ('ocean_model', 'aaniso_bih', Grd%vel_axes_uv(1:3), &
               Time%model_time,'U-cell  anisotropic bih visc', 'm^4/sec',              &
               missing_value=missing_value, range=(/-10.0,1.e20/))
  id_along  = register_diag_field ('ocean_model', 'along_bih', Grd%vel_axes_uv(1:3),   &
              Time%model_time,'U-cell along-stream bih visc', 'm^4/sec',               &
              missing_value=missing_value, range=(/-10.0,1.e20/))
  id_across  = register_diag_field ('ocean_model', 'across_bih', Grd%vel_axes_uv(1:3), &
               Time%model_time, 'U-cell cross-stream bih visc', 'm^4/sec',             &
               missing_value=missing_value, range=(/-10.0,1.e20/))
  id_bih_fric_u = register_diag_field ('ocean_model', 'bih_fric_u', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict on u-zonal',   &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))
  id_bih_fric_v = register_diag_field ('ocean_model', 'bih_fric_v', Grd%vel_axes_uv(1:3), &
                  Time%model_time, 'Thickness and rho wghtd horz bih frict on v-merid',   &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1.e20,1.e20/))

  id_aiso_back   = register_static_field('ocean_model', 'aiso_bih_back', Grd%vel_axes_uv(1:3),   &
                    'U-cell background aiso bih visc', 'm^4/sec',                                &
                    missing_value=missing_value, range=(/-10.0,1.e10/))
  id_aaniso_back = register_static_field('ocean_model', 'aaniso_bih_back', Grd%vel_axes_uv(1:3), &
                    'U-cell background aaniso bih visc', 'm^4/sec',                              &
                    missing_value=missing_value, range=(/-10.0,1.e10/))
  id_along_back   = register_static_field('ocean_model', 'along_bih_back', Grd%vel_axes_uv(1:3), &
                    'U-cell background along-stream bih visc', 'm^4/sec',                        &
                    missing_value=missing_value, range=(/-10.0,1.e10/))
  id_across_back = register_static_field('ocean_model', 'across_bih_back', Grd%vel_axes_uv(1:3), &
                   'U-cell background cross-stream bih visc', 'm^4/sec',                         &
                   missing_value=missing_value, range=(/-10.0,1.e10/))

  if (id_aiso_back > 0) then 
    used = send_data (id_aiso_back, aiso_back(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 
  if (id_aaniso_back > 0) then 
    used = send_data (id_aaniso_back, aaniso_back(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 
  if (id_along_back > 0) then 
    used = send_data (id_along_back, aiso_back(isc:iec,jsc:jec,:)+0.5*aaniso_back(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 
  if (id_across_back > 0) then 
    used = send_data (id_across_back, aiso_back(isc:iec,jsc:jec,:)-0.5*aaniso_back(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

  id_horz_bih_diss = register_diag_field ('ocean_model', 'horz_bih_diss', Grd%vel_axes_uv(1:3),      &
                  Time%model_time, 'Energy dissipation from horizontal biharmonic friction', 'W/m^2',&
                  missing_value=missing_value, range=(/-1.e20,1.e20/))


end subroutine ocean_bihgen_friction_init
! </SUBROUTINE>  NAME="ocean_bihgen_friction_init"


!#######################################################################
! <SUBROUTINE NAME="bihgen_friction">
!
! <DESCRIPTION>
! This subroutine computes the time tendency for horizontal 
! velocity (i.e., the acceleration) from horizontal biharmonic friction.  
! The algorithm is derived from a functional approach that ensures kinetic 
! energy is consistenty dissipated for all flow configurations. 
! The triad do-loops are expanded in order to enhance the 
! ability of cache-based machines to keep most of the variables 
! on-cache.  
! 
! Fundamental to the scheme are the rates of horizontal deformation  <BR/>
! horizontal tension = DT = (dy)(u/dy)_x - (dx)(v/dx)_y              <BR/>
! horizontal strain  = DS = (dx)(u/dx)_y + (dy)(v/dy)_x              <BR/>
! Units of the tension and strain are sec^-1.
!
! Four tensions and four strains are computed for each velocity point, <BR/>
! corresponding to the four triads surrounding the point.              <BR/>
! The following notation is used to distinguish the triads:            <BR/>   
! (0,1)=northwest triad  (1,1)=northeast triad,                        <BR/> 
! (0,0)=southwest triad, (1,0)=southeast triad                        
!
! A triad contributes when at least one of its velocities is             
! not a land point.  In order to obtain the correct tension           
! and strain next to boundaries, tension and strain should not be       
! masked with umask. 
!
! As shown in Griffies and Hallberg (2000), 
! a biharmonic operator with a nonconstant viscosity is guaranteed to 
! dissipate kinetic energy *only* when using the sqrt of the biharmonic
! viscosity at each of the two stages of the algorithm. 
! The sqrt approach is employed here.  
!
! </DESCRIPTION>
!
subroutine bihgen_friction(Time, Thickness, Velocity, bih_viscosity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  real, dimension(isd:,jsd:), intent(inout) :: bih_viscosity
  logical,                    intent(in)    :: energy_analysis_step

  real, dimension(isc:iec,jsc:jec,0:1,0:1) :: aiso_smag
  real, dimension(isc:iec,jsc:jec,0:1,0:1) :: aaniso_smag

  integer :: i, j, k, n, ip, jq, tau, taum1

  real :: u1_m10, u2_m10
  real :: u1_10,  u2_10
  real :: u1_00,  u2_00
  real :: u1_01,  u2_01
  real :: u1_0m1, u2_0m1 
  real :: dxuer_ip0, dxuer_ip1
  real :: dyunr_jq0, dyunr_jq1
  real :: usqrd, vsqrd, umagr
  real :: tension, strain, deform, delta 
  real :: tension_metric, strain_metric

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bihgen_friction_mod (horz_bih_friction): needs to be initialized')
  endif 

  if(.not. use_this_module) then 

      if(energy_analysis_step) then 
          do n=1,2
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      Velocity%wrkv(i,j,k,n) = 0.0
                   enddo
                enddo
             enddo
          enddo
      endif

      return 
  endif

  stress(:,:,:,:,:) = 0.0
  wrk1_v(:,:,:,:)   = 0.0
  wrk1(:,:,:)       = 0.0
  wrk2(:,:,:)       = 0.0

  taum1 = Time%taum1
  tau   = Time%tau
  
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           massqc(i,j,k) = 0.25*Grd%dau(i,j)*Thickness%rho_dzu(i,j,k,tau)
        enddo
     enddo
  enddo


  do k=1,nk


     ! Laplacian part of algorithm
     do j=jsc,jec
        do i=isc,iec

           tension_metric = -Velocity%u(i,j,k,1,taum1)*Grd%dh2dx(i,j) + Velocity%u(i,j,k,2,taum1)*Grd%dh1dy(i,j)  
           strain_metric  = -Velocity%u(i,j,k,1,taum1)*Grd%dh1dy(i,j) - Velocity%u(i,j,k,2,taum1)*Grd%dh2dx(i,j)  

           if(equatorial_zonal .and. abs(Grd%yu(i,j)) <= equatorial_zonal_lat) then  
             sin2theta(i,j) = 0.0
             cos2theta(i,j) = 1.0
           else 
             usqrd = Velocity%u(i,j,k,1,taum1)*Velocity%u(i,j,k,1,taum1)
             vsqrd = Velocity%u(i,j,k,2,taum1)*Velocity%u(i,j,k,2,taum1)
             umagr = 1.0/(epsln + usqrd + vsqrd)
             sin2theta(i,j) = 2.0*Velocity%u(i,j,k,1,taum1)*Velocity%u(i,j,k,2,taum1)*umagr
             cos2theta(i,j) = (usqrd-vsqrd)*umagr
           endif  

          ! compute stress tensor components for the four triads.  
          ! expanded do-loops are faster.
          ! dimensions[aiso_smag]=m^2/s^0.5
          ! dimensions[tension and strain]=1/s
          ! dimensions[stress]=m^2/s^1.5

           ip=0 ; dxuer_ip0 = Grd%dxuer(i+ip-1,j)
           ip=1 ; dxuer_ip1 = Grd%dxuer(i+ip-1,j)
           jq=0 ; dyunr_jq0 = Grd%dyunr(i,j+jq-1)
           jq=1 ; dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  ; u1_m10 = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_m10 = Velocity%u(i+ip,j+jq,k,2,taum1)
           ip=1  ; jq=0  ; u1_10  = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_10  = Velocity%u(i+ip,j+jq,k,2,taum1) 
           ip=0  ; jq=0  ; u1_00  = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_00  = Velocity%u(i+ip,j+jq,k,2,taum1) 
           ip=0  ; jq=1  ; u1_01  = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_01  = Velocity%u(i+ip,j+jq,k,2,taum1)
           ip=0  ; jq=-1 ; u1_0m1 = Velocity%u(i+ip,j+jq,k,1,taum1) ; u2_0m1 = Velocity%u(i+ip,j+jq,k,2,taum1)

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           deform  = sqrt(tension**2 + strain**2)
           aiso_smag(i,j,ip,jq)    = min(visc_crit(i,j), aiso_back(i,j,k)   + fsmag_iso(i,j)*deform)
           aaniso_smag(i,j,ip,jq)  = min(visc_crit(i,j), aaniso_back(i,j,k) + fsmag_aniso(i,j)*deform)
           aiso_smag(i,j,ip,jq)    = sqrt(aiso_smag(i,j,ip,jq))
           aaniso_smag(i,j,ip,jq)  = sqrt(aaniso_smag(i,j,ip,jq))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ! viscosities for diagnostics 
           aiso(i,j,k) = 0.25*(  aiso_smag(i,j,0,0)**2   + aiso_smag(i,j,1,0)**2    &
                               + aiso_smag(i,j,0,1)**2   + aiso_smag(i,j,1,1)**2)
           wrk1(i,j,k) = 0.25*(  aaniso_smag(i,j,0,0)**2 + aaniso_smag(i,j,1,0)**2  &
                               + aaniso_smag(i,j,0,1)**2 + aaniso_smag(i,j,1,1)**2)

        enddo
     enddo

     call mpp_update_domains (stress(:,:,:,:,:), Dom%domain2d)    

     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains
     ! dimensions[tmpfdelx and tmpfdely]=kg*m^2/s^1.5

     tmpfdelx(:,:,:) = 0.0
     do j=jsc,jec
        do i=isd,iec
           tmpfdelx(i,j,1) = (stress(i,j,1,0,1)   + stress(i,j,1,1,1)  )*massqc(i,j,k)    &
                +(stress(i+1,j,0,0,1) + stress(i+1,j,0,1,1))*massqc(i+1,j,k)
           tmpfdelx(i,j,2) = (stress(i,j,1,0,2)   + stress(i,j,1,1,2)  )*massqc(i,j,k)    &  
                +(stress(i+1,j,0,0,2) + stress(i+1,j,0,1,2))*massqc(i+1,j,k)              
        enddo
     enddo

     tmpfdely(:,:,:) = 0.0
     do j=jsd,jec
        do i=isc,iec  
           tmpfdely(i,j,1) = (stress(i,j,0,1,2)   + stress(i,j,1,1,2)  )*massqc(i,j,k)    &                      
                +(stress(i,j+1,0,0,2) + stress(i,j+1,1,0,2))*massqc(i,j+1,k) 
           tmpfdely(i,j,2) = (stress(i,j,0,1,1)   + stress(i,j,1,1,1)  )*massqc(i,j,k)    &                      
                +(stress(i,j+1,0,0,1) + stress(i,j+1,1,0,1))*massqc(i,j+1,k)
        enddo
     enddo

     ! compute laplacian operator
     ! dimensions [tmplap=m/s^1.5]
     do n=1,2
        tmp(:,:) = BDX_EU_smag(tmpfdelx(:,:,n))+(3-2*n)*BDY_NU_smag(tmpfdely(:,:,n))    

        ! do not use rho_dzur since that is at taup1 and need to divide by rho_dau(tau)
        do j=jsc,jec
           do i=isc,iec
              tmplap(i,j,n) = tmp(i,j)*Grd%umask(i,j,k)/Thickness%rho_dzu(i,j,k,tau)
              wrk2(i,j,k)   = wrk2(i,j,k) - 4.0*massqc(i,j,k)*tmplap(i,j,n)**2 
           enddo
        enddo

     enddo
     call mpp_update_domains (tmplap(:,:,:), Dom%domain2d)    

     ! Second part of the iteration
     ! tmplap[m/s^1.5] replaces velocity[m/s]
     ! dimensions[aiso_smag]=m^2/s^0.5
     ! dimensions[tension and strain]=1/s^1.5
     ! dimensions[stress]=m^2/s^2

     stress(:,:,:,:,:)   = 0.0

     do j=jsc,jec
        do i=isc,iec

           tension_metric = -tmplap(i,j,1)*Grd%dh2dx(i,j) + tmplap(i,j,2)*Grd%dh1dy(i,j)  
           strain_metric  = -tmplap(i,j,1)*Grd%dh1dy(i,j) - tmplap(i,j,2)*Grd%dh2dx(i,j)  

          ! compute stress tensor components for the four triads.  
          ! expanded do-loops are quicker. 

           ip=0 ; dxuer_ip0 = Grd%dxuer(i+ip-1,j)
           ip=1 ; dxuer_ip1 = Grd%dxuer(i+ip-1,j)
           jq=0 ; dyunr_jq0 = Grd%dyunr(i,j+jq-1)
           jq=1 ; dyunr_jq1 = Grd%dyunr(i,j+jq-1)

           ip=-1 ; jq=0  ; u1_m10 = tmplap(i+ip,j+jq,1) ; u2_m10 = tmplap(i+ip,j+jq,2)
           ip=1  ; jq=0  ; u1_10  = tmplap(i+ip,j+jq,1) ; u2_10  = tmplap(i+ip,j+jq,2) 
           ip=0  ; jq=0  ; u1_00  = tmplap(i+ip,j+jq,1) ; u2_00  = tmplap(i+ip,j+jq,2) 
           ip=0  ; jq=1  ; u1_01  = tmplap(i+ip,j+jq,1) ; u2_01  = tmplap(i+ip,j+jq,2)
           ip=0  ; jq=-1 ; u1_0m1 = tmplap(i+ip,j+jq,1) ; u2_0m1 = tmplap(i+ip,j+jq,2)

           ip=0 ; jq=0 
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=1 ; jq=0 
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq0*(u2_00-u2_0m1)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq0*(u1_00-u1_0m1)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           ip=0 ; jq=1
           tension =  dxuer_ip0*(u1_00-u1_m10)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip0*(u2_00-u2_m10)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

           jq=1 ; ip=1
           tension =  dxuer_ip1*(u1_10-u1_00)-dyunr_jq1*(u2_01-u2_00)+tension_metric
           strain  =  dxuer_ip1*(u2_10-u2_00)+dyunr_jq1*(u1_01-u1_00)+strain_metric
           delta   =  0.5*(strain*cos2theta(i,j) - tension*sin2theta(i,j))
           stress(i,j,ip,jq,1) = aiso_smag(i,j,ip,jq)*tension + aaniso_smag(i,j,ip,jq)*delta*sin2theta(i,j)
           stress(i,j,ip,jq,2) = aiso_smag(i,j,ip,jq)*strain  - aaniso_smag(i,j,ip,jq)*delta*cos2theta(i,j)

        enddo
     enddo

     call mpp_update_domains (stress(:,:,:,:,:), Dom%domain2d)    

     ! sum the stress tensor components over the triads
     ! loops designed to reduce calls to mpp_update_domains
     ! dimensions[tmpfdelx and tmpfdely]=kg*m^2/s^2

     tmpfdelx(:,:,:) = 0.0
     do j=jsc,jec
        do i=isd,iec
           tmpfdelx(i,j,1) = (stress(i,j,1,0,1)   + stress(i,j,1,1,1)  )*massqc(i,j,k)    &
                +(stress(i+1,j,0,0,1) + stress(i+1,j,0,1,1))*massqc(i+1,j,k)
           tmpfdelx(i,j,2) = (stress(i,j,1,0,2)   + stress(i,j,1,1,2)  )*massqc(i,j,k)    &  
                +(stress(i+1,j,0,0,2) + stress(i+1,j,0,1,2))*massqc(i+1,j,k)              
        enddo
     enddo

     tmpfdely(:,:,:) = 0.0
     do j=jsd,jec
        do i=isc,iec  
           tmpfdely(i,j,1) = (stress(i,j,0,1,2)   + stress(i,j,1,1,2)  )*massqc(i,j,k)    &                      
                +(stress(i,j+1,0,0,2) + stress(i,j+1,1,0,2))*massqc(i,j+1,k) 
           tmpfdely(i,j,2) = (stress(i,j,0,1,1)   + stress(i,j,1,1,1)  )*massqc(i,j,k)    &                      
                +(stress(i,j+1,0,0,1) + stress(i,j+1,1,0,1))*massqc(i,j+1,k)
        enddo
     enddo

    ! compute acceleration from horizontal biharmonic friction
    ! dimensions[tmp]=kg*/(m*s^2)
    ! dimensions[friction]=(kg/m^3)*(m^2/s^2)=N/m^2
     do n=1,2

        tmp(:,:) = BDX_EU_smag(tmpfdelx(:,:,n))+(3-2*n)*BDY_NU_smag(tmpfdely(:,:,n)) 
        do j=jsc,jec
           do i=isc,iec
              wrk1_v(i,j,k,n) = -tmp(i,j)*Grd%umask(i,j,k)
           enddo
        enddo

        ! reduce to 5-point laplacian at bottom to avoid problems with thin bottom partial cells
        if (bottom_5point) then
            tmp(:,:) = BDX_EU(visc_bottom(:,:)*FMX(Thickness%rho_dzu(:,:,k,tau))*FDX_U(Velocity%u(:,:,k,n,taum1))) &
                      +BDY_NU(visc_bottom(:,:)*FMY(Thickness%rho_dzu(:,:,k,tau))*FDY_U(Velocity%u(:,:,k,n,taum1)))
            do j=jsc,jec
               do i=isc,iec
                  if(k==Grd%kmu(i,j)) wrk1_v(i,j,k,n) = tmp(i,j)*Grd%umask(i,j,k)
               enddo
            enddo
        endif

     enddo !end of n-loop

  enddo  !end of k-loop


  if(energy_analysis_step) then 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%wrkv(i,j,k,n) = wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

  else 

      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      ! vertically averaged viscosity at U-cell centre
      bih_viscosity(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               bih_viscosity(i,j) = bih_viscosity(i,j)  &
               + Grd%umask(i,j,k)*Thickness%rho_dzu(i,j,k,tau)*aiso(i,j,k)
            enddo
         enddo
      enddo
      do j=jsc,jec
         do i=isc,iec
            bih_viscosity(i,j) = Grd%umask(i,j,1)*bih_viscosity(i,j)/(epsln+Thickness%mass_u(i,j,tau))
         enddo
      enddo
      call mpp_update_domains (bih_viscosity(:,:), Dom%domain2d)    


      if (id_aiso > 0)    used = send_data (id_aiso, aiso(:,:,:), &
           Time%model_time, rmask=Grd%umask(:,:,:),               &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

      if (id_aaniso > 0)  used = send_data (id_aaniso, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%umask(:,:,:),                 &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

      if (id_along > 0)   used = send_data (id_along,  aiso(:,:,:)+0.5*wrk1(:,:,:),  &
           Time%model_time, rmask=Grd%umask(:,:,:),                                  &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

      if (id_across > 0)  used = send_data (id_across, aiso(:,:,:)-0.5*wrk1(:,:,:),  &
           Time%model_time, rmask=Grd%umask(:,:,:),                                  &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

      if (id_horz_bih_diss > 0)  then
          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   wrk2(i,j,k) = wrk2(i,j,k)*Grd%daur(i,j)
                enddo
             enddo
          enddo
          used = send_data (id_horz_bih_diss, wrk2(:,:,:),&
          Time%model_time, rmask=Grd%umask(:,:,:),        &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif 

      if (id_bih_fric_u > 0) used = send_data(id_bih_fric_u, wrk1_v(isc:iec,jsc:jec,:,1), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))

      if (id_bih_fric_v > 0) used = send_data(id_bih_fric_v, wrk1_v(isc:iec,jsc:jec,:,2), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))


  endif

  if(debug_this_module) then
      write(stdoutunit,*) ' ' 
      write(stdoutunit,*) 'From ocean_bihgen_friction_mod: friction chksums'
      call write_timestamp(Time%model_time)
      write(stdoutunit,*) 'rho_dzu(tau)       = ', &
                         mpp_chksum(Thickness%rho_dzu(isc:iec,jsc:jec,:,tau)*Grd%umask(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'bihgen friction(1) = ', &
                         mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,1)*Grd%umask(isc:iec,jsc:jec,:))
      write(stdoutunit,*) 'bihgen friction(2) = ', &
                         mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,2)*Grd%umask(isc:iec,jsc:jec,:))
  endif 

end subroutine bihgen_friction
! </SUBROUTINE> NAME="bihgen_friction"




!#######################################################################
! <SUBROUTINE NAME="ncar_boundary_scale">
!
! <DESCRIPTION>
!
!     Recale the background viscosities to be larger in the western 
!     boundary regions.  The algorithm is taken directly from the 
!     anisotropic_ncar routine in ocean_lapgen_friction.F90.
!
!   NOTE: The nearest western boundary computations are done along the
!         model i-grid lines. Therefore, viscosity based on these are 
!         only approximate in the high Northern Hemisphere when using 
!         generalized coordinates with coordinate pole(s) shifted onto 
!         land. 
!
! </DESCRIPTION>
!
subroutine ncar_boundary_scale(Time)

  type(ocean_time_type), intent(in) :: Time

  integer :: i, j, k
  integer :: n, ncount
  integer :: is, ie, ip1, ii
  integer :: index, indexo
  integer, dimension(nx) :: iwp
  integer, dimension(nx) :: nwbp_global

  real, dimension(nx,ny)            :: dxtn_global
  real, dimension(nx,ny)            :: kmu_global
  real, dimension(isd:ied,jsd:jed)  :: kmu_tmp
  real, dimension(isd:ied,jsd:jed)  :: ncar_rescale2d
  real, dimension(nx)               :: dist_g 

  real :: dist_max
  real :: ncar_rescale_min
  real :: huge=1e20

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bihgen_friction_mod (ncar_velocity_scale): module must be initialized')
  endif 

  dist_max         = 1.e10  ! distance for ACC region (cm)
  dxtn_global(:,:) = 0.0
  kmu_global(:,:)  = 0.0
  kmu_tmp(:,:)     = Grd%kmu(:,:)
  ncar_rescale(:,:,:) = 0.0
  ncar_rescale2d(:,:) = 0.0

  call mpp_global_field(Dom%domain2d,kmu_tmp,kmu_global)  
  call mpp_global_field(Dom%domain2d,Grd%dxtn,dxtn_global)

  do k=1,nk
     do j=jsc,jec

        ! determine nearest western boundary
        ncount         = 0
        nwbp_global(:) = 0
        iwp(:)         = 0
        dist_g(:)      = 0.0

        do i=1,Grd%ni
           ip1 = i+1
           if(i==Grd%ni) then 
               if(Grd%cyclic_x) then 
                   ip1=1
               else
                   ip1=i  
               endif
           endif
           if(kmu_global(i,j+Dom%joff)<k .and. kmu_global(ip1,j+Dom%joff) >= k) then 
               ncount      = ncount+1
               iwp(ncount) = i
           endif
        enddo

        if (ncount > 0) then
            do n=1,ncount-1
               is = iwp(n)
               ie = iwp(n+1)-1
               do i=is,ie
                  nwbp_global(i)=is
               enddo
            enddo
            do i=1,Grd%ni
               if(nwbp_global(i)==0) nwbp_global(i) = iwp(ncount)
            enddo
        endif


        ! determine distance (cm) to nearest western boundary
        do i=1,Grd%ni

           index  = nwbp_global(i)
           indexo = index + ncar_vconst_5
           if ( index .eq. 0) then
               dist_g(i) = dist_max
           elseif ( i .ge. index  .and.  i .le. indexo ) then
               dist_g(i) = 0.0
           elseif ( (i .gt. indexo) ) then
               dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i-1)
           elseif ( i .lt. index ) then
               if (indexo .le. Grd%ni) then
                   if (i .eq. 1) then
                       dist_g(i) = 0.0
                       do ii=indexo+1,Grd%ni
                          dist_g(i)=dxtn_global(ii,j+Dom%joff) + dist_g(i)
                       enddo
                       dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i)
                   else
                       dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i-1)
                   endif
               else
                   if (i .le. (indexo - Grd%ni)) then
                       dist_g(i) = 0.0
                   else
                       dist_g(i) = dxtn_global(i,j+Dom%joff)+dist_g(i-1)
                   endif
               endif
           endif

        enddo

        ! convert distance in mom4 (m) to POP1.4 (cm)
        dist_g(:) = 100.0*dist_g(:)
        do i=isc,iec
           ncar_rescale(i,j,k) = Grd%umask(i,j,k)*(1.0+exp(-(ncar_vconst_4*dist_g(i+Dom%ioff))**2))
        enddo

     enddo ! j
  enddo    ! k

  ! the only depth dependence to ncar_rescale arises from absence or presence of 
  ! topography.  so to compute the min ncar_rescale, we only need to look at k=1.  
  ! also, to avoid getting ncar_rescale=0.0 due to land, we set ncar_rescale2d
  ! to a huge number if over land. 
  do j=jsc,jec
     do i=isc,iec
        if(Grd%umask(i,j,1)== 0.0) then 
           ncar_rescale2d(i,j) = huge
        else 
           ncar_rescale2d(i,j) = ncar_rescale(i,j,1) 
        endif 
     enddo
  enddo
  ncar_rescale_min = mpp_global_min(Dom%domain2d,ncar_rescale2d(:,:))

  ! check for error
  if(ncar_rescale_min==0.0) then
    call mpp_error(FATAL, &
    '==>Error from ocean_bihgen_friction_mod(ncar_boundary_scale): ncar_rescale_min=0.  Something wrong')
  else
    write(stdoutunit,'(a,e14.6)') 'From ncar_boundary_scale, minimum ncar_rescale =  ',ncar_rescale_min
  endif 

  ! rescale the background viscosities 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           ncar_rescale(i,j,k) = (ncar_rescale(i,j,k)/ncar_rescale_min)**ncar_rescale_power
           aiso_back(i,j,k)    = aiso_back(i,j,k)  *ncar_rescale(i,j,k)
           aaniso_back(i,j,k)  = aaniso_back(i,j,k)*ncar_rescale(i,j,k)
        enddo
     enddo
  enddo

  id_ncar_rescale = register_static_field('ocean_model', 'ncar_rescale', Grd%vel_axes_uv(1:3),&
                    'rescaling used for the background viscosity', 'dimensionless',           &
                     missing_value=missing_value, range=(/-10.0,1.e10/))
  if (id_ncar_rescale > 0) then 
    used = send_data (id_ncar_rescale, ncar_rescale(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 


  ! need these viscosities in the halo regions
  call mpp_update_domains (aiso_back(:,:,:),   Dom%domain2d)    
  call mpp_update_domains (aaniso_back(:,:,:), Dom%domain2d)    


end subroutine ncar_boundary_scale
! </SUBROUTINE>  NAME="ncar_boundary_scale"



!#######################################################################
! <FUNCTION NAME="BDX_EU_smag">
!
! <DESCRIPTION>
! Compute backwards Derivative in X of a quantity defined on the east 
! face of a U-cell. Slightly modified version of BDX_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i-1/2,j).
!
! BDX_EU_smag changes dimensions by m^-3 
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! field defined on the east face of a U-cell
! </IN>
!
function BDX_EU_smag(a)

  real, dimension(isd:,jsd:), intent(in) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDX_EU_smag
  integer :: i, j

  do j=jsd,jed
    do i=isd+1,ied
      BDX_EU_smag(i,j) = (Grd%dyue_dxuer(i,j)*a(i,j) - Grd%dyue_dxuer(i-1,j)*a(i-1,j))*daur_dyur(i,j)
    enddo
    BDX_EU_smag(isd,j) = 0.0
  enddo

end function BDX_EU_smag
! </FUNCTION> NAME="BDX_EU_smag"


!#######################################################################
! <FUNCTION NAME="BDY_NU_smag">
!
! <DESCRIPTION>
! Compute backwards Derivative in Y of a quantity defined on the north
! face of a U-cell. Slightly modified version of BDY_EU used in
! ocean_operators.F90. If input is a(i,j) then output is defined 
! at (i,j-1/2).
!
! BDY_EU_smag changes dimensions by m^-3 
!
! </DESCRIPTION>
!
! <IN NAME="a" TYPE="real" DIM="(isd:ied,jsd:jed)">
! field defined on the north face of a U-cell
! </IN>
!
function BDY_NU_smag(a)

  real, dimension(isd:,jsd:), intent(in) :: a
  real, dimension(isd:ied,jsd:jed)       :: BDY_NU_smag
  integer :: i, j

  do j=jsd+1,jed
    do i=isd,ied
      BDY_NU_smag(i,j) = (Grd%dxun_dyunr(i,j)*a(i,j) - Grd%dxun_dyunr(i,j-1)*a(i,j-1))*daur_dxur(i,j)
    enddo
  enddo
  BDY_NU_smag(:,jsd) = 0.0

end function BDY_NU_smag
! </FUNCTION> NAME="BDY_EU_smag"


!#######################################################################
! <SUBROUTINE NAME="bihgen_viscosity_check">
!
! <DESCRIPTION>
! Subroutine to perform linear stability check for the biharmonic
! operator given a value for the horizontal biharmonic viscosity.
! </DESCRIPTION>
!
subroutine bihgen_viscosity_check

  integer :: num, i, j, k
  integer, save :: unit=6, max_num=3
  real :: fudge

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bihgen_friction_mod (bihgen_viscosity_check): needs initialization')
  endif 

  fudge=1.001  ! to eliminate roundoffs causing aiso=visc_crit+epsln

  write (stdoutunit,'(//60x,a/)') ' Excessive horizontal biharmonic friction summary:'
  write (stdoutunit,'(1x,a/)')'Locations (if any) where viscosity exceeds conservative linear stability constraint'
  num = 0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if (Grd%umask(i,j,k) > 0.0) then            
          if (aiso(i,j,k) > fudge*visc_crit(i,j) .and. num < max_num ) then
            num = num + 1
            write (unit,9600) 'aiso(',aiso(i,j,k), visc_crit(i,j),i,j,k,Grd%xu(i,j),Grd%yu(i,j),Grd%zt(k),mpp_pe()
          endif
        endif 
      enddo
    enddo
  enddo
  call mpp_sum (num)
  if (num > max_num) write (stdoutunit,*) ' (a max of ',max_num,' violations per pe are listed)'

9600  format(/' Warning: ',a,es16.8,' m^4/s) exceeds max value (',es16.8,') at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m, pe=',i3,')')


end subroutine bihgen_viscosity_check
! </SUBROUTINE> NAME="bihgen_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="bihgen_reynolds_check">
!
! <DESCRIPTION>
! Subroutine to compute the biharmonic grid Reynolds number.  Large 
! Reynolds numbers indicate regions where solution may experience 
! some grid noise due to lack of enough horizontal friction. 
! </DESCRIPTION>
!
subroutine bihgen_reynolds_check(Time, Velocity)

  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  real :: rame, ramn, reyx, reyy
  real :: reynx0, reyny0, reynx,reyny, reynu,reynv, reynmu,reynmv
  integer :: ireynx,jreynx,kreynx, ireyny,jreyny,kreyny
  integer :: i, j, k, tau
  integer, save :: unit=6

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_bihgen_friction_mod (gen_bih_reynolds_check): nees to be initialized')
  endif 

  if(.not. use_this_module) return 

  tau = Time%tau

  ! look for max grid reynolds numbers using velocities at "tau".

  ireynx=isc; jreynx=jsc; kreynx=1; reynx=0.0
  ireyny=isc; jreyny=jsc; kreyny=1; reyny=0.0
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        rame = 1.0/(aiso(i,j,k) + epsln)
        ramn = 1.0/(aiso(i,j,k) + epsln)

        reyx = abs(Velocity%u(i,j,k,1,tau)*Grd%dxu(i,j)**3)*rame
        if (reyx > reynx) then
          ireynx = i
          jreynx = j
          kreynx = k
          reynx  = reyx
          reynu  = Velocity%u(i,j,k,1,tau)
          reynmu = 1.0/rame
        endif
        reyy = abs(Velocity%u(i,j,k,2,tau)*Grd%dyu(i,j)**3)*ramn
        if (reyy > reyny) then
          ireyny = i
          jreyny = j
          kreyny = k
          reyny  = reyy
          reynv  = Velocity%u(i,j,k,2,tau)
          reynmv = 1.0/ramn
        endif

      enddo
    enddo
  enddo
  write (stdoutunit,'(/60x,a/)') ' Horizontal biharmonic Reynolds number summary'

  reynx0 = reynx
  reyny0 = reyny
  call mpp_max(reynx)
  call mpp_max(reyny)
  
  if (reynx == reynx0) then
    write (unit,10300) reynx, ireynx, jreynx, kreynx, &
                       Grd%xu(ireynx,jreynx), Grd%yu(ireynx,jreynx), Grd%zt(kreynx), reynu, reynmu
  endif
  
  if (reyny == reyny0) then
    write (unit,10400) reyny, ireyny, jreyny, kreyny, &
                       Grd%xu(ireyny,jreyny), Grd%yu(ireyny,jreyny), Grd%zt(kreyny), reynv, reynmv
  endif

  10300 format (1x,'Maximum zonal Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), U =',es9.2, ', mix =',e9.2)
  10400 format (1x,'Maximum merid Re =',es9.2,' at (i,j,k) = (',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.2,'), V =',es9.2, ', mix =',e9.2)


end subroutine bihgen_reynolds_check
! </SUBROUTINE> NAME="bihgen_reynolds_check"

end module ocean_bihgen_friction_mod
      
      




