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
module ocean_bbc_mod
!  
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matthew Harrison
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Set bottom boundary conditions 
!</OVERVIEW>
!
!<DESCRIPTION>
! Set bottom boundary conditions 
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
! S.M. Griffies, 2006: Elements of mom4p1
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_bbc_nml">
!
!  <DATA NAME="cdbot" UNITS="dimensionless" TYPE="real">
! Dimensionless coefficient for quadratic bottom drag. 
!  </DATA> 
!
!  <DATA NAME="bmf_implicit" TYPE="logical">
!  For incorporating the bottom momentum drag implicitly in time.
!  This code is fresh and not fully tested.  
!  Default is bmf_implicit=.false. 
!  </DATA> 
!
!  <DATA NAME="cdbot_law_of_wall" TYPE="logical">
!  For determining bottom drag coefficient using a constant roughness length.
!  Will take maximum between cdbot and the computed value using law of
!  wall log-profile.  This option of use when have very very 
!  refined vertical resolution (say on order of meters) near the bottom.
!  Terrain following coordinates should use this option since they generally 
!  have very refined vertical grid spacing on topography. 
!  Default is cdbot_law_of_wall=.false. 
!  </DATA> 
!  <DATA NAME="law_of_wall_rough_length" UNITS="metre" TYPE="real">
!  Bottom roughness length.  Default is law_of_wall_rough_length=0.01m, following
!  the default used in the Princeton Ocean Model (POM). This value 
!  corresponds to "Law of Wall" physics.  
!  </DATA> 
!
!  <DATA NAME="cdbot_roughness_length" TYPE="logical">
!  For determining bottom drag coefficient using a map of the roughness length.
!  This approach is more relevant for coarse models
!  than the constant roughness length used in the cdbot_law_of_wall option. 
!  Default is cdbot_roughness_length=.false. 
!  </DATA> 
!
!  <DATA NAME="uresidual" UNITS="m/s" TYPE="real">
!  Residual bottom velocity due to unresolved fluctuations (e.g., waves and tides)
!  that contribute to bottom dissipation.  Should be set to zero when running 
!  with explicit representation of tidal forcing and when waves are well resolved.
!  Default is uresidual=.05.
!  </DATA> 
!
!  <DATA NAME="uvmag_max" UNITS="m/s" TYPE="real">
!  Maximum magnitude of the bottom velocity used to compute the bottom 
!  momentum drag.  Default is uvmag_max=10.0.
!  </DATA> 
!
!  <DATA NAME="bmf_max" UNITS="N/m2" TYPE="real">
!  Maximum magnitude of the bottom momentum drag.
!  Default is bmf_max=1.0.
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging purposes.
!  </DATA> 
!</NAMELIST>
!
use diag_manager_mod,     only: register_diag_field, register_static_field, send_data
use fms_mod,              only: open_namelist_file, check_nml_error, write_version_number
use fms_mod,              only: close_file, read_data
use mpp_mod,              only: mpp_error, FATAL, WARNING, stdout, stdlog
use mpp_domains_mod,      only: mpp_update_domains, BGRID_NE

use ocean_parameters_mod, only: missing_value, TERRAIN_FOLLOWING
use ocean_parameters_mod, only: onefourth, rho0, rho0r, cp_ocean 
use ocean_types_mod,      only: ocean_velocity_type, ocean_domain_type
use ocean_types_mod,      only: ocean_grid_type, ocean_prog_tracer_type
use ocean_types_mod,      only: ocean_time_type, ocean_thickness_type, ocean_options_type
use ocean_domains_mod,    only: get_local_indices

implicit none

private

#include <ocean_memory.h>

integer :: num_prog_tracers
integer :: id_gamma_bmf =-1
integer :: id_drag_coeff=-1
integer :: id_cdbot     =-1
integer :: id_bmf_u     =-1
integer :: id_bmf_v     =-1
logical :: used

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

character(len=128) :: version=&
     '$Id: ocean_bbc.F90,v 16.0.2.3.48.1.32.1 2009/10/10 00:41:51 nnz Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

public :: ocean_bbc_init
public :: get_ocean_bbc

logical :: module_is_initialized = .false.

real    :: vonkarman=0.4  ! dimensionless 
real    :: vonkarman2          
real    :: rho0_cdbot
real    :: law_of_wall_rough_length_r

real    :: cp_r
real    :: uresidual2
real    :: small=1.e-20

integer :: index_temp=-1
integer :: index_salt=-1

integer :: id_geo_heat 
integer :: bbc_geothermal 
real, allocatable, dimension(:,:) :: cdbot_array       ! dimensionless drag coefficient
real, allocatable, dimension(:,:) :: cdbot_lowall      ! dimensionless drag coefficient from law of wall
real, allocatable, dimension(:,:) :: drag_coeff        ! rho0 * dimensionless drag coefficient
real, allocatable, dimension(:,:) :: roughness_length  ! for computing bottom drag coefficient
real, allocatable, dimension(:,:) :: geo_heat          ! geothermal heating
real, allocatable, dimension(:,:) :: data       

! nml variables 
logical :: bmf_implicit             = .false. 
logical :: debug_this_module        = .false. 
logical :: cdbot_law_of_wall        = .false. 
logical :: cdbot_roughness_length   = .false.  
logical :: use_geothermal_heating   = .false. 
real    :: convert_geothermal       = .001     ! for converting units from mW to W
real    :: law_of_wall_rough_length = 0.01     ! metre
real    :: cdbot                    = 2.5e-3   ! dimensionless
real    :: uresidual                = .05      ! m/sec
real    :: cdbot_hi                 = 3.0e-3   ! hi-end of cdbot for cdbot_roughness_length 
real    :: cdbot_lo                 = 1.0e-3   ! lo-end of cdbot for cdbot_roughness_length 
real    :: cdbot_gamma              = 40.0     ! for setting exp decay for cdbot w/ cdbot_roughness_length 
real    :: uvmag_max                = 10.0     ! max velocity scale (m/s) for computing bmf
real    :: bmf_max                  = 1.0      ! max bmf (N/m^2)

 namelist /ocean_bbc_nml/ bmf_implicit, cdbot, uresidual, cdbot_law_of_wall, law_of_wall_rough_length, &
                         cdbot_roughness_length, use_geothermal_heating, convert_geothermal,           &
                         cdbot_hi, cdbot_lo, cdbot_gamma, uvmag_max, bmf_max, debug_this_module 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bbc_init">
!
! <DESCRIPTION>
! Initialize the bottom boundary condition module
! </DESCRIPTION>
!
subroutine ocean_bbc_init(Grid, Domain, Time, T_prog, Velocity, Ocean_options, vert_coordinate_type)

type(ocean_grid_type),   target, intent(in)    :: Grid
type(ocean_domain_type), target, intent(in)    :: Domain
type(ocean_time_type),           intent(in)    :: Time
type(ocean_prog_tracer_type),    intent(inout) :: T_prog(:)
type(ocean_velocity_type),       intent(inout) :: Velocity
type(ocean_options_type),        intent(inout) :: Ocean_options
integer,                         intent(in)    :: vert_coordinate_type

integer :: i,j,n, ioun, io_status, ierr
real    :: depth 

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 


module_is_initialized = .TRUE.

call write_version_number(version, tagname)

ioun = open_namelist_file()
read(ioun, ocean_bbc_nml, iostat=io_status)
write (stdoutunit,'(/)')
write (stdoutunit, ocean_bbc_nml)
write (stdlogunit, ocean_bbc_nml)
ierr = check_nml_error(io_status,'ocean_bbc_nml')
call close_file(ioun)

#ifndef MOM4_STATIC_ARRAYS
call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
nk = Grid%nk
#endif

Dom => Domain
Grd => Grid

num_prog_tracers = size(T_prog(:))
index_temp=-1;index_salt=-1
do n=1,num_prog_tracers
   if (T_prog(n)%name == 'temp') index_temp = n
   if (T_prog(n)%name == 'salt') index_salt = n
enddo

do n=1, num_prog_tracers
#ifndef MOM4_STATIC_ARRAYS
   allocate(T_prog(n)%btf(isd:ied,jsd:jed))
#endif
   T_prog(n)%btf = 0.0
enddo

allocate (cdbot_array(isd:ied,jsd:jed))
allocate (cdbot_lowall(isd:ied,jsd:jed))
allocate (drag_coeff(isd:ied,jsd:jed))
allocate (roughness_length(isd:ied,jsd:jed))
allocate (geo_heat(isd:ied,jsd:jed))
allocate (data(isd:ied,jsd:jed))

rho0_cdbot            = rho0*cdbot
vonkarman2            = vonkarman*vonkarman
drag_coeff(:,:)       = rho0*cdbot*Grd%umask(:,:,1) 
cdbot_array(:,:)      = cdbot*Grd%umask(:,:,1)
cdbot_lowall(:,:)     = cdbot*Grd%umask(:,:,1)
geo_heat(:,:)         = 0.0
data(:,:)             = 0.0
roughness_length(:,:) = 0.0

cp_r                       = 1.0/cp_ocean
uresidual2                 = uresidual**2
law_of_wall_rough_length_r = 1.0/law_of_wall_rough_length

if(cdbot_law_of_wall) then 
   write(stdoutunit,'(a)') &
   '==>Note from ocean_bbc_mod: Computing bottom drag coefficient from roughness length.'
   Ocean_options%bottom_roughness = 'Computed bottom drag coefficient from Law of Wall roughness length.'
elseif (cdbot_roughness_length) then 
   Ocean_options%bottom_roughness = 'Computed bottom drag from map of bottom roughness using law of wall.'
else 
   Ocean_options%bottom_roughness = 'Computed bottom drag from specified drag coefficient.'
endif 

if(bmf_implicit) then 
   call mpp_error(WARNING,'==>Warning from ocean_bbc_mod: bmf_implicit=.true. has not been fully tested.')
   write(stdoutunit,'(a)') &
   '==>Note from ocean_bbc_mod: Computing bottom drag implicitly in time.'
   Ocean_options%bmf_implicit = 'Computed bottom drag implicitly in time.'
   Velocity%bmf_implicit = .true.
else
   Ocean_options%bmf_implicit = 'Computed bottom drag explicitly in time.'
   Velocity%bmf_implicit = .false.
endif 


! compute drag coefficient using law of wall with a roughness length 
! read in from a map.  This approach is more relevant for large-scale 
! coarse models than the constant roughness length used in the 
! cdbot_law_of_wall option. 
if(cdbot_roughness_length) then

    ! read in topographic roughness; assumed on U-grid
    data=0.0
    call read_data('INPUT/roughness_cdbot.nc','roughness_cdbot', data, Domain%domain2d)
    write (stdoutunit,*) '==>ocean_bbc_mod: Completed read of topographic roughness length.'
    do j=jsc,jec
       do i=isc,iec
          roughness_length(i,j) = Grd%umask(i,j,1)*max(0.0,data(i,j))
       enddo
    enddo
    call mpp_update_domains(roughness_length(:,:), Dom%domain2d) 

    ! bottom drag from law of wall using roughness_length
    cdbot_lowall=0.0
    do j=jsd,jed
       do i=isd,ied
          depth = Grd%ht(i,j)                 
          if(depth > 0.0) then 
              cdbot_lowall(i,j)= (vonkarman/( alog(depth/(roughness_length(i,j)+small)) + small ) )**2
          endif
       enddo
    enddo

    ! bottom drag coefficient; the range of cdbot is
    ! cd_beta <= cdbot(i,j) <= cd_alpha
    do j=jsd,jed
       do i=isd,ied
          cdbot_array(i,j) = Grd%umask(i,j,1)*cdbot_hi*(1.0-exp(-cdbot_gamma*cdbot_lowall(i,j)))
          cdbot_array(i,j) = max(cdbot_lo,cdbot_array(i,j))
       enddo
    enddo

endif


if(use_geothermal_heating) then 
    write(stdoutunit,'(a)') &
         '==>Note from ocean_bbc_mod: Geothermal heating introduced at ocean bottom.'
    Ocean_options%geothermal_heating = 'Geothermal heating introduced at ocean bottom.'

    ! fill geothermal heating array 
    geo_heat = 0.0
    call read_data('INPUT/geothermal_heating.nc', 'geo_heat', data(:,:), &
                   Dom%domain2d, timelevel=1)

    ! make sure there are no spurious negative values.
    ! also convert to units appropriate for temperature btf. 
    do j=jsc,jec
       do i=isc,iec
          geo_heat(i,j) = Grd%tmask(i,j,1)*max(0.0,data(i,j))*cp_r*convert_geothermal
       enddo
    enddo

else
    Ocean_options%geothermal_heating = 'Did NOT introduce geothermal heating.'
endif


if(vert_coordinate_type == TERRAIN_FOLLOWING .and. .not. cdbot_law_of_wall) then 
   call mpp_error(WARNING,&
   '==>Warning from ocean_bbc_init: recommend law_of_wall=.true. for TERRAIN_FOLLOWING coordinates.')
endif 

id_drag_coeff = register_diag_field ('ocean_model', 'drag_coeff',         &
                Grd%vel_axes_uv(1:2), Time%model_time,                    &
                'Dimensionless bottom drag coefficient', 'dimensionless', &
                missing_value=missing_value, range=(/-1.0,1.e3/))
id_gamma_bmf      = register_diag_field ('ocean_model', 'gamma_bmf',   &
                Grd%vel_axes_uv(1:2), Time%model_time,                 &
                'Bottom drag factor rho0*cdbot*uvmag', 'kg/(m^2 sec)', &
                missing_value=missing_value, range=(/-1.0,1.e12/))
id_bmf_u      = register_diag_field ('ocean_model', 'bmf_u',  &
                Grd%vel_axes_uv(1:2), Time%model_time,        &
                'Bottom u-stress via bottom drag', 'N/m^2',   &
                missing_value=missing_value, range=(/-1.0,1.e2/))
id_bmf_v      = register_diag_field ('ocean_model', 'bmf_v',  &
                Grd%vel_axes_uv(1:2), Time%model_time,        &
                'Bottom v-stress via bottom drag', 'N/m^2',   &
                missing_value=missing_value, range=(/-1.0,1.e2/))

id_cdbot      = register_static_field ('ocean_model', 'cdbot',                    &
                Grd%vel_axes_uv(1:2), 'Bottom drag coefficient', 'dimensionless', &     
                missing_value=missing_value, range=(/-1.0,1e6/))
if (id_cdbot > 0) then 
    used = send_data (id_cdbot, cdbot_array(:,:), Time%model_time, &
           rmask=Grd%umask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif

id_geo_heat = register_static_field ('ocean_model', 'geo_heat',    &
              Grd%tracer_axes(1:2), 'Geothermal heating', 'W/m^2', &     
              missing_value=missing_value, range=(/-10.0,1e6/),    &
              standard_name='upward_geothermal_heat_flux_at_sea_floor')
if (id_geo_heat > 0) then 
    used = send_data (id_geo_heat, geo_heat(:,:)*cp_ocean, Time%model_time, &
           rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif


return

end subroutine ocean_bbc_init
! </SUBROUTINE> NAME="ocean_bbc_init"


!#######################################################################
! <SUBROUTINE NAME="get_ocean_bbc">
!
! <DESCRIPTION>
! Set bottom boundary conditions for velocity and tracer.
!
! Dimensions of bottom momentum flux are 
! N/m^2 = (kg/m^3)*(m^2/s^2).
!
! Note the use of rho0 for the conversion from m^2/s^2 to 
! (kg/m^3)*(m^2/s^2).  We do not know the precise value 
! of cdbot, so the rho0 approximate value is well within 
! our level of uncertainty.  No reason therefore to 
! use in situ rho for this conversion, even when using 
! non-Boussinesq version of mom4.   
!
!
! Note that bmf needs to be computed on the data domain since the 
! halo values are accessed in ocean_vert_gotm.F90.
!
! </DESCRIPTION>
!
subroutine get_ocean_bbc(Time, Thickness, Velocity, T_prog)

type(ocean_time_type),        intent(in)    :: Time
type(ocean_thickness_type),   intent(in)    :: Thickness
type(ocean_velocity_type),    intent(inout) :: Velocity
type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

integer :: i, j, kmu, n
integer :: taum1
real    :: uvmag, argument, distance 
real    :: umag, vmag

if (.not. module_is_initialized) then 
   call mpp_error(FATAL,'==>Error from ocean_bbc_mod (get_ocean_bbc): module must be initialized ')
endif 

taum1 = Time%taum1


! Bottom tracer flux defaulted to zero.
do n=1,num_prog_tracers
   do j=jsd,jed
      do i=isd,ied
         T_prog(n)%btf(i,j) = 0.0
      enddo
   enddo
enddo

! T_prog(index_temp)%btf<0 means heat enters ocean; hence the minus sign. 
if (use_geothermal_heating) then
    do j=jsc,jec
       do i=isc,iec
          T_prog(index_temp)%btf(i,j) = -geo_heat(i,j)
       enddo
    enddo
endif


! set bottom drag coefficient
drag_coeff(:,:) = 0.0
if(cdbot_law_of_wall) then 
    do j=jsd,jed
       do i=isd,ied
          kmu = Grd%kmu(i,j)
          if(kmu > 0) then 
              distance        = Thickness%depth_zwu(i,j,kmu) - Thickness%depth_zu(i,j,kmu)
              argument        = law_of_wall_rough_length_r*(distance + law_of_wall_rough_length)
              drag_coeff(i,j) = rho0*max(cdbot,vonkarman2/(log(argument))**2)
          endif
       enddo
    enddo
else 
    do j=jsd,jed
       do i=isd,ied
          drag_coeff(i,j) = rho0*cdbot_array(i,j)*Grd%umask(i,j,1) 
       enddo
    enddo
endif

! set quadratic bottom momentum drag
! units of bmf are N/m^2
! introduce ceilings for uvmag and bmf
! to minimize potential for instability. 
Velocity%bmf(:,:,:) = 0.0
do j=jsd,jed
   do i=isd,ied
      kmu = Grd%kmu(i,j)
      Velocity%gamma(i,j) = 0.0
      if (kmu > 0) then
         uvmag = sqrt(uresidual2 + Velocity%u(i,j,kmu,1,taum1)**2 + Velocity%u(i,j,kmu,2,taum1)**2)
         uvmag = min(uvmag,uvmag_max)
         Velocity%gamma(i,j)  = drag_coeff(i,j)*uvmag
         umag                 = abs(Velocity%u(i,j,kmu,1,taum1))
         vmag                 = abs(Velocity%u(i,j,kmu,2,taum1))
         Velocity%bmf(i,j,1)  = min(bmf_max, Velocity%gamma(i,j)*umag)*sign(1.0,Velocity%u(i,j,kmu,1,taum1))
         Velocity%bmf(i,j,2)  = min(bmf_max, Velocity%gamma(i,j)*vmag)*sign(1.0,Velocity%u(i,j,kmu,2,taum1))
      endif
   enddo
enddo


! send diagnostics 
if (id_drag_coeff > 0) then 
  used = send_data(id_drag_coeff, rho0r*drag_coeff(:,:), Time%model_time, &
         rmask=Grd%umask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 
if (id_gamma_bmf > 0) then 
  used = send_data(id_gamma_bmf, Velocity%gamma(:,:), Time%model_time, &
         rmask=Grd%umask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 
if (id_bmf_u > 0) then 
  used = send_data(id_bmf_u, Velocity%bmf(:,:,1), Time%model_time, &
         rmask=Grd%umask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 
if (id_bmf_v > 0) then 
  used = send_data(id_bmf_v, Velocity%bmf(:,:,2), Time%model_time, &
         rmask=Grd%umask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
endif 

return

end subroutine get_ocean_bbc
! </SUBROUTINE> NAME="get_ocean_bbc"

end module ocean_bbc_mod


