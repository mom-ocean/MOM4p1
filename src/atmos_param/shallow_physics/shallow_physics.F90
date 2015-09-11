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
module shallow_physics_mod

use  fms_mod, only: open_namelist_file, file_exist,   &
                    close_file, check_nml_error,      &
                    error_mesg, FATAL, WARNING,       &
                    write_version_number, stdlog,     &
                    mpp_pe, mpp_root_pe

use time_manager_mod, only: time_type

implicit none
private

!========================================================================

public :: shallow_physics_init,    &
          shallow_physics,         &
          shallow_physics_end

interface shallow_physics_init
   module procedure shallow_physics_init_1d, shallow_physics_init_2d
end interface
!========================================================================
! version information 
character(len=128) :: version = '$Id: shallow_physics.F90,v 17.0 2009/07/21 02:58:12 fms Exp $'
character(len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'
!========================================================================

real, allocatable, dimension(:,:) :: h_eq

real    :: kappa_m, kappa_t

logical :: module_is_initialized = .false.

!========================================================================
! namelist 

real    :: fric_damp_time  = -20.0
real    :: therm_damp_time = -10.0
real    :: h_0             = 3.e04
real    :: h_monsoon       = 2.e04
real    :: lon_monsoon     =  90.0
real    :: lat_monsoon     =  25.0
real    :: width_monsoon   =  15.0
real    :: h_itcz          = 1.e05
real    :: width_itcz      =  4.0
logical :: no_forcing      = .false.

namelist /shallow_physics_nml/ fric_damp_time, therm_damp_time, &
                               h_0, h_monsoon, width_monsoon,   &
                               lon_monsoon, lat_monsoon,        &
                               width_itcz, h_itcz, no_forcing

contains

!========================================================================

subroutine shallow_physics_init_2d (axes, Time, lon, lat) 
integer, intent(in) :: axes(4)
type(time_type), intent(in) :: Time
real, intent(in) :: lon(:,:), lat(:,:)  ! longitude and latitude in radians

integer :: i, j, unit, ierr, io, logunit
real    :: xm, ym, dm, di
real    :: lon_m, lat_m, width_m, width_i, deg2rad

! cannot initialize the module more than once
  if (module_is_initialized) then
    call error_mesg ('shallow_physics_init', &
                     'module has already been initialized ', FATAL)
  endif

! read the namelist

  if (file_exist('input.nml')) then
    unit = open_namelist_file ()
    ierr=1
    do while (ierr /= 0)
      read  (unit, nml=shallow_physics_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'shallow_physics_nml')
    enddo
    10 call close_file (unit)
  endif

! write version info and namelist to logfile

  call write_version_number (version, tagname)
  logunit = stdlog()
  write(logunit,nml=shallow_physics_nml)

! damping times < 0 are in days (convert to seconds)

  if (fric_damp_time  < 0.0)  fric_damp_time = -  fric_damp_time*86400
  if (therm_damp_time < 0.0) therm_damp_time = - therm_damp_time*86400

! compute damping coefficients

  kappa_m = 0.0
  kappa_t = 0.0
  if ( fric_damp_time .ne. 0.0) kappa_m = 1./fric_damp_time
  if (therm_damp_time .ne. 0.0) kappa_t = 1./therm_damp_time

! global storage

  allocate ( h_eq(size(lon,1),size(lon,2)) )

! convert namelist variables in degrees to radians

  deg2rad = acos(0.0)/90.
  lon_m = lon_monsoon * deg2rad
  lat_m = lat_monsoon * deg2rad
  width_m = width_monsoon * deg2rad
  width_i = width_itcz    * deg2rad

! compute constants

  do j = 1, size(lon,2)
  do i = 1, size(lon,1)
     xm = (lon(i,j) - lon_m)/(width_m*2.)
     ym = (lat(i,j) - lat_m)/width_m
     dm =  xm*xm + ym*ym
     di = (lat(i,j)/width_i)**2
     h_eq(i,j) = h_0 + h_monsoon*max(1.e-10, exp(-dm)) + h_itcz*exp(-di)
  enddo
  enddo

  module_is_initialized = .true.

end subroutine shallow_physics_init_2d

!=======================================================================

subroutine shallow_physics_init_1d (axes, Time, lon, lat) 
integer, intent(in) :: axes(4)
type(time_type), intent(in) :: Time
real, intent(in) :: lon(:), lat(:)  ! longitude and latitude in radians

real, dimension(size(lon),size(lat)) :: lon2, lat2

   lon2 = spread(lon,2,size(lat))
   lat2 = spread(lat,1,size(lon))
   call shallow_physics_init_2d (axes, Time, lon2, lat2)

end subroutine shallow_physics_init_1d

!=======================================================================

subroutine shallow_physics ( is, ie, js, je, timelev, dt, Time,    &
                             um, vm, hm, u, v, h, u_dt, v_dt, h_dt )

integer,         intent(in) :: is, ie, js, je, timelev
real,            intent(in) :: dt
type(time_type), intent(in) :: Time
real, intent(in)   , dimension(is:ie,js:je) :: um, vm, hm, u, v, h
real, intent(inout), dimension(is:ie,js:je) :: u_dt, v_dt, h_dt

integer :: i, j

  if (.not.module_is_initialized) then
    call error_mesg ('shallow_physics', &
                     'module has not been initialized ', FATAL)
  endif

  if (no_forcing) return

! choose which time level is used to compute forcing

  select case (timelev)
     case(-1)
         ! previous time level (tau-1)
         do j = js, je
         do i = is, ie
            u_dt(i,j) = u_dt(i,j) - kappa_m *  um(i,j)
            v_dt(i,j) = v_dt(i,j) - kappa_m *  vm(i,j)
            h_dt(i,j) = h_dt(i,j) - kappa_t * (hm(i,j) - h_eq(i,j))
         enddo
         enddo
     case(0)
         ! current time level (tau)
         do j = js, je
         do i = is, ie
            u_dt(i,j) = u_dt(i,j) - kappa_m *  u(i,j)
            v_dt(i,j) = v_dt(i,j) - kappa_m *  v(i,j)
            h_dt(i,j) = h_dt(i,j) - kappa_t * (h(i,j) - h_eq(i,j))
         enddo
         enddo
     case(+1)
         ! next time level (tau+1)
         do j = js, je
         do i = is, ie
            u_dt(i,j) = u_dt(i,j)*(1.-kappa_m*dt) - kappa_m *  um(i,j)
            v_dt(i,j) = v_dt(i,j)*(1.-kappa_m*dt) - kappa_m *  vm(i,j)
            h_dt(i,j) = h_dt(i,j)*(1.-kappa_t*dt) - kappa_t * (hm(i,j)-h_eq(i,j))
         enddo
         enddo
     case default
         call error_mesg ('shallow_physics', &
                          'invalid value for timelev argument', FATAL)
  end select

end subroutine shallow_physics

!======================================================================

subroutine shallow_physics_end

  if (.not.module_is_initialized) then
    call error_mesg ('shallow_physics_end', &
                     'module has not been initialized ', WARNING)
    return
  endif

! release global storage

  deallocate ( h_eq )

  module_is_initialized = .false.

end subroutine shallow_physics_end

!======================================================================

end module shallow_physics_mod
