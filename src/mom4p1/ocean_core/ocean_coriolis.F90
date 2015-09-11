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
module ocean_coriolis_mod
!  
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
!</CONTACT>
!
!<CONTACT EMAIL="Tony.Rosati@noaa.gov"> A. Rosati 
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Coriolis acceleration 
!</OVERVIEW>

!<DESCRIPTION>
! This module computes Coriolis acceleration on a B-grid. 
!</DESCRIPTION>
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
! Elements of MOM4p1 (2005)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_coriolis_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.  
!  </DATA> 
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to add contributions from Coriolis force.   
!  </DATA> 
!
!  <DATA NAME="acor" TYPE="real">
!  acor=0.0 means explicit Coriolis force.  0.5 le acor le 1.0 means semi-implicit.
!  Semi-implicit removes dtuv time step constraint associated with intertial oscillations,
!  but it leads to Coriolis force affecting energy balances.  
!  If use two-level tendency discretization, then acor=0 is NOT allowed since the 
!  model will be linearly unstable with growth rate going as f*(delta time). 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: radius, radian, pi, epsln, omega
use diag_manager_mod, only: register_static_field, register_diag_field, send_data
use fms_mod,          only: write_version_number, mpp_error, FATAL
use fms_mod,          only: check_nml_error, close_file, open_namelist_file
use mpp_mod,          only: mpp_max, stdout, stdlog, mpp_chksum

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: TWO_LEVEL
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,      only: ocean_velocity_type, ocean_thickness_type
use ocean_util_mod,       only: write_timestamp
use ocean_workspace_mod,  only: wrk1, wrk2, wrk1_v  

implicit none

private

! for diagnostics 
logical :: used
integer :: id_coriolis  =-1
integer :: id_beta      =-1
integer :: id_beta_eff  =-1
integer :: id_cor_u     =-1
integer :: id_cor_v     =-1
integer :: id_hrho_cor_u=-1
integer :: id_hrho_cor_v=-1
integer :: id_ucori_impl=-1
integer :: id_vcori_impl=-1

#include <ocean_memory.h>

public ocean_coriolis_init
public coriolis_force
public coriolis_force_implicit

type(ocean_grid_type), pointer :: Grd =>NULL()

character(len=128) :: version = &
     '$Id: ocean_coriolis.F90,v 16.0.108.1 2009/10/10 00:41:52 nnz Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

real    :: dtime 
logical :: module_is_initialized = .FALSE.

logical :: debug_this_module = .FALSE.
logical :: use_this_module   = .true.
real    :: acor              = 0.5

namelist /ocean_coriolis_nml/  debug_this_module, use_this_module, acor 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_coriolis_init">
!
! <DESCRIPTION>
! Initialize the Coriolis module.
! </DESCRIPTION>
!
subroutine ocean_coriolis_init(Grid, Domain, Time, Time_steps, debug)

  type(ocean_grid_type),       intent(inout),target :: Grid
  type(ocean_domain_type),     intent(in),   target :: Domain  
  type(ocean_time_type),       intent(in)           :: Time
  type(ocean_time_steps_type), intent(inout)        :: Time_steps
  logical,                     intent(in), optional :: debug

  real    :: deg2m, sin1, cos1, y, fmax
  real    :: max_dt_for_inertial_oscillation
  integer :: i, j
  integer :: ioun, io_status, ierr

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_coriolis_mod (ocean_coriolis_init): module already initialized.')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read  (ioun, ocean_coriolis_nml,iostat=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_coriolis_nml)  
  write (stdlogunit, ocean_coriolis_nml)
  ierr = check_nml_error(io_status, 'ocean_coriolis_nml')
  call close_file(ioun)

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Note: running ocean_coriolis_mod with debug_this_module=.true.'  
  endif 
  if(.not. use_this_module) then 
    write(stdoutunit,'(a)') &
    '==>Warning: NOT adding Coriolis to acceleration. Simulation has no rotational effects.'  
  endif 

  Grd=> Grid
  dtime = Time_steps%dtime_u

  ! parameter determining how to time discretize Coriolis force 
  Time_steps%acor = acor
  
  if(Time_steps%tendency==TWO_LEVEL .and. acor==0.0) then  
     call mpp_error(FATAL,&
     '==>ocean_coriolis_mod: acor=0.0 w/ 2-level tendency is unstable. must set 0.5<=acor<=1.0')
  endif 

#ifndef MOM4_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
  allocate (Grid%f(isd:ied,jsd:jed))
  allocate (Grid%fstar(isd:ied,jsd:jed))
  allocate (Grid%beta(isd:ied,jsd:jed))
  allocate (Grid%beta_eff(isd:ied,jsd:jed))
#endif

  if (Grid%beta_plane .or. Grid%f_plane) then

    ! beta plane with f     = f0     + beta*y where f0    is at f_plane_latitude
    ! beta plane with fstar = fstar0 + beta*y where fstar is at f_plane_latitude
    ! if f_plane then beta  = 0

    if (Grid%f_plane) then
      Grid%beta(:,:) = 0.0
    else
      Grid%beta(:,:) = 2.0*omega*cos(Grid%f_plane_latitude*pi/180.0)/radius
    endif

    deg2m = radius/radian
    sin1  = sin(Grid%f_plane_latitude*pi/180.0)
    cos1  = cos(Grid%f_plane_latitude*pi/180.0)
    do j=jsd,jed
      do i=isd,ied
        y      = (Grid%yu(i,j)-Grid%f_plane_latitude)*deg2m
        Grid%f(i,j)     = 2.0*omega*sin1  + y*Grid%beta(i,j)
        Grid%fstar(i,j) = 2.0*omega*cos1  + y*Grid%beta(i,j)
      enddo
    enddo
    if (Grid%f_plane) then
      write (stdoutunit,'(//,a,f6.2,a//)') &
      ' Note: "f plane" set using f0 at latitude =', Grid%f_plane_latitude,'deg'
    else
      write (stdoutunit,'(//,a,f6.2,a,g14.7//)') &
      ' Note: "beta plane" set using f0 at latitude=',&
       Grid%f_plane_latitude,'deg and beta =',Grid%beta(isc,jsc)
    endif

  else

    Grid%f(:,:)    = 2.0*omega*sin(Grid%phiu(:,:))
    Grid%fstar(:,:)= 2.0*omega*cos(Grid%phiu(:,:))
    Grid%beta(:,:) = 2.0*omega*cos(Grid%phiu(:,:))/radius

  endif

  ! beta_eff centered on u-cell
  ! beta_eff = H*|grad(f/H)| --> |(f/hu)*dht_dx| + |beta-(f/hu)*dht_dy|
  Grid%beta_eff(:,:) = 0.0
  do j=jsc,jec
    do i=isc,iec
      Grid%beta_eff(i,j) = sqrt( (Grid%f(i,j)/(epsln+Grid%hu(i,j))*Grid%dht_dx(i,j))**2 &
                                +(Grid%beta(i,j)-Grid%f(i,j)/(epsln+Grid%hu(i,j))*Grid%dht_dy(i,j))**2 )
    enddo
  enddo 

  ! check for marginally resolved inertial oscillation
  fmax = 2.0*omega*0.0002 ! ~ 0.01 deg latitude
  do j=jsc,jec
    do i=isc,iec
      if (Grid%kmu(i,j) /= 0) then
        fmax = max(fmax,abs(Grid%f(i,j)))
      endif
    enddo
  enddo
  call mpp_max (fmax)
  max_dt_for_inertial_oscillation = int(1.0/fmax) ! assume 2*pi timesteps to resolve oscillation
  write (stdoutunit,'(/1x, a,f8.0,a)')&
  ' ==> Note: 2*pi timesteps/(min inertial period) implies a maximum dtuv =',&
  max_dt_for_inertial_oscillation,' sec.'

  if(acor == 0.0) then 
    write (stdoutunit,'(6x,a)') &
    'Coriolis force is treated explicitly in time since the acor=0.'
  endif 
  if(acor > 0.0 .and. acor < 0.5) then 
    call mpp_error(FATAL,&
    '==> Error in ocean_coriolis_mod: "acor" must be set either to acor=0.0 or 0.5<=acor<=1.0')
  endif 
  if(acor>=0.5 .and. acor<=1.0) then 
    write (stdoutunit,'(6x,a,f6.3)') &
    ' ==> Note: Coriolis force is treated semi-implicitly in time with acor= ',acor
    write (stdoutunit,'(6x,a)') &
    'This approach removes the constraint on dtuv associated with inertial oscillations.'
  endif 

  if (Time_steps%dtuv > max_dt_for_inertial_oscillation .and. acor==0.0) then
    call mpp_error(FATAL,&
    '==> Error in ocean_coriolis_mod: inertial oscillation not resolved. Reduce "dtuv" or set 0.5<=acor<=1.0.')
  endif

  ! diagnostic manager registers and sends for static fields 

  id_coriolis  = register_static_field ('ocean_model', 'f_coriolis', Grd%vel_axes_uv(1:2), &
                                        'Coriolis frequency on U-cell', '1/s',             &
                                         missing_value=missing_value, range=(/-10.0,10.0/))
  if (id_coriolis > 0) used = send_data (id_coriolis, Grid%f(:,:), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  id_beta  = register_static_field ('ocean_model', 'beta', Grd%vel_axes_uv(1:2), &
                                    'planetary beta', '1/(m*s)',&
                                    missing_value=missing_value, range=(/-10.0,10.0/))
  if (id_beta > 0) used = send_data (id_beta, Grd%beta(:,:), &
                          Time%model_time, rmask=Grd%umask(:,:,1), &
                          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  id_beta_eff  = register_static_field ('ocean_model', 'beta_eff', Grd%vel_axes_uv(1:2), &
                                        'effective beta', '1/(m*s)',&
                                        missing_value=missing_value, range=(/-10.0,10.0/))
  if (id_beta_eff > 0) used = send_data (id_beta_eff, Grd%beta_eff(:,:), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! diagnostic manager registers for dynamic fields 
  id_cor_u =  register_diag_field ('ocean_model', 'cor_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'explicit coriolis accel in i-direct', 'm/s^2', missing_value=missing_value, range=(/-1e9,1e9/))
  id_cor_v =  register_diag_field ('ocean_model', 'cor_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'explicit coriolis accel in j-direct', 'm/s^2', missing_value=missing_value, range=(/-1e9,1e9/))

  id_hrho_cor_u =  register_diag_field ('ocean_model', 'hrho_cor_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'rho*dz weighted explicit coriolis accel in i-direct', '(kg/m^3)*m^2/s^2', &
     missing_value=missing_value, range=(/-1e9,1e9/))
  id_hrho_cor_v =  register_diag_field ('ocean_model', 'hrho_cor_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'rho*dz weighted explicit coriolis accel in j-direct', '(kg/m^3)*m^2/s^2', &
     missing_value=missing_value, range=(/-1e9,1e9/))

  id_ucori_impl = register_diag_field ('ocean_model', 'ucori_impl', Grd%vel_axes_uv(1:3), Time%model_time, &
     'implicit Coriolis force in i-direction', '(kg/m^3)*(m/s^2)', missing_value=missing_value, range=(/-1e9,1e9/))
  id_vcori_impl = register_diag_field ('ocean_model', 'vcori_impl', Grd%vel_axes_uv(1:3), Time%model_time, &
     'implicit Coriolis force in j-direction', '(kg/m^3)*(m/s^2)', missing_value=missing_value, range=(/-1e9,1e9/))

end subroutine ocean_coriolis_init
! </SUBROUTINE> NAME="ocean_coriolis_init"


!#######################################################################
! <SUBROUTINE NAME="coriolis_force">
!
! <DESCRIPTION>
! Compute thickness and density weighted acceleration due to Coriolis
! force on a B-grid.
! </DESCRIPTION>
!
subroutine coriolis_force(Time, Thickness, Velocity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step

  integer :: itime, tau, taum1, taup1
  integer :: i, j, k, n
  
  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_coriolis_mod (coriolis_force): module must be initialized')
  endif 

  if(.not. use_this_module) then 
      if(energy_analysis_step) then 
          Velocity%wrkv=0.0
      endif
      return
  endif

  wrk1_v = 0.0

  tau    = Time%tau
  taum1  = Time%taum1
  taup1  = Time%taup1

  if(acor == 0.0) then 
    itime = tau
  else 
    itime = taum1
  endif 


  if(energy_analysis_step) then 

      if(acor > 0) then 

          ! Coriolis force alters kinetic energy when acor>0.  To determine 
          ! effects via the energy analysis, compute Coriolis force as here.  
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   Velocity%wrkv(i,j,k,1) = Grd%f(i,j)  &
                        *(Velocity%u(i,j,k,2,taum1)*(1.0-acor) + acor*Velocity%u(i,j,k,2,taup1)) &
                        *Thickness%rho_dzu(i,j,k,tau)
                   Velocity%wrkv(i,j,k,2) =-Grd%f(i,j) &
                        *(Velocity%u(i,j,k,1,taum1)*(1.0-acor) + acor*Velocity%u(i,j,k,1,taup1)) &
                        *Thickness%rho_dzu(i,j,k,tau)
                enddo
             enddo
          enddo

      else 

          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   Velocity%wrkv(i,j,k,1) =  Grd%f(i,j)*Velocity%u(i,j,k,2,itime)*Thickness%rho_dzu(i,j,k,tau)
                   Velocity%wrkv(i,j,k,2) = -Grd%f(i,j)*Velocity%u(i,j,k,1,itime)*Thickness%rho_dzu(i,j,k,tau)
                enddo
             enddo
          enddo
      endif

  else 

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) =  Grd%f(i,j)*Velocity%u(i,j,k,2,itime)
               wrk1_v(i,j,k,2) = -Grd%f(i,j)*Velocity%u(i,j,k,1,itime)
            enddo
         enddo
      enddo

      if (id_cor_u > 0) used = send_data(id_cor_u, wrk1_v(isc:iec,jsc:jec,:,1), &
                               Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
      if (id_cor_v > 0) used = send_data(id_cor_v, wrk1_v(isc:iec,jsc:jec,:,2), &
                               Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))

      ! weight acceleration with thickness and density of velocity cell
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v(i,j,k,n)         = wrk1_v(i,j,k,n)*Thickness%rho_dzu(i,j,k,tau)
                  Velocity%accel(i,j,k,n) = Velocity%accel(i,j,k,n) + wrk1_v(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

      if (id_hrho_cor_u > 0) used = send_data( id_hrho_cor_u, wrk1_v(isc:iec,jsc:jec,:,1), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
      if (id_hrho_cor_v > 0) used = send_data( id_hrho_cor_v, wrk1_v(isc:iec,jsc:jec,:,2), &
           Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif


end subroutine coriolis_force
! </SUBROUTINE> NAME="coriolis_force"


!#######################################################################
! <SUBROUTINE NAME="coriolis_force_implicit">
!
! <DESCRIPTION>
! Contributions to thickness weighted and density weighted 
! acceleration from time-implicit Coriolis force.
! </DESCRIPTION>
!
subroutine coriolis_force_implicit(Time, Velocity)

  type(ocean_time_type),     intent(in)    :: Time
  type(ocean_velocity_type), intent(inout) :: Velocity

  integer :: i, j, k
  real    :: dtimeacor
  real    :: lambda,factor

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) then 
      return
  endif


  wrk1_v = 0.0
  
  if(acor > 0.0) then 

      dtimeacor = dtime*acor 

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,:) = Velocity%accel(i,j,k,:)  
               lambda          = dtimeacor*Grd%f(i,j)
               factor          = 1.0/(1.0 + lambda*lambda) 
               wrk1(i,j,k)     = (Velocity%accel(i,j,k,1) + lambda*Velocity%accel(i,j,k,2))*factor
               wrk2(i,j,k)     = (Velocity%accel(i,j,k,2) - lambda*Velocity%accel(i,j,k,1))*factor
            enddo
         enddo
      enddo

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec       
               Velocity%accel(i,j,k,1) = wrk1(i,j,k) 
               Velocity%accel(i,j,k,2) = wrk2(i,j,k) 
               wrk1_v(i,j,k,:)         = Velocity%accel(i,j,k,:) - wrk1_v(i,j,k,:) 
            enddo
         enddo
      enddo

      if (debug_this_module) then
          write(stdoutunit,*) ' ' 
          write(stdoutunit,*) 'From ocean_coriolis_mod: chksums after acor>0 Coriolis' 
          call write_timestamp(Time%model_time)
          write(stdoutunit,*) 'accel(1) = ',mpp_chksum(Velocity%accel(isc:iec,jsc:jec,:,1))
          write(stdoutunit,*) 'accel(2) = ',mpp_chksum(Velocity%accel(isc:iec,jsc:jec,:,2))
          write(stdoutunit,*) 'vel(1)   = ',mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,1))
          write(stdoutunit,*) 'vel(2)   = ',mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,2))
      endif

  endif

  ! send to diagnostics manager  
  if (id_ucori_impl > 0) used = send_data(id_ucori_impl, wrk1_v(:,:,:,1), &
                                Time%model_time, rmask=Grd%umask(:,:,:),  &
                                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)  

  if (id_vcori_impl > 0) used = send_data(id_vcori_impl, wrk1_v(:,:,:,2), &
                                Time%model_time, rmask=Grd%umask(:,:,:),  &
                                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)  


end subroutine coriolis_force_implicit
! </SUBROUTINE> NAME="coriolis_force_implicit"




end module ocean_coriolis_mod
