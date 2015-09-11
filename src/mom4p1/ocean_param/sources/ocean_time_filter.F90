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
module ocean_time_filter_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<OVERVIEW>
! Perform Robert-Asselin time filter when using leap-frog based
! time stepping scheme along with GEOPOTENTIAL vertical coordinate.
!</OVERVIEW>
!
!<DESCRIPTION>
! Module to perform Robert-Asselin time filter on prognostic fields.
! This time filter is only applicable when using the older leap-frog
! based time stepping scheme AND when vertical coordinate is GEOPOTENTIAL. 
! The leap-frog time stepping scheme is not recommended. 
! It remaines in mom4 for legacy purposes only. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_time_filter_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module.  Default is false. 
!  </DATA> 
!
!  <DATA NAME="robert_asselin_param" TYPE="real">
!  Parameter setting the strength of the time filtering applied via the 
!  Robert-Asselin filter when using threelevel_time_filter=.true.
!  Typically taken as robert ~ 0.05 for global models. Notably, however,
!  this value may not be sufficient to suppress time step noise in regions
!  of high frequency variability, such as along the equator.  
!  </DATA> 
!
!  <DATA NAME="velocity_change_diag_step" TYPE="integer">
!  perform velocity_change_check every n timesteps (1 == every_tstep). 
!  </DATA> 
!  <DATA NAME="velocity_change_check" TYPE="logical">
!  For checking to see whether the abs(vel(tau)-vel(tam1))
!  is greater than velocity_delta_max.  If so, then we likely have 
!  problems with leap-frog noise.  
!  </DATA> 
!  <DATA NAME="velocity_change_max" TYPE="real" UNITS="meter/sec" >
!  For checking to see whether the abs(vel(tau)-vel(taum1)) >velocity_change_max.
!  If so, then will bring the model down if have problems in too many places.  
!  </DATA> 
!  <DATA NAME="velocity_change_max_num" TYPE="integer">
!  Maximum number of points allowed where abs(vel(tau)-vel(taum1)) >velocity_change_max.
!  </DATA> 

!</NAMELIST>

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_domains_mod,  only: mpp_update_domains
use mpp_mod,          only: mpp_error, FATAL, NOTE, stdout, stdlog

use ocean_domains_mod,       only: get_local_indices
use ocean_operators_mod,     only: REMAP_BT_TO_BU
use ocean_parameters_mod,    only: GEOPOTENTIAL, missing_value, rho0
use ocean_types_mod,         only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,         only: ocean_velocity_type, ocean_thickness_type
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_external_mode_type 
use ocean_velocity_diag_mod, only: velocity_change
use ocean_workspace_mod,     only: wrk1, wrk1_v

implicit none

private

public time_filter 
public ocean_time_filter_init

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk
integer :: num_prog_tracers

character(len=128) :: version = &
     '$Id: ocean_time_filter.F90,v 16.0.108.1 2009/10/10 00:42:55 nnz Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

! time step related variables 
real :: p5dttsr  ! 1/(2dtts) 
real :: p5dtuvr  ! 1/(2dtuv)


! for diagnosing leap-frog noise in velocity field (look at the equator for amazingly large noise levels)
logical :: velocity_change_check     =.false. ! for checking abs of the single time step change in velocity
real    :: velocity_change_max       = 10.0   ! (m/s)
integer :: velocity_change_max_num   = 10     ! max times permit abs(vel(tau)-vel(taum1)) > velocity_change_max
integer :: velocity_change_diag_step = -1     ! velocity_change_check every n timesteps (1==every_tstep)


! for diagnostics sent to diagnostic manager 
integer, allocatable, dimension(:) :: id_tfilter_tracer
integer :: id_tfilter_velocity(2)=-1
integer :: id_tfilter_etat=-1
integer :: id_tfilter_etau=-1
integer :: id_forward_u
integer :: id_forward_v
logical :: used

real    :: robert_asselin_param  = 0.05
logical :: module_is_initialized =.false.
logical :: use_this_module       =.false.

namelist /ocean_time_filter_nml/ use_this_module, robert_asselin_param,           &
                                 velocity_change_max, velocity_change_max_num,    &
                                 velocity_change_check, velocity_change_diag_step

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_time_filter_init">
!
! <DESCRIPTION>
! Initialize the time filter module. 
! </DESCRIPTION>
!
subroutine ocean_time_filter_init(Grid, Domain, Time, Time_steps, T_prog, vert_coordinate, time_tendency)

    type(ocean_grid_type),        intent(in), target :: Grid
    type(ocean_domain_type),      intent(in), target :: Domain
    type(ocean_time_type),        intent(in)         :: Time
    type(ocean_time_steps_type),  intent(inout)      :: Time_steps
    type(ocean_prog_tracer_type), intent(in)         :: T_prog(:)
    integer,                      intent(in)         :: vert_coordinate
    character(len=32),            intent(in)         :: time_tendency 

    integer :: ioun, io_status, ierr
    integer :: n

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error from ocean_time_filter_mod (ocean_time_filter_init): module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    Grd => Grid
    Dom => Domain

    call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grd%nk

    num_prog_tracers = size(T_prog(:))

    ! provide for namelist over-ride of defaults 
    ioun = open_namelist_file()
    read  (ioun, ocean_time_filter_nml,iostat=io_status)
    write (stdoutunit,'(/)')
    write (stdoutunit, ocean_time_filter_nml)
    write (stdlogunit, ocean_time_filter_nml)
    ierr = check_nml_error(io_status,'ocean_time_filter_nml')
    call close_file (ioun)

    if(use_this_module) then 
      call mpp_error(NOTE,'==>USING ocean_time_filter_mod.')
    else 
      call mpp_error(NOTE,'==>NOT using ocean_time_filter_mod.')
      return
    endif 

    if(robert_asselin_param==0.0 .and. time_tendency=='threelevel') then 
      call mpp_error(FATAL,'==>Note from ocean_time_filter_mod: robert_asselin_param=0.0 is unstable with leap-frog.')
    endif 
    if(robert_asselin_param /= 0.0 .and. time_tendency=='threelevel') then 
      write(stdoutunit,'(/1x,a,1x,f6.2)')'==>Note: using Robert-Asselin time filter with param =',robert_asselin_param
    endif 

    if(time_tendency=='twolevel') then 
       write(stdoutunit,'(/1x,a)') '==>Note: using twolevel time tendency, so Robert-Asselin time filter is not used.'
       robert_asselin_param = 0.0  
    endif 

    if(vert_coordinate /= GEOPOTENTIAL) then 
       write(stdoutunit,'(/1x,a)') '==>Note: not using GEOPOTENTIAL, so Robert-Asselin time filter is not used.'
       robert_asselin_param = 0.0
    endif 

    Time_steps%robert_asselin_param = robert_asselin_param
    p5dttsr = 0.5/(epsln+Time_steps%dtts)
    p5dtuvr = 0.5/(epsln+Time_steps%dtuv)

    if (velocity_change_diag_step == 0) velocity_change_diag_step = 1


    ! register diagnostics 

    id_tfilter_velocity(1) = register_diag_field ('ocean_model', 'tfilter_u', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'u-accel from time filter', 'm/s^2',&
           missing_value=missing_value, range=(/-1e6,1e6/))
    id_tfilter_velocity(2) = register_diag_field ('ocean_model', 'tfilter_v', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'v-accel from time filter', 'm/s^2',&
           missing_value=missing_value, range=(/-1e6,1e6/))

    id_tfilter_etat = register_diag_field ('ocean_model', 'tfilter_etat', &
          Grid%tracer_axes(1:2), Time%model_time, &
          'eta_t change from time filter', 'm',&
           missing_value=missing_value, range=(/-1e6,1e6/))

    id_tfilter_etau = register_diag_field ('ocean_model', 'tfilter_etau', &
          Grid%vel_axes_uv(1:2), Time%model_time, &
          'eta_u change from time filter', 'm',&
           missing_value=missing_value, range=(/-1e6,1e6/))

    id_forward_u = register_diag_field ('ocean_model', 'forward_u', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'u(tau)-u(taum1)', 'm/s',&
           missing_value=missing_value, range=(/-1e6,1e6/))
    id_forward_v = register_diag_field ('ocean_model', 'forward_v', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'v(tau)-v(taum1)', 'm/s',&
           missing_value=missing_value, range=(/-1e6,1e6/))

     allocate( id_tfilter_tracer(num_prog_tracers) )
     id_tfilter_tracer(:) = -1
     do n=1,num_prog_tracers
        if (T_prog(n)%name=='temp') then
           id_tfilter_tracer(n) = register_diag_field ('ocean_model', 'tfilter_'//trim(T_prog(n)%name), &
                Grid%tracer_axes(1:3), Time%model_time, &
                'thk wghtd heating from time filter', 'Watts/m^2',&
                missing_value=missing_value, range=(/-1e6,1e6/))
        else 
           id_tfilter_tracer(n) = register_diag_field ('ocean_model', 'tfilter_'//trim(T_prog(n)%name), &
                Grid%tracer_axes(1:3), Time%model_time,                                                 &
                'thk wghtd tendency from time filter for '//trim(T_prog(n)%longname), 'kg/(m^2*sec)',   &
                missing_value=missing_value, range=(/-1e6,1e6/))
        endif
     enddo


end subroutine ocean_time_filter_init
! </SUBROUTINE> NAME="ocean_time_filter_init"



!#######################################################################
! <SUBROUTINE NAME="time_filter">
!
! <DESCRIPTION>
! When threelevel time tendency is used, perform Robert-Asselin time 
! filtering on velocity and tracer.
!
! Also perform time filter on surface height by replacing eta_t(tau)
! with eta_t_bar(tau).
!
! This scheme is only implemented for vertical_coordinate='geopotential'.
! Other vertical coordinates must use the default TWO_LEVEL tendency.
!
! </DESCRIPTION>
!
  subroutine time_filter(Time, Grid, Thickness, T_prog, Velocity, Ext_mode)

    type(ocean_time_type),          intent(in)    :: Time
    type(ocean_grid_type),          intent(in)    :: Grid
    type(ocean_thickness_type),     intent(inout) :: Thickness
    type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
    type(ocean_velocity_type),      intent(inout) :: Velocity
    type(ocean_external_mode_type), intent(inout) :: Ext_mode

    real, dimension(isd:ied,jsd:jed) :: tmp              
    integer                          :: i, j, k, n
    integer                          :: tau, taum1, taup1 

    if(.not. use_this_module) return 

    taum1 = Time%taum1
    tau   = Time%tau
    taup1 = Time%taup1    
 

!---surface height and thickness of surface cell 

    tmp=0.0
    if(id_tfilter_etat > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp(i,j) = Ext_mode%eta_t_bar(i,j,tau) - Ext_mode%eta_t(i,j,tau)
           enddo
        enddo
        used = send_data(id_tfilter_etat, tmp(:,:), &
               Time%model_time, rmask=Grid%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif
    if(id_tfilter_etau > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp(i,j) = Ext_mode%eta_u(i,j,tau)
           enddo
        enddo
    endif

    Ext_mode%eta_t(:,:,tau) = Ext_mode%eta_t_bar(:,:,tau)
    Ext_mode%eta_u(:,:,tau) = Grid%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%eta_t(:,:,tau))
    call mpp_update_domains (Ext_mode%eta_u(:,:,tau), Dom%domain2d)

    Thickness%rho_dzt(:,:,1,tau) = rho0*(Grid%dzt(1) + Ext_mode%eta_t(:,:,tau))
    Thickness%rho_dzu(:,:,1,tau) = rho0*(Grid%dzt(1) + Ext_mode%eta_u(:,:,tau))

    if(id_tfilter_etau > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp(i,j) = Ext_mode%eta_u(i,j,tau) - tmp(i,j)
           enddo
        enddo
        used = send_data(id_tfilter_etau, tmp(:,:), &
               Time%model_time, rmask=Grid%umask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif


!---tracers 

    wrk1 = 0.0
    do n=1,num_prog_tracers

       if(id_tfilter_tracer(n) > 0) then 
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1(i,j,k) = T_prog(n)%field(i,j,k,tau)
                 enddo
              enddo
           enddo
       endif

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                T_prog(n)%field(i,j,k,tau) = T_prog(n)%field(i,j,k,tau) + &
                robert_asselin_param*(  &
                0.5*(T_prog(n)%field(i,j,k,taup1) + T_prog(n)%field(i,j,k,taum1)) - T_prog(n)%field(i,j,k,tau))
             enddo
          enddo
       enddo

       if(id_tfilter_tracer(n) > 0) then 
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1(i,j,k) = (T_prog(n)%field(i,j,k,tau)-wrk1(i,j,k)) &
                                   *p5dttsr*Thickness%rho_dzt(i,j,k,tau)*T_prog(n)%conversion 
                 enddo
              enddo
           enddo
           used = send_data(id_tfilter_tracer(n), wrk1(:,:,:), &
                  Time%model_time, rmask=Grid%tmask(:,:,:), &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
       endif

    enddo

!---velocity 

    wrk1_v=0.0 

    if(id_tfilter_velocity(1) > 0 .or. id_tfilter_velocity(2) > 0) then 
        do n=1,2
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1_v(i,j,k,n) = Velocity%u(i,j,k,n,tau)
                 enddo
              enddo
           enddo
        enddo
    endif

    do n=1,2
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                Velocity%u(i,j,k,n,tau) = Velocity%u(i,j,k,n,tau) + &
                robert_asselin_param*( &
                0.5*(Velocity%u(i,j,k,n,taup1) + Velocity%u(i,j,k,n,taum1)) - Velocity%u(i,j,k,n,tau))
             enddo
          enddo
       enddo
    enddo

    if(id_tfilter_velocity(1) > 0 .or. id_tfilter_velocity(2) > 0) then 
        do n=1,2
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1_v(i,j,k,n) = p5dtuvr*(Velocity%u(i,j,k,n,tau)-wrk1_v(i,j,k,n))
                 enddo
              enddo
           enddo
        enddo
        if(id_tfilter_velocity(1) > 0) then 
           used = send_data(id_tfilter_velocity(1), wrk1_v(:,:,:,1), &
                  Time%model_time, rmask=Grid%umask(:,:,:), &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
        endif 
        if(id_tfilter_velocity(2) > 0) then 
           used = send_data(id_tfilter_velocity(2), wrk1_v(:,:,:,2), &
                  Time%model_time, rmask=Grid%umask(:,:,:), &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
        endif 
    endif

    ! check leap-frog noise remaining after time filter applied 
    if(id_forward_u > 0) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1_v(i,j,k,1) = Velocity%u(i,j,k,1,tau)-Velocity%u(i,j,k,1,taum1)
              enddo
           enddo
        enddo
        used = send_data(id_forward_u, wrk1_v(:,:,:,1), &
               Time%model_time, rmask=Grid%umask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
    
    if(id_forward_v > 0 ) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1_v(i,j,k,2) = Velocity%u(i,j,k,2,tau)-Velocity%u(i,j,k,2,taum1)
              enddo
           enddo
        enddo
        used = send_data(id_forward_v, wrk1_v(:,:,:,2), &
               Time%model_time, rmask=Grid%umask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif

    if (velocity_change_check .and. velocity_change_diag_step > 0 .and. &
       mod(Time%itt, velocity_change_diag_step) == 0) then 
      call velocity_change(Time, Velocity, velocity_change_max, velocity_change_max_num)
    endif 


  end subroutine time_filter
! </SUBROUTINE> NAME="time_filter"


end module ocean_time_filter_mod
