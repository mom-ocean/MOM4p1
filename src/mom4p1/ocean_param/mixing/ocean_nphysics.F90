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
module ocean_nphysics_mod
! 
!<CONTACT EMAIL="stephen.griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Driver for ocean neutral physics.
!</OVERVIEW>
!
!<DESCRIPTION>
! Driver for ocean neutral physics.
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_nphysics_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. 
!  Default use_this_module=.false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  Default debug_this_module=.false.
!  </DATA> 
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="use_nphysicsA" TYPE="logical">
!  For using the nphysicsA method of neutral physics, based on that 
!  developed in MOM4p0.  This scheme is more robust and recommended for 
!  general use.  Default use_nphysicsA=.true. 
!  </DATA> 
!  <DATA NAME="use_nphysicsB" TYPE="logical">
!  For using the nphysicsB method of neutral physics.  This method is 
!  experimental, and is not recommended for general use.
!  Default use_nphysicsB=.false. 
!  </DATA> 
!  <DATA NAME="use_nphysicsC" TYPE="logical">
!  For using the nphysicsC method of neutral physics.  This method is 
!  experimental, and is not recommended for general use.  
!  Default use_nphysicsC=.false. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln, pi, grav, rho0r
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: FATAL, WARNING, NOTE
use fms_mod,          only: file_exist
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_mod,          only: mpp_error, stdout, stdlog
use mpp_domains_mod,  only: mpp_update_domains

use ocean_domains_mod,       only: get_local_indices
use ocean_nphysics_util_mod, only: ocean_nphysics_util_init
use ocean_nphysicsA_mod,     only: ocean_nphysicsA_init, ocean_nphysicsA_end, nphysicsA
use ocean_nphysicsA_mod,     only: ocean_nphysicsA_restart
use ocean_nphysicsB_mod,     only: ocean_nphysicsB_init, ocean_nphysicsB_end, nphysicsB
use ocean_nphysicsB_mod,     only: ocean_nphysicsB_restart
use ocean_nphysicsC_mod,     only: ocean_nphysicsC_init, ocean_nphysicsC_end, nphysicsC
use ocean_nphysicsC_mod,     only: ocean_nphysicsC_restart
use ocean_parameters_mod,    only: TERRAIN_FOLLOWING, missing_value, onefourth
use ocean_types_mod,         only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_thickness_type, ocean_density_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type, ocean_options_type


implicit none

public ocean_nphysics_init
public ocean_nphysics_end
public neutral_physics 
public ocean_nphysics_restart

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

logical :: used

logical :: use_nphysicsA     = .true.   
logical :: use_nphysicsB     = .false.   
logical :: use_nphysicsC     = .false.   
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.
logical :: write_a_restart   = .true. 

#include <ocean_memory.h>

character(len=128) :: version=&
     '$Id: ocean_nphysics.F90,v 1.1.2.3.2.30.22.1.38.1.54.1 2009/10/10 00:42:24 nnz Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: module_is_initialized = .FALSE.

namelist /ocean_nphysics_nml/ use_this_module, debug_this_module, write_a_restart, &
                              use_nphysicsA, use_nphysicsB, use_nphysicsC

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module. 
! </DESCRIPTION>
!
subroutine ocean_nphysics_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,    &
                               Ocean_options, vert_coordinate_type, vert_coordinate_class, &
                               debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_thickness_type),   intent(in)           :: Thickness
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_options_type),     intent(inout)        :: Ocean_options
  integer,                      intent(in)           :: vert_coordinate_type
  integer,                      intent(in)           :: vert_coordinate_class
  logical,                      intent(in), optional :: debug

  integer :: ioun, io_status, ierr  
  integer :: num_schemes 
  real    :: agm_closure_lower_depth
  real    :: agm_closure_upper_depth
  real    :: agm_closure_buoy_freq
  real    :: smax 
  real    :: swidth 
 
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  
  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysics_mod (ocean_nphysics_init):already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  Dom => Domain
  Grd => Grid

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_nphysics_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_nphysics_nml)  
  write (stdlogunit,ocean_nphysics_nml)
  ierr = check_nml_error(io_status,'ocean_nphysics_nml')
  call close_file (ioun)

  if(use_this_module) then 
    call mpp_error(NOTE, &
    '==> from ocean_nphysics_mod: USING ocean_nphysics_mod.')
    write(stdoutunit,'(1x,a)')    &
    '==> Note from ocean_nphysics_mod: USING ocean_nphysics.'
    if(vert_coordinate_type==TERRAIN_FOLLOWING) then 
    call mpp_error(WARNING, &
    '==>Warning: ocean_nphysics is NOT supported with TERRRAIN_FOLLOWING vertical coordinates.')
    endif 
  else 
    call mpp_error(NOTE, &
    '==> from ocean_nphysics_mod: NOT using ocean_nphysics_mod.')
    write(stdoutunit,'(1x,a)')    &
    '==> Note from ocean_nphysics_mod: NOT using ocean_nphysics.'
    Ocean_options%neutral_physics = 'Did NOT use neutral physics option.'
    return
  endif 

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_nphysics_mod with debug_this_module=.true.'  
  endif

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_nphysics_mod: NO restart written.'
    call mpp_error(WARNING,&
    '==>Warning from ocean_nphysics_mod: NO restart written.')
  endif 

#ifndef MOM4_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  call ocean_nphysics_util_init(Grid, Domain, Time, Time_steps, Dens, T_prog,  &
       agm_closure_lower_depth, agm_closure_upper_depth, agm_closure_buoy_freq,&
       smax, swidth, debug)

  num_schemes=0
  if(use_nphysicsA) then 
    num_schemes = num_schemes+1
    call mpp_error(NOTE, &
    '==> from ocean_nphysics_mod: USING ocean_nphysicsA.')
    write(stdoutunit,'(1x,a)') &
    '==> Note from ocean_nphysics_mod: Using ocean_nphysicsA.'
    Ocean_options%neutral_physics = 'Used neutral physics using nphysicsA algorithm.'
    call ocean_nphysicsA_init(Grid, Domain, Time, Time_steps, Thickness, T_prog, &
         vert_coordinate_class, agm_closure_lower_depth, agm_closure_upper_depth,&
         smax, swidth, debug)

  elseif(use_nphysicsB) then 
    num_schemes = num_schemes+1
    call mpp_error(NOTE, &
    '==> from ocean_nphysics_mod: USING ocean_nphysicsB.')
    write(stdoutunit,'(1x,a)') &
    '==> Note from ocean_nphysics_mod: Using ocean_nphysicsB.'
    write(stdoutunit,'(1x,a)') &
    '    the ocean_nphysicsB module is experimental, and not recommended for general use.'
    Ocean_options%neutral_physics = 'Used neutral physics using nphysicsB algorithm.'
    call ocean_nphysicsB_init(Grid, Domain, Time, Time_steps, Thickness, T_prog, &
         vert_coordinate_class, agm_closure_lower_depth, agm_closure_upper_depth,&
         smax, swidth, debug)

  elseif(use_nphysicsC) then 
    num_schemes = num_schemes+1
    call mpp_error(NOTE, &
    '==> from ocean_nphysics_mod: USING ocean_nphysicsC.')
    write(stdoutunit,'(1x,a)') &
    '==> Note from ocean_nphysics_mod: Using ocean_nphysicsC.'
    write(stdoutunit,'(1x,a)') &
    '    the ocean_nphysicsC module is experimental, and not recommended for general use.'
    Ocean_options%neutral_physics = 'Used neutral physics using nphysicsC algorithm.'
    call ocean_nphysicsC_init(Grid, Domain, Time, Time_steps, Thickness, T_prog, &
         vert_coordinate_class, agm_closure_lower_depth, agm_closure_upper_depth,&
         smax, swidth, debug)
  endif 

  if(num_schemes > 1) then 
     call mpp_error(FATAL,  &
    '==>ocean_nphysics_mod: Can only enable one of the nphysics schemes: A, B, or C.')
     write(stdoutunit,'(1x,a)') &
    '==>ocean_nphysics_mod: Can only enable one of the nphysics schemes: A, B, or C.'
  endif 
  if(num_schemes==0) then 
     call mpp_error(WARNING,  &
    '==>ocean_nphysics_mod: no nphysics scheme enabled. Choose one of nphysicsA, nphysicsB, or nphysicsC.')
     write(stdoutunit,'(1x,a)') &
    '==>ocean_nphysics_mod: no nphysics scheme enabled. Choose one of nphysicsA, nphysicsB, or nphysicsC.'
  endif 


end subroutine ocean_nphysics_init
! </SUBROUTINE>  NAME="ocean_nphysics_init"


!#######################################################################
! <SUBROUTINE NAME="neutral_physics">
! <DESCRIPTION>
!
! Call the relevant neutral physics scheme.
!
! </DESCRIPTION>
!
!
subroutine neutral_physics (Time, Thickness, Dens, rho, T_prog, &
               gm_diffusivity, surf_blthick, bott_blthick, rossby_radius_raw)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  real, dimension(isd:,jsd:,:),   intent(in)    :: rho
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(in)    :: surf_blthick
  real, dimension(isd:,jsd:),     intent(in)    :: bott_blthick
  real, dimension(isd:,jsd:,:),   intent(inout) :: gm_diffusivity
  real, dimension(isd:,jsd:),     intent(inout) :: rossby_radius_raw

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysics (neutral_physics): needs initialization')
  endif 

  if(use_nphysicsA) then 
      call nphysicsA(Time, Thickness, Dens, rho, T_prog, &
           surf_blthick, gm_diffusivity, rossby_radius_raw)
  elseif(use_nphysicsB) then 
      call nphysicsB(Time, Thickness, Dens, rho, T_prog, &
           surf_blthick, gm_diffusivity, rossby_radius_raw)
  elseif(use_nphysicsC) then 
      call nphysicsC(Time, Thickness, Dens, rho, T_prog, &
           surf_blthick, gm_diffusivity, rossby_radius_raw)
  endif

end subroutine neutral_physics
! </SUBROUTINE> NAME="neutral_physics"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_restart">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysics_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  if(.not. use_this_module) return

  if(use_nphysicsA) then 
     call ocean_nphysicsA_restart(time_stamp)
  endif 
  if(use_nphysicsB) then 
     call ocean_nphysicsB_restart(time_stamp)
  endif 
  if(use_nphysicsC) then 
     call ocean_nphysicsC_restart(time_stamp)
  endif   

end subroutine ocean_nphysics_restart
! </SUBROUTINE> NAME="ocean_nphysics_restart"



!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysics_end(Time)

  type(ocean_time_type), intent(in) :: Time
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysics (ocean_nphysics_end): needs initialization')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_nphysics_mod: NO restart written.'
    call mpp_error(WARNING,&
    '==>Warning from ocean_nphysics_mod: NO restart written.')
    return 
  endif 

  call ocean_nphysics_restart

  if(use_nphysicsA) then 
     call ocean_nphysicsA_end(Time)
  endif 
  if(use_nphysicsB) then 
     call ocean_nphysicsB_end(Time)
  endif 
  if(use_nphysicsC) then 
     call ocean_nphysicsC_end(Time)
  endif 

end subroutine ocean_nphysics_end
! </SUBROUTINE> NAME="ocean_nphysics_end"


end module ocean_nphysics_mod
