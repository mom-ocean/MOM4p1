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
module ocean_lap_friction_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! This module calls the appropriate lateral laplacian friction modules. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module serves as an interface to the chosen lateral 
! laplacian friction modules.  
!</DESCRIPTION>
!
! <INFO>
!
! <NOTE>
! The model can generally run with both Laplacian and biharmonic friction
! enabled at the same time.  Such has been found useful for some eddying 
! ocean simulations. 
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_lap_friction_nml">
!
! <DATA NAME="lap_friction_scheme" TYPE="character">
! To determine the laplacian friction scheme: "const" or "general"
! </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging.  
!  </DATA> 
!
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!</NAMELIST>

use diag_manager_mod,only: register_diag_field, send_data
use fms_mod,         only: open_namelist_file, check_nml_error, write_version_number, close_file
use fms_mod,         only: FATAL, WARNING, stdout, stdlog
use fms_mod,         only: file_exist
use fms_io_mod,      only: register_restart_field, save_restart, restore_state
use fms_io_mod,      only: restart_file_type
use mpp_domains_mod, only: mpp_update_domains
use mpp_mod,         only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE
use mpp_mod,         only: mpp_error, mpp_chksum

use ocean_domains_mod,         only: get_local_indices
use ocean_lapcst_friction_mod, only: ocean_lapcst_friction_init, lapcst_friction
use ocean_lapcst_friction_mod, only: lapcst_viscosity_check, lapcst_reynolds_check
use ocean_lapgen_friction_mod, only: ocean_lapgen_friction_init, lapgen_friction
use ocean_lapgen_friction_mod, only: lapgen_viscosity_check, lapgen_reynolds_check 
use ocean_operators_mod,       only: BAX, BAY, BDX_EU, BDY_NU, FDX_U, FDY_U, FMX, FMY
use ocean_parameters_mod,      only: missing_value
use ocean_types_mod,           only: ocean_time_type, ocean_grid_type
use ocean_types_mod,           only: ocean_domain_type, ocean_adv_vel_type
use ocean_types_mod,           only: ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,           only: ocean_options_type 
use ocean_util_mod,            only: write_timestamp

implicit none

private

public ocean_lap_friction_init
public ocean_lap_friction_end
public lap_friction 
public lap_viscosity_check
public lap_reynolds_check
public lap_friction_barotropic
public ocean_lap_friction_restart

! for choosing the mixing scheme
character(len=10) :: lap_friction_scheme='general'  ! 'const' or 'general'
integer :: MIX_SCHEME=2  ! 1=const, 2=general
logical :: debug_this_module = .false.

! for deciding whether we are using laplacian friction 
logical :: use_lapcst_friction=.false. 
logical :: use_lapgen_friction=.false. 
logical :: use_this_module    =.false. 

! vertically averaged viscosity on east and north face of U cells (m^2/s)
real, dimension(:,:), allocatable :: lap_viscosity
real, dimension(:,:), allocatable :: lap_visc_ceu
real, dimension(:,:), allocatable :: lap_visc_cnu

! for diagnostics 
integer :: id_lap_viscosity=-1 
logical :: used

! for clocks
integer :: id_clock_lapcst
integer :: id_clock_lapgen

! for restart
type(restart_file_type), save :: Lap_restart

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
     '=>Using: ocean_lap_friction.f90 ($Id: ocean_lap_friction.F90,v 16.0.2.1.38.1.54.1 2009/10/10 00:42:19 nnz Exp $)'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: module_is_initialized = .FALSE.
logical :: write_a_restart       = .true.

namelist /ocean_lap_friction_nml/ lap_friction_scheme, debug_this_module, write_a_restart 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_lap_friction_init">
!
! <DESCRIPTION>
! Initialize the horizontal Laplacian friction module.
! </DESCRIPTION>
!
subroutine ocean_lap_friction_init(Grid, Domain, Time, Ocean_options, dtime, obc, debug)

  type(ocean_grid_type),   target,  intent(in)    :: Grid
  type(ocean_domain_type), target,  intent(in)    :: Domain
  type(ocean_time_type),            intent(in)    :: Time
  type(ocean_options_type),         intent(inout) :: Ocean_options
  real,                             intent(in)    :: dtime
  logical,                          intent(in)    :: obc
  logical,                optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: num_methods=0 
  integer :: id_restart

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
     '==>Error in ocean_lap_friction_mod (ocean_lap_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_lap_friction_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_lap_friction_nml)  
  write (stdlogunit,ocean_lap_friction_nml)
  ierr = check_nml_error(io_status,'ocean_lap_friction_nml')
  call close_file(ioun)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  Grd => Grid
  Dom => Domain

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 

  ! for clocks 
  id_clock_lapcst = mpp_clock_id('(Ocean const   lap frict) ' ,grain=CLOCK_MODULE)
  id_clock_lapgen = mpp_clock_id('(Ocean general lap frict) ' ,grain=CLOCK_MODULE)

  ! vertically averaged vicsosities 
  allocate (lap_viscosity(isd:ied,jsd:jed))
  allocate (lap_visc_ceu(isd:ied,jsd:jed))
  allocate (lap_visc_cnu(isd:ied,jsd:jed))
  lap_viscosity(:,:) = 0.0
  lap_visc_ceu(:,:)  = 0.0
  lap_visc_cnu(:,:)  = 0.0


  ! initialize laplacian friction scheme
 if(lap_friction_scheme == 'const') then 
     num_methods=num_methods+1
     write(stdoutunit,'(a)') &
     '==>Note from ocean_lap_friction_init: constant laplacian friction scheme used.'
     MIX_SCHEME=1
     call ocean_lapcst_friction_init (Grid, Domain, Time, Ocean_options, dtime, &
                                      obc, use_lapcst_friction, debug_this_module)
 elseif(lap_friction_scheme == 'general') then 
     num_methods=num_methods+1
     write(stdoutunit,'(a)') &
     '==>Note from ocean_lap_friction_init: general laplacian friction scheme used.'
     MIX_SCHEME=2
     call ocean_lapgen_friction_init (Grid, Domain, Time, Ocean_options, dtime, &
                                      obc, use_lapgen_friction, debug_this_module)
 endif

 if(num_methods==0) then 
     write(stdoutunit,'(a)') &
     '==>Error: ocean_lap_friction_mod: pick a laplacian friction scheme: "const" or "general".'
     call mpp_error(FATAL, &
     '==>Error: ocean_lap_friction_mod: pick a laplacian friction scheme: "const" or "general".')
 elseif(num_methods>1) then 
     write(stdoutunit,'(a)') &
     '==>Error: ocean_lap_friction_mod: choose only one laplacian friction scheme.'
     call mpp_error(FATAL, &
     '==>Error: ocean_lap_friction_mod: choose only one laplacian friction scheme.')
 endif

  ! to determine whether we are doing any laplacian friction 
  if(use_lapcst_friction .or. use_lapgen_friction) then 
     use_this_module = .true. 
  endif 


  if(use_this_module) then 
       id_restart = register_restart_field(Lap_restart, 'ocean_lap_friction.res.nc', 'lap_viscosity', &
                    lap_viscosity, domain=Dom%domain2d)
      if(.NOT.file_exist('INPUT/ocean_lap_friction.res.nc')) then
          if (.NOT. Time%init) then 
              call mpp_error(FATAL,'Expecting file INPUT/ocean_lap_friction.res.nc to exist.&
                   &This file was not found and Time%init=.false.')
          endif
      else
          call restore_state(Lap_restart)
          call mpp_update_domains(lap_viscosity,Dom%domain2d)
          write (stdoutunit,'(1x,a)') 'lap_viscosity being read from restart.'
          write (stdoutunit, *) 'checksum start lap_viscosity = ', &
               mpp_chksum(lap_viscosity(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,1))

          lap_visc_ceu(:,:) = BAY(lap_viscosity(:,:))
          lap_visc_cnu(:,:) = BAX(lap_viscosity(:,:))
          call mpp_update_domains (lap_visc_ceu, Dom%domain2d)
          call mpp_update_domains (lap_visc_cnu, Dom%domain2d)

      endif

      if(.not. write_a_restart) then 
          write(stdoutunit,'(a)') '==>Note: running ocean_lap_friction with write_a_restart=.false.'
          write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
      endif

  endif


 id_lap_viscosity = register_diag_field ('ocean_model', 'lap_viscosity_2d', Grd%vel_axes_uv(1:2), &
                     Time%model_time,'Vertically averaged laplacian viscosity',                   &
                    'm^2/s',missing_value=missing_value, range=(/-10.0,1.e20/))


end subroutine ocean_lap_friction_init
! </SUBROUTINE>  NAME="ocean_lap_friction_init"


!#######################################################################
! <SUBROUTINE NAME="lap_friction">
!
! <DESCRIPTION>
! Compute the thickness weighted and density weighted accel due to 
! lateral laplacian friction.  Add this contribution to Velocity%accel. 
! </DESCRIPTION>
!
subroutine lap_friction(Time, Thickness, Adv_vel, Velocity, energy_analysis_step)

  type(ocean_time_type),       intent(in)    :: Time
  type(ocean_thickness_type),  intent(in)    :: Thickness
  type(ocean_adv_vel_type),    intent(in)    :: Adv_vel
  type(ocean_velocity_type),   intent(inout) :: Velocity
  logical,                     intent(in)    :: energy_analysis_step 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lap_friction_mod (lap_friction): module needs to be initialized')
  endif 

  if(MIX_SCHEME==1) then 
    call mpp_clock_begin(id_clock_lapcst)
    call lapcst_friction(Time, Thickness, Velocity, &
                         lap_viscosity, energy_analysis_step)
    call mpp_clock_end(id_clock_lapcst)
  elseif(MIX_SCHEME==2) then 
    call mpp_clock_begin(id_clock_lapgen)
    call lapgen_friction(Time, Thickness, Adv_vel, Velocity, &
                         lap_viscosity, energy_analysis_step)
    call mpp_clock_end(id_clock_lapgen)
  endif 

  ! vertically averaged viscosity for north and east face of U-cells 
  if(use_this_module) then 
     lap_visc_ceu(:,:) = BAY(lap_viscosity(:,:))
     lap_visc_cnu(:,:) = BAX(lap_viscosity(:,:))
     call mpp_update_domains (lap_visc_ceu, Dom%domain2d)
     call mpp_update_domains (lap_visc_cnu, Dom%domain2d)
  endif 

  if (id_lap_viscosity > 0) then 
      used = send_data(id_lap_viscosity, lap_viscosity(:,:), &
             Time%model_time, rmask=Grd%umask(:,:,1),        &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


end subroutine lap_friction
! </SUBROUTINE> NAME="lap_friction"


!#######################################################################
! <SUBROUTINE NAME="lap_viscosity_check">
!
! <DESCRIPTION>
! To check that the viscosity is not too large. 
! </DESCRIPTION>
!
subroutine lap_viscosity_check

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_lap_friction_mod (lap_viscosity_check): needs initialization')
  endif 

  if(MIX_SCHEME==1) then 
    call lapcst_viscosity_check
  elseif(MIX_SCHEME==2) then 
    call lapgen_viscosity_check
  endif 

end subroutine lap_viscosity_check
! </SUBROUTINE> NAME="lap_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="lap_reynolds_check">
!
! <DESCRIPTION>
! To check that the Reynolds number is not too large. 
! </DESCRIPTION>
!
subroutine lap_reynolds_check(Time, Velocity)

  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_lap_friction_mod (lap_reynolds_check): needs initialization')
  endif 

  if(MIX_SCHEME==1) then 
    call lapcst_reynolds_check(Time, Velocity)
  elseif(MIX_SCHEME==2) then 
    call lapgen_reynolds_check(Time, Velocity)
  endif 

end subroutine lap_reynolds_check
! </SUBROUTINE> NAME="lap_reynolds_check"

      
!#######################################################################
! <SUBROUTINE NAME="lap_friction_barotropic">
!
! <DESCRIPTION>
!
! This routine computes the laplacian friction acting on a two-dim
! array.  It uses the two-dimensional vertically averaged viscosity 
! used in the laplacian friction module.  The intent is to apply this
! 2d operator to the vertically integrated horizontal momentum. We 
! ignore the spherical metric terms in this form of the operator, 
! since we are aiming for a fast smoothing operator to be applied
! during each of the many barotropic time steps.  We also apply 
! just the isotropic portion of the more general anisotropic 
! laplacian operator.  
!
! </DESCRIPTION>
!
subroutine lap_friction_barotropic(lap_ceu_back, lap_cnu_back, array, friction)

  real, dimension(isd:,jsd:),   intent(in)    :: lap_ceu_back
  real, dimension(isd:,jsd:),   intent(in)    :: lap_cnu_back 
  real, dimension(isd:,jsd:,:), intent(in)    :: array
  real, dimension(isd:,jsd:,:), intent(inout) :: friction

  integer :: i, j, n
  real    :: tmp(isd:ied,jsd:jed)

  do n=1,2
     tmp(:,:) =  BDX_EU( (lap_ceu_back(:,:)+lap_visc_ceu(:,:))*FDX_U(array(:,:,n)) ) &
                +BDY_NU( (lap_cnu_back(:,:)+lap_visc_cnu(:,:))*FDY_U(array(:,:,n)) )  
     do j=jsc,jec
        do i=isc,iec
           friction(i,j,n) = tmp(i,j)*Grd%umask(i,j,1)
        enddo
     enddo
  enddo

  call mpp_update_domains (friction(:,:,:), Dom%domain2d)


end subroutine lap_friction_barotropic
! </SUBROUTINE> NAME="lap_friction_barotropic"


!#######################################################################
! <SUBROUTINE NAME="ocean_lap_friction_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_lap_friction_restart(time_stamp)
  character(len=*),           intent(in), optional :: time_stamp

  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_lap_friction (ocean_lap_friction_end): needs initialization')
  endif 

  call save_restart(Lap_restart, time_stamp)

end subroutine ocean_lap_friction_restart
! </SUBROUTINE> NAME="ocean_lap_friction_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_lap_friction_end">
!
! <DESCRIPTION>
! Write to restart of the vertically averaged viscosity. 
! </DESCRIPTION>
!
subroutine ocean_lap_friction_end(Time)

  type(ocean_time_type), intent(in) :: Time
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_lap_friction (ocean_lap_friction_end): needs initialization')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_lap_friction (ocean_lap_friction_end): NO restart written.'
    call mpp_error(WARNING,&
    '==>Warning from ocean_lap_friction_mod (ocean_lap_friction_end): NO restart written.')
    return
  endif 

  call ocean_lap_friction_restart

  write(stdoutunit,*)' ' 
  write(stdoutunit,*) 'From ocean_lap_friction_end: ending chksum'
  call write_timestamp(Time%model_time)

   write (stdoutunit, *) &
  'checksum ending lap_viscosity', mpp_chksum(lap_viscosity(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,1))


end subroutine ocean_lap_friction_end
! </SUBROUTINE> NAME="ocean_lap_friction_end"


end module ocean_lap_friction_mod
      
