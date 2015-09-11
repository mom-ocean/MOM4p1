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
module ocean_bih_friction_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! This module calls the appropriate lateral biharmonic friction modules. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module serves as an interface to the chosen lateral 
! biharmonic friction modules.  
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
!<NAMELIST NAME="ocean_bih_friction_nml">
! <DATA NAME="bih_friction_scheme" TYPE="character">
! To determine the biharmonic friction scheme: "const" or "general"
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

use constants_mod,    only: epsln
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: open_namelist_file, check_nml_error, write_version_number, close_file
use fms_mod,          only: FATAL, WARNING, stdout, stdlog
use fms_mod,          only: file_exist
use fms_io_mod,       only: register_restart_field, save_restart, restore_state
use fms_io_mod,       only: restart_file_type
use mpp_domains_mod,  only: mpp_update_domains
use mpp_mod,          only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_MODULE
use mpp_mod,          only: mpp_error, mpp_chksum

use ocean_bihcst_friction_mod, only: ocean_bihcst_friction_init, bihcst_friction
use ocean_bihcst_friction_mod, only: bihcst_viscosity_check, bihcst_reynolds_check
use ocean_bihgen_friction_mod, only: ocean_bihgen_friction_init, bihgen_friction
use ocean_bihgen_friction_mod, only: bihgen_viscosity_check, bihgen_reynolds_check
use ocean_domains_mod,         only: get_local_indices
use ocean_operators_mod,       only: BAX, BAY, BDX_EU, BDY_NU, FDX_U, FDY_U, FMX, FMY
use ocean_parameters_mod,      only: missing_value
use ocean_types_mod,           only: ocean_time_type, ocean_grid_type
use ocean_types_mod,           only: ocean_domain_type
use ocean_types_mod,           only: ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,           only: ocean_options_type 
use ocean_util_mod,            only: write_timestamp

implicit none

private

public ocean_bih_friction_init
public ocean_bih_friction_end
public bih_friction 
public bih_viscosity_check
public bih_reynolds_check
public bih_friction_barotropic
public ocean_bih_friction_restart

! for choosing the mixing scheme
character(len=10) :: bih_friction_scheme='general'  ! 'const' or 'general'
integer           :: MIX_SCHEME         = 1         ! 1=const, 2=general
logical           :: debug_this_module  = .false.

! for deciding whether we are using biharmonic friction 
logical :: use_bihcst_friction=.false. 
logical :: use_bihgen_friction=.false. 
logical :: use_this_module    =.false. 

! vertically averaged viscosity on east and north face of U cells (m^4/s)
real, dimension(:,:),   allocatable :: bih_viscosity
real, dimension(:,:),   allocatable :: bih_visc_ceu
real, dimension(:,:),   allocatable :: bih_visc_cnu
real, dimension(:,:,:), allocatable :: del2_vel

! for diagnostics 
integer :: id_bih_viscosity=-1 
logical :: used

! for clocks
integer :: id_clock_bihcst
integer :: id_clock_bihgen

! for restart
type(restart_file_type), save :: Bih_restart

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=256) :: version=&
     '=>Using: ocean_bih_friction.f90 ($Id: ocean_bih_friction.F90,v 16.0.42.1.54.1 2009/10/10 00:42:13 nnz Exp $)'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: module_is_initialized = .false.
logical :: write_a_restart       = .true.

namelist /ocean_bih_friction_nml/ bih_friction_scheme, debug_this_module, write_a_restart

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_bih_friction_init">
!
! <DESCRIPTION>
! Initialize the horizontal biharmonic friction module.
! </DESCRIPTION>
!
subroutine ocean_bih_friction_init(Grid, Domain, Time, Ocean_options, dtime, obc, debug)

  type(ocean_grid_type),   target, intent(in)    :: Grid
  type(ocean_domain_type), target, intent(in)    :: Domain
  type(ocean_time_type),           intent(in)    :: Time
  type(ocean_options_type),        intent(inout) :: Ocean_options
  real,                            intent(in)    :: dtime
  logical,                         intent(in)    :: obc
  logical,               optional, intent(in)    :: debug

  integer :: ioun, io_status, ierr
  integer :: num_methods=0 
  integer :: id_restart

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bih_friction_mod (ocean_bih_friction_init): already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_bih_friction_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_bih_friction_nml)  
  write (stdlogunit,ocean_bih_friction_nml)
  ierr = check_nml_error(io_status,'ocean_bih_friction_nml')
  call close_file(ioun)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  Grd => Grid
  Dom => Domain

  if (PRESENT(debug) .and. .not. debug_this_module) then 
    debug_this_module = debug
  endif 

  ! for clocks 
  id_clock_bihcst = mpp_clock_id('(Ocean const   bih frict) ' ,grain=CLOCK_MODULE)
  id_clock_bihgen = mpp_clock_id('(Ocean general bih frict) ' ,grain=CLOCK_MODULE)

  ! vertically averaged vicsosities 
  allocate (bih_viscosity(isd:ied,jsd:jed))
  allocate (bih_visc_ceu(isd:ied,jsd:jed))
  allocate (bih_visc_cnu(isd:ied,jsd:jed))
  allocate (del2_vel(isd:ied,jsd:jed,2))
  bih_viscosity(:,:) = 0.0
  bih_visc_ceu(:,:)  = 0.0
  bih_visc_cnu(:,:)  = 0.0
  del2_vel(:,:,:)    = 0.0

  ! initialize biharmonic friction scheme
  if(bih_friction_scheme == 'const') then 
      num_methods=num_methods+1
      write(stdoutunit,'(a)') &
      '==>Note from ocean_bih_friction_init: constant biharmonic friction scheme used.'
      MIX_SCHEME=1
      call ocean_bihcst_friction_init (Grid, Domain, Time, Ocean_options, dtime, &
                                       obc, use_bihcst_friction, debug_this_module)
  elseif(bih_friction_scheme == 'general') then 
      num_methods=num_methods+1
      write(stdoutunit,'(a)') &
      '==>Note from ocean_bih_friction_init: general biharmonic friction scheme used.'
      MIX_SCHEME=2
      call ocean_bihgen_friction_init (Grid, Domain, Time, Ocean_options, dtime, &
                                       obc, use_bihgen_friction, debug_this_module)
  endif

  if(num_methods==0) then 
      write(stdoutunit,'(a)') &
      '==>Error: ocean_bih_friction_mod: pick a biharmonic friction scheme: "const" or "general".'
      call mpp_error(FATAL,&
      '==>Error: ocean_bih_friction_mod: pick a biharmonic friction scheme: "const" or "general".')
  elseif(num_methods>1) then 
      write(stdoutunit,'(a)')&
      '==>Error: ocean_bih_friction_mod: choose only one biharmonic friction scheme.'
      call mpp_error(FATAL,&
      '==>Error: ocean_bih_friction_mod: choose only one biharmonic friction scheme.')
  endif

  ! to determine whether we are doing any biharmonic friction 
  if(use_bihcst_friction .or. use_bihgen_friction) then 
     use_this_module = .true. 
  endif 

  if(use_this_module) then 
      id_restart = register_restart_field(Bih_restart, 'ocean_bih_friction.res.nc','bih_viscosity', &
                   bih_viscosity,domain=Dom%domain2d)
      if(.NOT.file_exist('INPUT/ocean_bih_friction.res.nc')) then
          if (.NOT. Time%init) then 
              call mpp_error(FATAL,'Expecting file INPUT/ocean_bih_friction.res.nc to exist.&
                   &This file was not found and Time%init=.false.')
          endif
      else
          call restore_state(Bih_restart)
          call mpp_update_domains(bih_viscosity,Dom%domain2d)
          write (stdoutunit,'(1x,a)') 'bih_viscosity being read from restart.'
          write (stdoutunit, *) 'checksum start bih_viscosity = ', &
               mpp_chksum(bih_viscosity(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,1))

          bih_visc_ceu(:,:) = BAY(bih_viscosity(:,:))
          bih_visc_cnu(:,:) = BAX(bih_viscosity(:,:))
          call mpp_update_domains (bih_visc_ceu, Dom%domain2d)
          call mpp_update_domains (bih_visc_cnu, Dom%domain2d)
          bih_visc_ceu(:,:) = sqrt(bih_visc_ceu(:,:))
          bih_visc_cnu(:,:) = sqrt(bih_visc_cnu(:,:))

      endif

      if(.not. write_a_restart) then 
          write(stdoutunit,'(a)') '==>Note: running ocean_bih_friction with write_a_restart=.false.'
          write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
      endif

  endif

  id_bih_viscosity = register_diag_field ('ocean_model', 'bih_viscosity', Grd%vel_axes_uv(1:2), &
                     Time%model_time,'Vertically averaged biharmonic viscosity',                &
                    'm^4/s',missing_value=missing_value, range=(/-10.0,1.e20/))


end subroutine ocean_bih_friction_init
! </SUBROUTINE>  NAME="ocean_bih_friction_init"


!#######################################################################
! <SUBROUTINE NAME="bih_friction">
!
! <DESCRIPTION>
! Compute the thickness weighted and density weighted accel due to 
! lateral biharmonic friction.  Add this contribution to Velocity%accel. 
! </DESCRIPTION>
!
subroutine bih_friction(Time, Thickness, Velocity, energy_analysis_step)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_thickness_type), intent(in)    :: Thickness
  type(ocean_velocity_type),  intent(inout) :: Velocity
  logical,                    intent(in)    :: energy_analysis_step

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bih_friction_mod (bih_friction): module needs to be initialized')
  endif 

  if(MIX_SCHEME==1) then 
    call mpp_clock_begin(id_clock_bihcst)
    call bihcst_friction(Time, Thickness, Velocity, bih_viscosity, energy_analysis_step)
    call mpp_clock_end(id_clock_bihcst)
  elseif(MIX_SCHEME==2) then 
    call mpp_clock_begin(id_clock_bihgen)
    call bihgen_friction(Time, Thickness, Velocity, bih_viscosity, energy_analysis_step)
    call mpp_clock_end(id_clock_bihgen)
  endif 

  ! vertically averaged viscosity for north and east face of U-cells.
  ! apply the sqrt so to apply bih_viscosity on each iteration 
  ! of the laplacian.
  if(use_this_module) then 
     bih_visc_ceu(:,:) = BAY(bih_viscosity(:,:))
     bih_visc_cnu(:,:) = BAX(bih_viscosity(:,:))
     call mpp_update_domains (bih_visc_ceu, Dom%domain2d)
     call mpp_update_domains (bih_visc_cnu, Dom%domain2d)
     bih_visc_ceu(:,:) = sqrt(bih_visc_ceu(:,:))
     bih_visc_cnu(:,:) = sqrt(bih_visc_cnu(:,:))
  endif 

   if (id_bih_viscosity > 0) then 
      used = send_data(id_bih_viscosity, bih_viscosity(:,:), &
             Time%model_time, rmask=Grd%umask(:,:,1),        &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
   endif 


end subroutine bih_friction
! </SUBROUTINE> NAME="bih_friction"


!#######################################################################
! <SUBROUTINE NAME="bih_viscosity_check">
!
! <DESCRIPTION>
! To check that the viscosity is not too large. 
! </DESCRIPTION>
!
subroutine bih_viscosity_check

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_bih_friction_mod (bih_viscosity_check): needs initialization')
  endif 

  if(MIX_SCHEME==1) then 
    call bihcst_viscosity_check
  elseif(MIX_SCHEME==2) then 
    call bihgen_viscosity_check
  endif 

end subroutine bih_viscosity_check
! </SUBROUTINE> NAME="bih_viscosity_check"


!#######################################################################
! <SUBROUTINE NAME="bih_reynolds_check">
!
! <DESCRIPTION>
! To check that the Reynolds number is not too large. 
! </DESCRIPTION>
!
subroutine bih_reynolds_check(Time, Velocity)

  type(ocean_time_type),     intent(in) :: Time
  type(ocean_velocity_type), intent(in) :: Velocity

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_bih_friction_mod (bih_reynolds_check): needs initialization')
  endif 

  if(MIX_SCHEME==1) then 
    call bihcst_reynolds_check(Time, Velocity)
  elseif(MIX_SCHEME==2) then 
    call bihgen_reynolds_check(Time, Velocity)
  endif 

end subroutine bih_reynolds_check
! </SUBROUTINE> NAME="bih_reynolds_check"


!#######################################################################
! <SUBROUTINE NAME="bih_friction_barotropic">
!
! <DESCRIPTION>
!
! This routine computes the biharmonic friction acting on a two-dim
! array.  It uses the two-dimensional vertically averaged viscosity 
! used in the biharmonic friction module.  The intent is to apply this
! 2d operator to the vertically integrated horizontal momentum. We 
! ignore the spherical metric terms in this form of the operator, 
! since we are aiming for a fast smoothing operator to be applied
! during each of the many barotropic time steps.  We also apply 
! just the isotropic portion of the more general anisotropic 
! biharmonic operator.  
!
! </DESCRIPTION>
!
subroutine bih_friction_barotropic(bih_ceu_back, bih_cnu_back, array, friction)

  real, dimension(isd:,jsd:),   intent(in)    :: bih_ceu_back
  real, dimension(isd:,jsd:),   intent(in)    :: bih_cnu_back
  real, dimension(isd:,jsd:,:), intent(in)    :: array
  real, dimension(isd:,jsd:,:), intent(inout) :: friction

  integer :: i, j, n
  real    :: tmp(isd:ied,jsd:jed,2)

  ! first iteration 
  tmp(:,:,:) = 0.0
  do n=1,2
     tmp(:,:,n) =  BDX_EU( (bih_ceu_back(:,:)+bih_visc_ceu(:,:))*FDX_U(array(:,:,n)) ) &
                  +BDY_NU( (bih_cnu_back(:,:)+bih_visc_cnu(:,:))*FDY_U(array(:,:,n)) )  
     do j=jsc,jec
        do i=isc,iec
           del2_vel(i,j,n) = tmp(i,j,n)*Grd%umask(i,j,1)
        enddo
     enddo
  enddo
  call mpp_update_domains (del2_vel(:,:,:), Dom%domain2d)

  ! second iteration 
  tmp(:,:,:) = 0.0
  do n=1,2
     tmp(:,:,n) =  BDX_EU( (bih_ceu_back(:,:)+bih_visc_ceu(:,:))*FDX_U(del2_vel(:,:,n)) ) &
                  +BDY_NU( (bih_cnu_back(:,:)+bih_visc_cnu(:,:))*FDY_U(del2_vel(:,:,n)) )  
     do j=jsc,jec
        do i=isc,iec
           friction(i,j,n) = -tmp(i,j,n)*Grd%umask(i,j,1)
        enddo
     enddo
  enddo

  call mpp_update_domains (friction(:,:,:), Dom%domain2d)


end subroutine bih_friction_barotropic
! </SUBROUTINE> NAME="bi_friction_barotropic"


!#######################################################################
! <SUBROUTINE NAME="ocean_bih_friction_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_bih_friction_restart(time_stamp)
  character(len=*),           intent(in), optional :: time_stamp

  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_bih_friction (ocean_bih_friction_end): needs initialization')
  endif 

  call save_restart(Bih_restart, time_stamp)

end subroutine ocean_bih_friction_restart
! </SUBROUTINE> NAME="ocean_bih_friction_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_bih_friction_end">
!
! <DESCRIPTION>
! Write to restart of the vertically averaged viscosity. 
! </DESCRIPTION>
!
subroutine ocean_bih_friction_end(Time)

  type(ocean_time_type), intent(in) :: Time
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. use_this_module) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_bih_friction (ocean_bih_friction_end): needs initialization')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_bih_friction (ocean_bih_friction_end): NO restart written.'
    call mpp_error(WARNING,&
    '==>Warning from ocean_bih_friction_mod (ocean_bih_friction_end): NO restart written.')
    return
  endif 

  call ocean_bih_friction_restart

  write(stdoutunit,*)' ' 
  write(stdoutunit,*) 'From ocean_bih_friction_end: ending chksum'
  call write_timestamp(Time%model_time)

   write (stdoutunit, *) &
  'checksum ending bih_viscosity', mpp_chksum(bih_viscosity(isc:iec,jsc:jec)*Grd%umask(isc:iec,jsc:jec,1))


end subroutine ocean_bih_friction_end
! </SUBROUTINE> NAME="ocean_bih_friction_end"


end module ocean_bih_friction_mod
      
      
