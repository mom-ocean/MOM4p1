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
module ocean_rivermix_mod
!  
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Keith.Dixon@noaa.gov">  K.W. Dixon 
!</CONTACT>
!
!<OVERVIEW>
! Tracer source from discharging river with depth or 
! mixing rivers with depth. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness weighted tendency [tracer*rho*meter/sec]
! associated with discharge of river tracer content 
! over a user defined column of ocean points. Points are
! selected based on whether river flow into a point is nonzero.
! Contribution added to tracer source array.
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R. C. Pacanowski, and A. Rosati
! A Guide to MOM4 (2003)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <NOTE>
! Algorithm ensures total tracer is conserved.  Note that volume/mass is 
! modified by river water within the eta-equation using the big leap-frog.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_rivermix_nml">
!  <DATA NAME="use_this_module"  TYPE="logical">
!  Must be true to enable this module.  Default=.true., since
!  this is the only way that tracer in river water enters the ocean. 
!  </DATA>
!
! <DATA NAME="discharge_combine_runoff_calve"  TYPE="logical">
!  Set discharge_combine_runoff_calve=.true. to discharge combined tracer carried
!  by liquid and solid runoff. This approach is sensible when ocean 
!  assigns a tracer content to the liquid and solid runoff fields.  
!  The alternative is to have a land model that provides information about
!  the tracer coming into the ocean from land water, in which case it is 
!  preferable to set discharge_combine_runoff_calve=.false., so to do the runoff 
!  and calving separately.  
!  Default discharge_combine_runoff_calve=.true. 
!  </DATA>
!
!  <DATA NAME="river_insertion_thickness" TYPE="real" UNITS="meter">
!  Thickness of the column over which to insert tracers from 
!  rivers. Default river_insertion_thickness=0.0 (all in top).
!  </DATA>
!  <DATA NAME="runoff_insertion_thickness" TYPE="real" UNITS="meter">
!  Thickness of the column over which to insert tracers carried by 
!  liquid runoff. Default runoff_insertion_thickness=0.0 (all in top).
!  </DATA>
!  <DATA NAME="calving_insertion_thickness" TYPE="real" UNITS="meter">
!  Thickness of the column over which to insert tracers carried by 
!  solid runoff. Default calving_insertion_thickness=0.0 (all in top).
!  </DATA>
!
!  <DATA NAME="river_diffusion_thickness" TYPE="real" UNITS="meter">
!  Thickness of the column over which to diffuse tracers from 
!  rivers. 
!  </DATA> 
!  <DATA NAME="river_diffusivity" TYPE="real" UNITS="m^2/s">
!  Vertical diffusivity enhancement at river mouths which is applied 
!  to a depth of river_diffusion_thickness, with linear tapering to zero
!  enhancement from the ocean surface to river_diffusion_thickness. 
!  </DATA> 
!  <DATA NAME="river_diffuse_temp" TYPE="logical">
!  Logical to determine if enhance vertical diffusion of temp at river mouths 
!  </DATA> 
!  <DATA NAME="river_diffuse_salt" TYPE="logical">
!  Logical to determine if enhance vertical diffusion of salt and all other 
!  passive tracers at river mouths 
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging 
!  </DATA> 
!  <DATA NAME="debug_all_in_top_cell" TYPE="logical">
!  For debugging, by placing all in top cell, regardless value of 
!  river_insertion_thickness.
!  </DATA> 
!  <DATA NAME="debug_this_module_heat" TYPE="logical">
!  For debugging, print global sum of heating rate by river water. 
!  </DATA> 
!</NAMELIST>
!

use axis_utils_mod,      only: frac_index, nearest_index
use constants_mod,       only: epsln
use diag_manager_mod,    only: register_diag_field, send_data
use fms_mod,             only: write_version_number, open_namelist_file, check_nml_error, close_file
use fms_mod,             only: FATAL, NOTE, stdout, stdlog 
use mpp_domains_mod,     only: domain2D, mpp_global_sum, BITWISE_EXACT_SUM
use mpp_mod,             only: mpp_error

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_thickness_type 
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_external_mode_type, ocean_density_type
use ocean_workspace_mod,  only: wrk1, wrk1_2d 

implicit none

private 

public ocean_rivermix_init
public rivermix
private river_discharge_tracer
private runoff_calving_discharge_tracer
private river_kappa
private river_diagnostics


type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

real    :: river_insertion_thickness=0.0   ! static thickness (m) of ocean column where discharge rivers.  
                                           ! actual thickness is based on model grid spacing. min thickness=dtz(1).
real    :: runoff_insertion_thickness=0.0  ! static thickness (m) of ocean column where discharge liquid runoff.
                                           ! actual thickness is based on model grid spacing. min thickness=dtz(1).
real    :: calving_insertion_thickness=0.0 ! static thickness (m) of ocean column where discharge solid runoff.
                                           ! actual thickness is based on model grid spacing. min thickness=dtz(1).

real    :: river_diffusion_thickness=0.0 ! static thickness (m) of ocean column where diffuse tracer at river mouths.    
                                         ! actual thickness is based on model grid spacing. min thickness=dtz(1).      
real    :: river_diffusivity=5.0e-3      ! enhancement to the vertical diffusivity (m^2/s) at river mouths 
logical :: river_diffuse_temp=.false.    ! to enhance diffusivity of temp at river mouths over river_thickness column
logical :: river_diffuse_salt=.false.    ! to enhance diffusivity of salt at river mouths over river_thickness column

logical :: discharge_combine_runoff_calve=.true. ! to discharge liquid+solid together  


real, dimension(:,:), allocatable :: tracer_flux


! for diagnostics 
logical :: used
integer, dimension(:), allocatable :: id_rivermix
integer, dimension(:), allocatable :: id_runoffmix
integer, dimension(:), allocatable :: id_calvingmix
integer :: id_diff_cbt_river_t       =-1
integer :: id_diff_cbt_river_s       =-1
integer :: id_neutral_rho_rivermix   =-1
integer :: id_neutral_rho_runoffmix  =-1
integer :: id_neutral_rho_calvingmix =-1
integer :: id_wdian_rho_rivermix     =-1
integer :: id_wdian_rho_runoffmix    =-1
integer :: id_wdian_rho_calvingmix   =-1

integer :: unit=6       ! processor zero writes to unit 6

! time step (sec)
real    :: dtime 

integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1

character(len=128) :: version=&
       '$Id: ocean_rivermix.F90,v 16.0.2.1.10.1.2.1.52.1 2009/10/10 00:42:46 nnz Exp $'
character (len=128) :: tagname=&
     '$Name: mom4p1_pubrel_dec2009_nnz $'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

logical :: debug_this_module     = .false. 
logical :: debug_all_in_top_cell = .false.
logical :: debug_this_module_heat=.false.
logical :: module_is_initialized = .FALSE.
logical :: use_this_module       = .true.

namelist /ocean_rivermix_nml/ use_this_module, debug_this_module,   &
         debug_all_in_top_cell, debug_this_module_heat,             &
         river_diffuse_temp, river_diffuse_salt,                    & 
         river_diffusion_thickness, river_diffusivity,              &
         discharge_combine_runoff_calve, river_insertion_thickness, &
         runoff_insertion_thickness, calving_insertion_thickness

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_rivermix_init">
!
! <DESCRIPTION>
! Initial set up for mixing of tracers at runoff points. 
! </DESCRIPTION>
!
  subroutine ocean_rivermix_init(Grid, Domain, Time, Time_steps, T_prog, Ocean_options, debug)
  
    type(ocean_grid_type),        target, intent(in) :: Grid
    type(ocean_domain_type),      target, intent(in) :: Domain
    type(ocean_time_type),        intent(in)         :: Time 
    type(ocean_time_steps_type),  intent(in)         :: Time_steps
    type(ocean_prog_tracer_type), intent(inout)      :: T_prog(:)
    type(ocean_options_type),     intent(inout)      :: Ocean_options
    logical, optional,            intent(in)         :: debug

    integer :: n, nz_insert, nz_diffuse
    integer :: io_status, ioun, ierr

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
      call mpp_error(FATAL, &
       '==>Error from ocean_rivermix_mod (ocean_rivermix_init): module already initialized')
    endif 

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
    nk = Grid%nk

    Dom => Domain
    Grd => Grid

    dtime = Time_steps%dtts

    num_prog_tracers = size(T_prog(:))
    do n=1,num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo
    
    if (PRESENT(debug)) debug_this_module = debug

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_rivermix_nml,IOSTAT=io_status)
    write (stdlogunit,ocean_rivermix_nml)
    write (stdoutunit,'(/)') 
    write (stdoutunit,ocean_rivermix_nml)
    ierr = check_nml_error(io_status,'ocean_rivermix_nml')
    call close_file (ioun)

    if(use_this_module) then 
      call mpp_error(NOTE,&
      '==>From ocean_rivermix_mod: Using rivermix module to mix liquid and/or solid runoff into the ocean.')
      Ocean_options%rivermix = 'Used rivermix module to mix liquid and/or solid runoff into ocean.'
    else
      call mpp_error(NOTE,&
      '==>From ocean_rivermix_mod: NOT using rivermix module. Will NOT mix river tracers into ocean. Is this what you wish?')
      Ocean_options%rivermix = 'Did NOT use rivermix module. River tracers NOT mixed into the ocean. '
      return
    endif 


    ! allocate for tracer flux in runoff or calving 
    allocate (tracer_flux(isd:ied,jsd:jed))
    tracer_flux=0.0

 
    if(discharge_combine_runoff_calve) then 
        write(stdoutunit,'(a)')'==>Note: discharging calving+runoff together. The alternative is to separately discharge. '

        ! get nominal thickness 
        nz_insert  = nearest_index(river_insertion_thickness,Grd%zw)
        nz_insert  = max(1,nz_insert)
        if(river_insertion_thickness==0.0) then 
            call mpp_error(NOTE, &
            '==>Note: resetting river_insertion_thickness to dzt(1). Discharge all river into top cell')  
        else 
            write(stdoutunit,'(a)')'==>Note: if using waterflux and rivers, then will'
            write(stdoutunit,'(a,i3,a)')'   discharge river tracer over ',nz_insert,' grid points in vertical'
        endif

    else 

        ! get nominal thickness 
        nz_insert  = nearest_index(runoff_insertion_thickness,Grd%zw)
        nz_insert  = max(1,nz_insert)
        if(runoff_insertion_thickness==0.0) then 
            call mpp_error(NOTE, &
            '==>Note: resetting runoff_insertion_thickness to dzt(1). Discharge all liquid runoff into top cell')  
        else 
            write(stdoutunit,'(a)')'==>Note: if using waterflux and liquid runoff, then will'
            write(stdoutunit,'(a,i3,a)')'   discharge liquid runoff tracer over ',nz_insert,' grid points in vertical'
        endif


        ! get nominal thickness 
        nz_insert  = nearest_index(calving_insertion_thickness,Grd%zw)
        nz_insert  = max(1,nz_insert)
        if(calving_insertion_thickness==0.0) then 
            call mpp_error(NOTE, &
            '==>Note: resetting calving_insertion_thickness to dzt(1). Discharge all solid runoff into top cell')  
        else 
            write(stdoutunit,'(a)')'==>Note: if using waterflux and solid runoff, then will'
            write(stdoutunit,'(a,i3,a)')'   discharge solid runoff tracer over ',nz_insert,' grid points in vertical'
        endif

    endif


    ! get nominal thickness over which to diffuse river
    nz_diffuse = nearest_index(river_diffusion_thickness,Grd%zw)    
    nz_diffuse = max(1,nz_diffuse)    

    if(river_diffuse_temp) then 
        if(Time_steps%aidif==0.0) then
            call mpp_error(FATAL, &
             '==>Error: aidif=1.0 must be set to allow stable time stepping with river_diffuse_temp=.true.')  
        endif
        if(river_diffusion_thickness==0.0) then 
            call mpp_error(NOTE, &
            '==>Note: resetting river_diffusion_thickness=dzt(1) for enhanced diff_cbt of temperature.')  
        else  
            write(stdoutunit,'(a,i3,a)') &
            ' ==>Note: enhance temp diff_cbt at rivers over ',nz_diffuse,' points in vertical'
            write(stdoutunit,'(a,e12.5,a)')&
            '          Do so by adding diffusivity of ',river_diffusivity,' m^2/s to diff_cbt(temp)'
        endif
    endif

    if(river_diffuse_salt) then 
        if(Time_steps%aidif==0.0) then
            call mpp_error(FATAL, &
             ' ==>Error: aidif=1.0 must be set to allow stable time stepping with river_diffuse_salt=.true.')  
        endif
        if(river_diffusion_thickness==0.0) then 
            call mpp_error(NOTE, &
            '==>Note: resetting river_diffusion_thickness=dzt(1) for enhanced diff_cbt of salt and passive')  
        else 
            write(stdoutunit,'(a,i3,a)') &
            ' ==>Note: enhance salt/passive diff_cbt at rivers over ',nz_diffuse,' points in vertical'
            write(stdoutunit,'(a,e12.5,a)') &
            '          Do so by adding diffusivity of ',river_diffusivity,' m^2/s to diff_cbt(salt)'
        endif
    endif

    ! register for diag_manager 
    allocate (id_rivermix(num_prog_tracers))
    allocate (id_runoffmix(num_prog_tracers))
    allocate (id_calvingmix(num_prog_tracers))
    id_rivermix   = -1
    id_runoffmix  = -1
    id_calvingmix = -1

    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp') then 
           id_rivermix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_rivermix', &
                            Grd%tracer_axes(1:3), Time%model_time,                                 &
                            'cp*rivermix*rho_dzt*temp', 'Watt/m^2',                                & 
                            missing_value=missing_value, range=(/-1.e10,1.e10/))
           id_runoffmix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_runoffmix',&
                            Grd%tracer_axes(1:3), Time%model_time,                                  &
                            'cp*runoffmix*rho_dzt*temp', 'Watt/m^2',                                & 
                            missing_value=missing_value, range=(/-1.e10,1.e10/))
           id_calvingmix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_calvingmix',&
                            Grd%tracer_axes(1:3), Time%model_time,                                    &
                            'cp*calvingmix*rho_dzt*temp', 'Watt/m^2',                                 & 
                            missing_value=missing_value, range=(/-1.e10,1.e10/))
       else
           id_rivermix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_rivermix', &
                            Grd%tracer_axes(1:3), Time%model_time,                                 &
                            'rivermix*rho_dzt*tracer for '//trim(T_prog(n)%name),                  &
                            trim(T_prog(n)%flux_units),                                            &
                            missing_value=missing_value, range=(/-1.e10,1.e10/))
           id_runoffmix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_runoffmix',&
                            Grd%tracer_axes(1:3), Time%model_time,                                  &
                            'runoffmix*rho_dzt*tracer for '//trim(T_prog(n)%name),                  &
                            trim(T_prog(n)%flux_units),                                             &
                            missing_value=missing_value, range=(/-1.e10,1.e10/))
           id_calvingmix(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_calvingmix',&
                            Grd%tracer_axes(1:3), Time%model_time,                                    &
                            'calvingmix*rho_dzt*tracer for '//trim(T_prog(n)%name),                   &
                            trim(T_prog(n)%flux_units),                                               &
                            missing_value=missing_value, range=(/-1.e10,1.e10/))
       endif
    enddo

    if(debug_all_in_top_cell) then 
        call mpp_error(NOTE, &
        '==>Note: in rivermix module, debug_all_in_top_cell=.true., so placing all river water to top cell.')
    endif

   id_diff_cbt_river_t = -1
   id_diff_cbt_river_t = register_diag_field ('ocean_model', 'diff_cbt_river_t', &
                         Grid%tracer_axes(1:3), Time%model_time,                 &
                         'diff_cbt(temp) enhancement at rivers', 'm^2/s',        &
                         missing_value=missing_value, range=(/-10.0,1e6/))

   id_diff_cbt_river_s = -1
   id_diff_cbt_river_s = register_diag_field ('ocean_model', 'diff_cbt_river_s', &
                         Grid%tracer_axes(1:3), Time%model_time,                 &
                         'diff_cbt(salt) enhancement at rivers', 'm^2/s',        &
                         missing_value=missing_value,range=(/-10.0,1e6/))

  id_neutral_rho_rivermix = register_diag_field ('ocean_model', 'neutral_rho_rivermix', &
                            Grd%tracer_axes(1:3), Time%model_time,                      &
                            'update of neutral density from rivermix scheme',           &
                            'rho*rho_dz/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_neutral_rho_runoffmix = register_diag_field ('ocean_model', 'neutral_rho_runoffmix', &
                            Grd%tracer_axes(1:3), Time%model_time,                        &
                            'update of neutral density from runoffmix scheme',            &
                            'rho*rho_dz/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_neutral_rho_calvingmix = register_diag_field ('ocean_model', 'neutral_rho_calvingmix', &
                            Grd%tracer_axes(1:3), Time%model_time,                          &
                            'update of neutral density from calvingmix scheme',             &
                            'rho*rho_dz/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_wdian_rho_rivermix = register_diag_field ('ocean_model', 'wdian_rho_rivermix', &
                          Grd%tracer_axes(1:3), Time%model_time,                    &
                          'dianeutral velocity component due to rivermix scheme',   &
                          'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_wdian_rho_runoffmix = register_diag_field ('ocean_model', 'wdian_rho_runoffmix', &
                          Grd%tracer_axes(1:3), Time%model_time,                      &
                          'dianeutral velocity component due to runoffmix scheme',    &
                          'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_wdian_rho_calvingmix = register_diag_field ('ocean_model', 'wdian_rho_calvingmix',&
                          Grd%tracer_axes(1:3), Time%model_time,                       &
                          'dianeutral velocity component due to calvingmix scheme',    &
                          'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))


  end subroutine ocean_rivermix_init
! </SUBROUTINE> NAME="ocean_rivermix_init"


!#######################################################################
! <SUBROUTINE NAME="rivermix">
!
! <DESCRIPTION>
! This subroutine computes one or all of the following: 
!
! (1) Thickness weighted and rho weighted tracer source associated 
! with river tracer content discharged into a vertical column of ocean 
! tracer cells. This is done if river_discharge=.true.
!
! (2) Enhance vertical diffusivity at river mouths. 
! This is done if river_diffuse_temp=.true. or 
! river_diffuse_salt=.true. 
!
! Doing one or both are useful for models with fine vertical  
! resolution, where discharging river content to top cell 
! is often not numerically suitable nor physically relevant.
!
! </DESCRIPTION>
!
subroutine rivermix (Time, Thickness, Dens, T_prog, river, runoff, calving, &
                     diff_cbt, index_temp, index_salt)

  type(ocean_time_type),          intent(in)     :: Time
  type(ocean_thickness_type),     intent(in)     :: Thickness
  type(ocean_density_type),       intent(in)     :: Dens
  type(ocean_prog_tracer_type),   intent(inout)  :: T_prog(:)
  integer,                        intent(in)     :: index_temp
  integer,                        intent(in)     :: index_salt
  real, dimension(isd:,jsd:),     intent(in)     :: river
  real, dimension(isd:,jsd:),     intent(in)     :: runoff
  real, dimension(isd:,jsd:),     intent(in)     :: calving 
  real, dimension(isd:,jsd:,:,:), intent(inout)  :: diff_cbt

  integer :: i,j,n
  logical :: river_discharge  =.false.
  logical :: runoff_discharge =.false.
  logical :: calving_discharge=.false.

  if(.not. use_this_module) return

  if(.not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_rivermix_mod (rivermix): module must be initialized')
  endif 


  ! for discharging the combined runoff+calving fields  
  if(discharge_combine_runoff_calve) then 

      do n=1,num_prog_tracers  
        T_prog(n)%wrk1(:,:,:) = 0.0
      enddo  

      river_discharge=.false. 
      do j=jsc,jec
         do i=isc,iec
            if(river(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) river_discharge=.true.
         enddo
      enddo
      if(river_discharge) then 
          call river_discharge_tracer(Time, Thickness, Dens, T_prog(1:num_prog_tracers), river)
      endif

      do n=1,num_prog_tracers 
         if(id_rivermix(n) > 0) then 
             used = send_data (id_rivermix(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion, &
                  Time%model_time, rmask=Grd%tmask(:,:,:),                                 &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
         endif
      enddo
      call river_diagnostics(Time, Thickness, Dens, &
      T_prog(index_temp)%wrk1(:,:,:),T_prog(index_salt)%wrk1(:,:,:),3)


  ! for separately discharging liquid runoff and solid calving 
  else

      do n=1,num_prog_tracers  
        T_prog(n)%wrk1(:,:,:) = 0.0
      enddo  
      runoff_discharge=.false. 
      do j=jsc,jec
         do i=isc,iec
            if(runoff(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) runoff_discharge=.true.
         enddo
      enddo
      if(runoff_discharge) then
          call runoff_calving_discharge_tracer(Time, Thickness, Dens, &
          T_prog(1:num_prog_tracers), runoff, runoff_insertion_thickness, 1)
      endif 
      do n=1,num_prog_tracers 
         if(id_runoffmix(n) > 0) then 
             used = send_data (id_runoffmix(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion, &
                    Time%model_time, rmask=Grd%tmask(:,:,:),                                &
                    is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
         endif
      enddo
      call river_diagnostics(Time, Thickness, Dens, &
      T_prog(index_temp)%wrk1(:,:,:),T_prog(index_salt)%wrk1(:,:,:),1)


      do n=1,num_prog_tracers  
        T_prog(n)%wrk1(:,:,:) = 0.0
      enddo  
      calving_discharge=.false. 
      do j=jsc,jec
         do i=isc,iec
            if(calving(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) calving_discharge=.true.
         enddo
      enddo
      if(calving_discharge) then
          call runoff_calving_discharge_tracer(Time, Thickness, Dens, &
          T_prog(1:num_prog_tracers), calving, calving_insertion_thickness, 2)
      endif 
      do n=1,num_prog_tracers 
         if(id_calvingmix(n) > 0) then 
             used = send_data (id_calvingmix(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion, &
                  Time%model_time, rmask=Grd%tmask(:,:,:),                                   &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
         endif
      enddo
      call river_diagnostics(Time, Thickness, Dens, &
      T_prog(index_temp)%wrk1(:,:,:),T_prog(index_salt)%wrk1(:,:,:),2)

  endif  ! discharge_combine_runoff_calve


  if(river_diffuse_temp) then 
     call river_kappa(Time, Thickness, T_prog(index_temp), diff_cbt(isd:ied,jsd:jed,:,1))
  endif 
  if(river_diffuse_salt) then 
     call river_kappa(Time, Thickness, T_prog(index_salt), diff_cbt(isd:ied,jsd:jed,:,2))
  endif 


end subroutine rivermix
! </SUBROUTINE> NAME="rivermix"


!#######################################################################
! <SUBROUTINE NAME="river_discharge_tracer">
!
! <DESCRIPTION>
! Compute thickness weighted tracer source [tracer*m/s]
! associated with the discharge of tracer from a river over 
! a vertical column whose thickness is set by River_insertion_thickness 
! and whose horizontal location is given by the river array. 
!
! Jan 2005: converted to mass weighting for use with non-Boussinesq
! pressure-like coodinates. 
!
! This subroutine is maintained for legacy purposes.
!
! </DESCRIPTION>
!
subroutine river_discharge_tracer (Time, Thickness, Dens, T_prog, river)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(in)    :: river

  integer :: i, j, k, n, nz
  integer :: tau
  real    :: depth, thkocean
  real    :: delta(nk), delta_rho_tocean(nk), delta_rho0_triver(nk)
  real    :: zextra, zinsert, tracerextra, tracernew(nk)
  real    :: temporary, tracer_input 
  
  integer :: stdoutunit 
  stdoutunit=stdout() 


  if(.not. use_this_module) return
  
  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_rivermix_mod (river_discharge_tracer): module must be initialized')
  endif 

  tau   = Time%tau
  
  do n=1,num_prog_tracers  
    T_prog(n)%wrk1(:,:,:) = 0.0
  enddo  
  delta             = 0.0
  delta_rho_tocean  = 0.0
  delta_rho0_triver = 0.0
  wrk1              = 0.0
               
  do j=jsc,jec
     do i=isc,iec

        if (river(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then


        ! the array "river" contains the volume rate (m/s) or mass
        ! rate (kg/m2/sec) of fluid with tracer concentration triver
        ! that is to be distributed in the vertical. 

            depth = min(Grd%ht(i,j),river_insertion_thickness)        ! be sure not to discharge river content into rock 
            nz    = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw))) ! number of k-levels into which discharge rivers
            nz    = max(1,nz)                                         ! make sure have at least one cell to discharge into

            ! determine fractional thicknesses of grid cells 
            thkocean = 0.0
            do k=1,nz
               thkocean = thkocean + Thickness%rho_dzt(i,j,k,tau)
            enddo
            do k=1,nz
               delta(k) = Thickness%rho_dzt(i,j,k,tau)/(epsln+thkocean)
            enddo


            do n=1,num_prog_tracers

! mix a portion of the river mass(volume) inserted over the course of a tracer timestep
! into each model level from k=1:nz, in proportion to the fractional cell thickness.
! The mixing scheme used here follows scheme devised by Keith Dixon to handle hydraulic
! control across xland mixing points
!
! The following schematic illustrates the scheme, assuming constant layer thickness
! neglecting free surface height variations and assuming the river insertion depth
! extends to the base of the third model level.  Each of these assumptions are for 
! purposes of this schematic only.  The scheme that is implemented works in general
! for non-constant cell thicknesses.  
!               
!                  Total depth = H; cell thickness = dz; River mass flux/timestep = R
!  x--------x
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of 
!  !        !        Triver + 2*zinsert meters of water at T*(2)
!  !        !        T*(1) = (T(1).dz + 3.T*(2).zinsert + Triver.zinsert ) / (dz + 3.zinsert)               
!  x--------x               
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of 
!  !        !        Triver + zinsert meters of water at T*(3)           
!  !        !        T*(2) = (T(2).dz + 2.T*(3).zinsert + Triver.zinsert ) / (dz + 2.zinsert)               
!  !        !
!  x--------x                              
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of Triver 
!  !        !        T*(3) = (T(3).dz + Triver.zinsert ) / (dz + zinsert)
!  !        !
!  xxxxxxxxxx               

               zextra=0.0
               do k=nz,1,-1
                  tracernew(k) = 0.0

                  if (k.eq.nz) then
                      tracerextra=0.0
                  else
                      tracerextra = tracernew(k+1)
                  endif

                  zinsert = river(i,j)*dtime*delta(k)
                  tracernew(k) = (tracerextra*zextra + T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau) + &
                                  T_prog(n)%triver(i,j)*zinsert) / (zextra+Thickness%rho_dzt(i,j,k,tau)+zinsert)

                  zextra=zextra+zinsert
               enddo

               k=1
               T_prog(n)%wrk1(i,j,k) = (tracernew(k)*(Thickness%rho_dzt(i,j,k,tau)+river(i,j)*dtime) -&
                                        T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau))/dtime
               do k=2,nz
                  T_prog(n)%wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)*(tracernew(k) - T_prog(n)%field(i,j,k,tau))/dtime
               enddo

               if(debug_all_in_top_cell) then 
                   k=1
                   T_prog(n)%wrk1(i,j,:) = 0.0
                   T_prog(n)%wrk1(i,j,k) = Grd%tmask(i,j,k)*river(i,j)*T_prog(n)%triver(i,j)
               endif

               do k=1,nz
                  T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
               enddo

               if (debug_this_module) then
                   write(6,*) 'i,j,n,river((kg/m^3)*m/sec),triver= ',i,j,n,river(i,j),T_prog(n)%triver(i,j)
                   do k=1,nz
                      write(6,*) 'i,j,k,tracer(k),tracernew(k),tendency(tracer*rho*thickness/sec)= ',&
                           i,j,k,T_prog(n)%field(i,j,k,tau),tracernew(k),T_prog(n)%wrk1(i,j,k)
                   enddo
               endif
                            
               
            enddo ! n end  

        endif ! river > 0

     enddo   ! i end
  enddo      ! j end


  if (debug_this_module_heat) then
          wrk1_2d(:,:) = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*T_prog(index_temp)%conversion &
                               *river(i,j)*T_prog(index_temp)%triver(i,j)*dtime
             enddo
          enddo
          tracer_input = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), BITWISE_EXACT_SUM)
          write(stdoutunit,'(a)')
          write (stdoutunit,'(a,es24.17,a)') ' Heat input via river (runoff+calving) in rivermix  = ',&
                                           tracer_input,' Joule'
          write(stdoutunit,'(a)')
  endif 


end subroutine river_discharge_tracer
! </SUBROUTINE> NAME="river_discharge_tracer"


!#######################################################################
! <SUBROUTINE NAME="runoff_calving_discharge_tracer">
!
! <DESCRIPTION>
! Compute thickness weighted tracer source [tracer*m/s]
! associated with the discharge of tracer from runoff or calving over 
! a vertical column whose thickness is set by either runoff_insertion_thickness 
! or calving_insertion_thickness, and whose horizontal location is given 
! by the runoff or calving array. 
!
! Jan 2005: converted to mass weighting for use with non-Boussinesq
! pressure-like coodinates. 
!
! Feb 2009: now use calving_tracer_flux and runoff_tracer_flux, as the 
! land model carries information about the tracer content in the 
! liquid and solid runoff.  
!
! </DESCRIPTION>
!
subroutine runoff_calving_discharge_tracer (Time, Thickness, Dens, T_prog, &
           river, insertion_thickness, runoff_type)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:),     intent(in)    :: river
  real,                           intent(in)    :: insertion_thickness
  integer,                        intent(in)    :: runoff_type

  integer :: i, j, k, n, nz
  integer :: tau
  real    :: depth, thkocean
  real    :: delta(nk), delta_rho_tocean(nk), delta_rho0_triver(nk)
  real    :: zextra, zinsert, tracerextra, tracernew(nk)
  real    :: temporary, tracer_input 
  
  integer :: stdoutunit 
  stdoutunit=stdout() 


  if(.not. use_this_module) return
  
  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_rivermix_mod (runoff_calving_discharge_tracer): module must be initialized')
  endif 

  tau   = Time%tau
  
  do n=1,num_prog_tracers  
    T_prog(n)%wrk1(:,:,:) = 0.0
  enddo  
  delta             = 0.0
  delta_rho_tocean  = 0.0
  delta_rho0_triver = 0.0
  wrk1              = 0.0
  tracer_flux       = 0.0
               
  do j=jsc,jec
     do i=isc,iec

        if (river(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then


        ! the array "river" contains the volume rate (m/s) or mass
        ! rate (kg/m2/sec) of fluid with tracer 
        ! that is to be distributed in the vertical. 

            depth = min(Grd%ht(i,j),insertion_thickness)              ! be sure not to discharge river content into rock 
            nz    = min(Grd%kmt(i,j),floor(frac_index(depth,Grd%zw))) ! number of k-levels into which discharge rivers
            nz    = max(1,nz)                                         ! make sure have at least one cell to discharge into

            ! determine fractional thicknesses of grid cells 
            thkocean = 0.0
            do k=1,nz
               thkocean = thkocean + Thickness%rho_dzt(i,j,k,tau)
            enddo
            do k=1,nz
               delta(k) = Thickness%rho_dzt(i,j,k,tau)/(epsln+thkocean)
            enddo


            do n=1,num_prog_tracers

! mix a portion of the river mass(volume) inserted over the course of a tracer timestep
! into each model level from k=1:nz, in proportion to the fractional cell thickness.
! The mixing scheme used here follows scheme devised by Keith Dixon to handle hydraulic
! control across xland mixing points
!
! The following schematic illustrates the scheme, assuming constant layer thickness
! neglecting free surface height variations and assuming the river insertion depth
! extends to the base of the third model level.  Each of these assumptions are for 
! purposes of this schematic only.  The scheme that is implemented works in general
! for non-constant cell thicknesses.  
!               
!                  Total depth = H; cell thickness = dz; River mass flux/timestep = R
!  x--------x
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of 
!  !        !        Triver + 2*zinsert meters of water at T*(2)
!  !        !        T*(1) = (T(1).dz + 3.T*(2).zinsert + Triver.zinsert ) / (dz + 3.zinsert)               
!  x--------x               
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of 
!  !        !        Triver + zinsert meters of water at T*(3)           
!  !        !        T*(2) = (T(2).dz + 2.T*(3).zinsert + Triver.zinsert ) / (dz + 2.zinsert)               
!  !        !
!  x--------x                              
!  !        !
!  !        !  < === Insert  zinsert = R*dz/H meters of riverwater at a temperature of Triver 
!  !        !        T*(3) = (T(3).dz + Triver.zinsert ) / (dz + zinsert)
!  !        !
!  xxxxxxxxxx               

               if(runoff_type==1) then 
                  tracer_flux(i,j) = T_prog(n)%runoff_tracer_flux(i,j)
               else 
                  tracer_flux(i,j) = T_prog(n)%calving_tracer_flux(i,j)
               endif

               zextra=0.0
               do k=nz,1,-1
                  tracernew(k) = 0.0

                  if (k.eq.nz) then
                      tracerextra=0.0
                  else
                      tracerextra = tracernew(k+1)
                  endif

                  zinsert = river(i,j)*dtime*delta(k)
                  tracernew(k) = (tracerextra*zextra + T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau) + &
                                  tracer_flux(i,j)*dtime*delta(k)) / (zextra+Thickness%rho_dzt(i,j,k,tau)+zinsert)

                  zextra=zextra+zinsert
               enddo

               k=1
               T_prog(n)%wrk1(i,j,k) = (tracernew(k)*(Thickness%rho_dzt(i,j,k,tau)+river(i,j)*dtime) -&
                                        T_prog(n)%field(i,j,k,tau)*Thickness%rho_dzt(i,j,k,tau))/dtime
               do k=2,nz
                  T_prog(n)%wrk1(i,j,k) = Thickness%rho_dzt(i,j,k,tau)*(tracernew(k) - T_prog(n)%field(i,j,k,tau))/dtime
               enddo

               if(debug_all_in_top_cell) then 
                   k=1
                   T_prog(n)%wrk1(i,j,:) = 0.0
                   T_prog(n)%wrk1(i,j,k) = Grd%tmask(i,j,k)*river(i,j)*T_prog(n)%triver(i,j)
               endif

               do k=1,nz
                  T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
               enddo

               if (debug_this_module) then
                   write(6,*) 'i,j,n,river((kg/m^3)*m/sec),tracer_flux= ',i,j,n,river(i,j),tracer_flux(i,j)
                   do k=1,nz
                      write(6,*) 'i,j,k,tracer(k),tracernew(k),tendency(tracer*rho*thickness/sec)= ',&
                           i,j,k,T_prog(n)%field(i,j,k,tau),tracernew(k),T_prog(n)%wrk1(i,j,k)
                   enddo
               endif
                            
               
            enddo ! n end for number of prognostic tracers 

        endif ! river > 0

     enddo   ! i end
  enddo      ! j end



  if (debug_this_module_heat) then 

      if(runoff_type==1) then 
          wrk1_2d(:,:) = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*T_prog(index_temp)%conversion &
                              *T_prog(index_temp)%runoff_tracer_flux(i,j)*dtime
             enddo
          enddo
          tracer_input = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), BITWISE_EXACT_SUM)
          write(stdoutunit,'(a)')
          write (stdoutunit,'(a,es24.17,a)') ' Heat input via runoff in rivermix  = ',&
               tracer_input,' Joule'
          write(stdoutunit,'(a)')
      endif

      if(runoff_type==2) then 
          wrk1_2d(:,:) = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*T_prog(index_temp)%conversion &
                              *T_prog(index_temp)%calving_tracer_flux(i,j)*dtime
             enddo
          enddo
          tracer_input = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), BITWISE_EXACT_SUM)
          write(stdoutunit,'(a)')
          write (stdoutunit,'(a,es24.17,a)') ' Heat input via calving in rivermix  = ',&
               tracer_input,' Joule'
          write(stdoutunit,'(a)')
      endif

  endif

end subroutine runoff_calving_discharge_tracer
! </SUBROUTINE> NAME="runoff_calving_discharge_tracer"



!#######################################################################
! <SUBROUTINE NAME="river_kappa">
!
! <DESCRIPTION>
! This subroutine enhances the vertical diffusivity kappa over 
! a vertical column whose thickness is set by river_diffusion_thickness 
! and whose horizontal location is given by the rmask array.
! Note that rmask can be > 0 even if river=0 in the case when 
! use virtual salt flux.   
! The enhanced diffusivity is maximum at the top cell and is linearly 
! interpolated to the normal diffusivity at the depth set by 
! river_diffusion_thickness
! </DESCRIPTION>
!
subroutine river_kappa (Time, Thickness, Tracer, kappa)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_prog_tracer_type),  intent(in)    :: Tracer
  real, dimension(isd:,jsd:,:),  intent(inout) :: kappa

  integer :: i, j, k, nz
  real    :: depth 
  real, dimension(nk) :: zw_ij

  if(.not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_rivermix_mod (river_kappa): module must be initialized')
  endif 

  wrk1=0.0

  do j=jsc,jec
     do i=isc,iec

        do k=1,nk
           zw_ij(k) = Thickness%depth_zwt(i,j,k) 
        enddo

        if (Tracer%riverdiffuse(i,j) > 0.0 .and. Grd%kmt(i,j) > 0) then
            depth = min(Grd%ht(i,j),river_diffusion_thickness)
            nz    = min(Grd%kmt(i,j),floor(frac_index(depth,zw_ij)))
            nz    = max(1,nz)                                
            do k=1,nz 
              wrk1(i,j,k)  = river_diffusivity*(1.0 - Thickness%depth_zt(i,j,k)/Thickness%depth_zt(i,j,nk))
              kappa(i,j,k) = kappa(i,j,k) + wrk1(i,j,k)
            enddo
        endif
     enddo
  enddo

  if (id_diff_cbt_river_t > 0 .and. Tracer%name=='temp') then 
      used = send_data (id_diff_cbt_river_t, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_diff_cbt_river_s > 0 .and. Tracer%name=='salt') then 
      used = send_data (id_diff_cbt_river_s, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 

end subroutine river_kappa
! </SUBROUTINE> NAME="river_kappa"


!#######################################################################
! <SUBROUTINE NAME="river_diagnostics">
!
! <DESCRIPTION>
! For sending some diagnostics. 
! </DESCRIPTION>
!
subroutine  river_diagnostics (Time, Thickness, Dens, temp_wrk, salt_wrk, mix_type)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_thickness_type),   intent(in)  :: Thickness
  type(ocean_density_type),     intent(in)  :: Dens
  real, dimension(isd:,jsd:,:), intent(in)  :: temp_wrk
  real, dimension(isd:,jsd:,:), intent(in)  :: salt_wrk
  integer,                      intent(in)  :: mix_type 

  integer :: i,j,k,n,tau
  real    :: temporary

  ! runoff
  if(mix_type==1) then 

      if (id_neutral_rho_runoffmix > 0) then
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = Grd%tmask(i,j,k)            &  
                        *(Dens%drhodT(i,j,k)*temp_wrk(i,j,k) &
                         +Dens%drhodS(i,j,k)*salt_wrk(i,j,k))
                enddo
             enddo
          enddo
          used = send_data (id_neutral_rho_runoffmix, wrk1(:,:,:), &
               Time%model_time,rmask=Grd%tmask(:,:,:),             &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif

      if (id_wdian_rho_runoffmix > 0) then 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   temporary   = Dens%drhodz_zt(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
                   wrk1(i,j,k) = Grd%tmask(i,j,k)             &  
                        *(Dens%drhodT(i,j,k)*temp_wrk(i,j,k)  &
                         +Dens%drhodS(i,j,k)*salt_wrk(i,j,k)) &
                        /(epsln+temporary)
                enddo
             enddo
          enddo
          used = send_data (id_wdian_rho_runoffmix, wrk1(:,:,:),&
               Time%model_time,rmask=Grd%tmask(:,:,:),          &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif


  ! calving
  elseif(mix_type==2) then 

      if (id_neutral_rho_calvingmix > 0) then
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = Grd%tmask(i,j,k)            &  
                        *(Dens%drhodT(i,j,k)*temp_wrk(i,j,k) &
                         +Dens%drhodS(i,j,k)*salt_wrk(i,j,k))
                enddo
             enddo
          enddo
          used = send_data (id_neutral_rho_calvingmix, wrk1(:,:,:),&
               Time%model_time,rmask=Grd%tmask(:,:,:),             &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif

      if (id_wdian_rho_calvingmix > 0) then 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   temporary   = Dens%drhodz_zt(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
                   wrk1(i,j,k) = Grd%tmask(i,j,k)             &  
                        *(Dens%drhodT(i,j,k)*temp_wrk(i,j,k)  &
                         +Dens%drhodS(i,j,k)*salt_wrk(i,j,k)) &
                        /(epsln+temporary)
                enddo
             enddo
          enddo
          used = send_data (id_wdian_rho_calvingmix, wrk1(:,:,:),&
               Time%model_time,rmask=Grd%tmask(:,:,:),           &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif


  ! river=runoff+calving
  elseif(mix_type==3) then 

      if (id_neutral_rho_rivermix > 0) then
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   wrk1(i,j,k) = Grd%tmask(i,j,k)            &  
                        *(Dens%drhodT(i,j,k)*temp_wrk(i,j,k) &
                         +Dens%drhodS(i,j,k)*salt_wrk(i,j,k))
                enddo
             enddo
          enddo
          used = send_data (id_neutral_rho_rivermix, wrk1(:,:,:),&
               Time%model_time,rmask=Grd%tmask(:,:,:),           &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif

      if (id_wdian_rho_rivermix > 0) then 
          wrk1(:,:,:) = 0.0
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   temporary   = Dens%drhodz_zt(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
                   wrk1(i,j,k) = Grd%tmask(i,j,k)             &  
                        *(Dens%drhodT(i,j,k)*temp_wrk(i,j,k)  &
                         +Dens%drhodS(i,j,k)*salt_wrk(i,j,k)) &
                         /(epsln+temporary)
                enddo
             enddo
          enddo
          used = send_data (id_wdian_rho_rivermix, wrk1(:,:,:), &
               Time%model_time,rmask=Grd%tmask(:,:,:),          &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      endif


  endif 

end subroutine river_diagnostics
! </SUBROUTINE> NAME="river_diagnostics"



end module ocean_rivermix_mod
      
      




