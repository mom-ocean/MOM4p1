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
module ocean_convect_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<CONTACT EMAIL="russell.fiedler@csiro.au"> R. Fiedler 
!</CONTACT>
!
!<OVERVIEW>
! Vertically adjusts gravitationally unstable columns of ocean fluid.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module vertically adjusts gravitationally unstable columns of 
! ocean fluid.  Three algorithms are available:
!
! 1. Full convection from Rahmstorf.  The algorithm produces a fully 
! stable fluid column.  Since most convection propagates downward, the
! scheme looks downward first and follows any instability (upward or
! downward) before checking the other direction. The routine mixes 
! passive tracers only after the entire instability is found. 
!
! 2. Full convection from Rahmstorf as optimized for vector machines
! by Russ Fiedler. 
!
! 3. The Cox (1984) NCON-scheme.  This scheme is recommended only 
! for those wishing to maintain legacy code. 
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Stefan Rahmstorf (Ocean Modelling, 1993 vol 101 pages 9-11)
! </REFERENCE>
!
! <NOTE>
! Implementation of the full convection scheme is 
! based on mom2/3 code by Stefan Rahmstorf 
! (rahmstorf@pik-potsdam.de). But  modified slightly 
! for efficiency purposes in mom3.1
! by M. Eby (eby@uvic.ca) in June 2000. Notably, Eby 
! eliminated goto statements.     
! </NOTE>
!
! <NOTE>
! The Eby code was ported to mom4 by Griffies (Stephen.Griffies@noaa.gov). 
! </NOTE>
!
! <NOTE>
! To recover the exact same numerical values as the original
! Rahmstorf code, look for the two "Rahmstorf" comments in the code.  
! </NOTE>
!
! </INFO>
!<NAMELIST NAME="ocean_convect_nml">
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false. 
!  </DATA> 
!
!  <DATA NAME="convect_ncon" TYPE="logical">
!  If true, will use the old NCON convection scheme from Cox.  Recommended only 
!  if wishing to reproduce old results. 
!  </DATA> 
!  <DATA NAME="ncon" TYPE="integer">
!  Number of passes through the NCON-scheme. 
!  </DATA> 
!  <DATA NAME="convect_full_scalar" TYPE="logical">
!  If true, will use the full convection scheme as implemented at GFDL for scalar
!  machines. 
!  </DATA> 
!  <DATA NAME="convect_full_vector" TYPE="logical">
!  If true, will use the full convection scheme as optimized for vector machines
!  by Russ Fiedler. 
!  </DATA> 
!</NAMELIST>
!
use diag_manager_mod,    only: register_diag_field, send_data
use fms_mod,             only: open_namelist_file, check_nml_error
use fms_mod,             only: write_version_number, close_file
use fms_mod,             only: mpp_error, FATAL, stdout, stdlog

use ocean_density_mod,    only: density
use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: missing_value
use ocean_types_mod,      only: ocean_time_type, ocean_grid_type, ocean_domain_type, ocean_options_type
use ocean_types_mod,      only: ocean_density_type, ocean_prog_tracer_type, ocean_thickness_type
use ocean_workspace_mod,  only: wrk1


implicit none

private

public ocean_convect_init
public convection

private convection_ncon
private convection_full_scalar 
private convection_full_vector 

! for diagnostics 
integer :: id_ktot=-1
integer :: id_kven=-1
integer, allocatable, dimension(:) :: id_convect
logical :: used

integer :: nkm1         
integer :: index_temp=-1
integer :: index_salt=-1
real    :: dtimer

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

character(len=128) :: version=&
  '$Id: ocean_convect.F90,v 16.0.108.1 2009/10/10 00:42:17 nnz Exp $'
character (len=128) :: tagname = &
  '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: module_is_initialized=.FALSE.
logical :: use_this_module      =.false.
logical :: convect_ncon         =.false.  ! for Cox ncon scheme. if false, use Rahmstorf full convection 
logical :: convect_full_scalar  =.false.  ! for full convection scheme implemented at GFDL on scalar machines
logical :: convect_full_vector  =.false.  ! for full convection scheme implemented at CSIRO on vector machines
integer :: ncon=7                         ! number of passes through a column for convect_ncon=.true.   

namelist /ocean_convect_nml/ use_this_module, convect_ncon, ncon, convect_full_scalar, convect_full_vector 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_convect_init">
!
! <DESCRIPTION>
!
! Initialize the convection module.
! 
! For the full convection module, we register two fields
! for diagnostic output.  
! <BR/>
! ktot  = number of levels convected  in a vertical column
! <BR/> 
! kven  = number of levels ventilated in a vertical column
!
! Note that ktot can in rare cases count some levels twice, if they
! get involved in two originally separate, but then
! overlapping convection areas in the water column. The field 
! kven is 0 on land, 1 on ocean points with no
! convection, and any value up to nk on convecting points. 
! </DESCRIPTION>
!
subroutine ocean_convect_init (Grid, Domain, Time, T_prog, Ocean_options, dtime)

  type(ocean_grid_type),        intent(in), target :: Grid
  type(ocean_domain_type),      intent(in), target :: Domain
  type(ocean_time_type),        intent(in)         :: Time
  type(ocean_prog_tracer_type), intent(inout)      :: T_prog(:)
  type(ocean_options_type),     intent(inout)      :: Ocean_options
  real,                         intent(in)         :: dtime

  integer :: n, num_prog_tracers, ioun, ierr, io_status

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_convect_mod (ocean_convect_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  num_prog_tracers = size(T_prog(:))

  allocate( id_convect(num_prog_tracers) )
  id_convect(:) = -1

  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  Grd => Grid
  Dom => Domain

  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grd%nk
  nkm1 = nk-1

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_convect_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_convect_nml)  
  write (stdlogunit,ocean_convect_nml)
  ierr = check_nml_error(io_status,'ocean_convect_nml')
  call close_file (ioun)

  if(use_this_module) then 
     write(stdoutunit,'(a)') &
     '==>Note: Using convective adjustment to stabilize gravitationally unstable water columns.'
     Ocean_options%convective_adjustment = 'Used convective adjustment.'
  else 
     write(stdoutunit,'(a)') &
     '==>Note: NOT using convective adjustment in gravitationally unstable water columns.'
     Ocean_options%convective_adjustment = 'Did NOT use convective adjustment.'
     return 
  endif 

  if(convect_ncon) then 
      write(stdoutunit,'(a,i3)')' ==>Note: Using the Cox NCON convection scheme with ncon= ', ncon
  elseif(convect_full_vector .or. convect_full_scalar) then  
      write(stdoutunit,'(a)')' ==>Note: Using the Rahmstorf full convection scheme'
      if(convect_full_scalar .and. convect_full_vector) then 
          call mpp_error(FATAL,&
          '==>Error in ocean_convect_mod: "convect_full_vector" and "convect_full_scalar" cannot both be true.')
      endif
      if(convect_full_vector) then 
          write(stdoutunit,'(a)')'   as implemented on vector machines at CSIRO/Hobart.'
      endif
      if(convect_full_scalar) then 
          write(stdoutunit,'(a)')'   as implemented on scalar machines at GFDL.'
      endif
  else
      write(stdoutunit,'(a)')&
      '==>Error in ocean_convect_mod: no convective adjustment selected, yet use_this_module=.true.. '
      write(stdoutunit,'(a)')&
      '   Choose one of the following: "convect_ncon",  "convect_full_vector", or "convect_full_scalar"'
      call mpp_error(FATAL,&
      '==>ocean_convect_mod: convective adjustment not chosen.  Must choose a scheme since use_this_module=.true.')
  endif 
  

  ! register fields for diagnostic output
  dtimer = 1.0/dtime

  do n=1,num_prog_tracers 
     id_convect(n) = -1
     if (n == index_temp) then
         id_convect(n) = register_diag_field ('ocean_model', 'convect_heating',               &
                         Grid%tracer_axes(1:3), Time%model_time, 'heating due to convection', &
                         'Watts/m^2', missing_value=missing_value, range=(/-1.e10,1.e10/))
     else
         id_convect(n) = register_diag_field ('ocean_model', 'convect_'//trim(T_prog(n)%name),            &
                         Grid%tracer_axes(1:3), Time%model_time, 'rho*dzt*th_tendency due to convection', &
                        'm*kg/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
     endif
  enddo 

  id_ktot   = register_diag_field ('ocean_model', 'ktot', Grid%tracer_axes(1:2), Time%model_time, &
              'number convecting levels', 'number', missing_value=missing_value, range=(/-10.0,1.e20/))
  id_kven   = register_diag_field ('ocean_model', 'kven', Grid%tracer_axes(1:2), Time%model_time, &
              'number ventilated levels', 'number', missing_value=missing_value, range=(/-10.0,1.e20/))


end subroutine ocean_convect_init
! </SUBROUTINE>  NAME="ocean_convect_init"


!#######################################################################
! <SUBROUTINE NAME="convection">
!
! <DESCRIPTION>
! Subroutine calls one of the two possible convection schemes. 
! </DESCRIPTION>
!
subroutine convection (Time, Thickness, T_prog, Dens)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens

  if (.not. use_this_module) return 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error in ocean_convect_mod (convection): module needs to be initialized')
  endif 

  if(convect_ncon) then 
    call convection_ncon(Time, Thickness, T_prog, Dens)
  else 
    if(convect_full_scalar) call convection_full_scalar(Time, Thickness, T_prog, Dens)
    if(convect_full_vector) call convection_full_vector(Time, Thickness, T_prog, Dens)
  endif 

end subroutine convection 
! </SUBROUTINE> NAME="convection"


!#######################################################################
! <SUBROUTINE NAME="convection_full_scalar">
!
! <DESCRIPTION>
! Subroutine to vertically adjust gravitationally unstable columns of ocean fluid.
! Produces updated values for all the tracers.  Code implemented on scalar 
! machines at GFDL. Has been found to be slow on vector machines.  Use 
! convection_full_vector for vector machines.  
!
! internal variables:
!
!     chk_la = logical flag to check level above kt
!
!     chk_lb = logical flag to check level below kb
!
!     kb     = bottom level of (potential) instability
!
!     kbo    = bottom level of ocean
!
!     kt     = top level of (potential) instability
!
!     la     = test level above kt
!
!     lb     = test level below kb
!
!     rl     = lower level density referenced to lower level
!
!     ru     = upper level density referenced to lower level
!
!     tmx    = mixed tracer (1=temp, 2=salt, 3=other)
!
!     tsm    = sum of tracers (weighted by thickness) in the instability
!
!     zsm    = total thickness of the instability
!
! </DESCRIPTION>
!
subroutine convection_full_scalar (Time, Thickness, T_prog, Dens)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens

  real, dimension(3) :: tmx, tsm
  real, dimension(isc:iec,jsc:jec) :: ktot, kven

  integer :: i, j, k, n, kb, kbo, kt, la, lb, taup1, tau
  real :: rl, ru, zsm
  logical :: chk_la, chk_lb

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error in ocean_convect_mod (convection_full_scalar): module needs to be initialized')
  endif 

  ! for diagnostics, save initial value of the tracers 

  taup1 = Time%taup1
  tau   = Time%tau
  
  do n=1, size(T_prog(:))
     T_prog(n)%wrk1(:,:,:) = T_prog(n)%field(:,:,:,taup1)    
  enddo

  do j=jsc,jec
    do i=isc,iec
      kbo       = Grd%kmt(i,j)
      ktot(i,j) = 0.0
      kven(i,j) = 0.0

      ! search for unstable regions starting from the top
      kt = 1
      kb = 2
      do while (kt < kbo)

        ru = density (T_prog(index_salt)%field(i,j,kt,taup1),&
                      T_prog(index_temp)%field(i,j,kt,taup1),&
                      Dens%pressure_at_depth(i,j,kb))
        rl = density (T_prog(index_salt)%field(i,j,kb,taup1),&
                      T_prog(index_temp)%field(i,j,kb,taup1),&
                      Dens%pressure_at_depth(i,j,kb))

        ! sum the first pair found in an unstable region
        if (ru > rl) then
          chk_la = .true.
          chk_lb = .true.
          zsm    = Thickness%rho_dzt(i,j,kt,taup1) + Thickness%rho_dzt(i,j,kb,taup1)
          tsm(1) = T_prog(index_temp)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1) &
                 + T_prog(index_temp)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
          tmx(1) = tsm(1)/zsm
          tsm(2) = T_prog(index_salt)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1) &
                 + T_prog(index_salt)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
          tmx(2) = tsm(2)/zsm

          do while (chk_lb .or. chk_la)

            ! check for an unstable level (lb) below kb
            if (kb >= kbo) chk_lb = .false.
            do while (chk_lb)
              chk_lb = .false.
              lb = kb + 1 
              ru = density (tmx(2)     ,tmx(1)     ,Dens%pressure_at_depth(i,j,lb))
              rl = density (T_prog(index_salt)%field(i,j,lb,taup1),&
                            T_prog(index_temp)%field(i,j,lb,taup1),&
                            Dens%pressure_at_depth(i,j,lb))
              if (ru > rl) then
                ! add new level to sums
                kb     = lb
                zsm    = zsm + Thickness%rho_dzt(i,j,kb,taup1)
                tsm(1) = tsm(1) + T_prog(index_temp)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
                tmx(1) = tsm(1)/zsm
                tsm(2) = tsm(2) + T_prog(index_salt)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
                tmx(2) = tsm(2)/zsm
                chk_la = .true.
                if (kb < kbo) chk_lb = .true.
              endif
            enddo

            ! check for an unstable level (la) above kt
            ! To get equivalent of original Rahmstorf code, uncomment next line
            ! chk_la = .true.
            if (kt <= 1) chk_la = .false.
            do while (chk_la)
              chk_la = .false.
              la = kt - 1
              ru = density (T_prog(index_salt)%field(i,j,la,taup1),&
                            T_prog(index_temp)%field(i,j,la,taup1),&
                            Dens%pressure_at_depth(i,j,kt))
              rl = density (tmx(2), tmx(1), Dens%pressure_at_depth(i,j,kt))
              if (ru > rl) then
                ! add new level to sums
                kt     = la
                zsm    = zsm + Thickness%rho_dzt(i,j,kt,taup1) 
                tsm(1) = tsm(1) + T_prog(index_temp)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1)
                tmx(1) = tsm(1)/zsm
                tsm(2) = tsm(2) + T_prog(index_salt)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1)
                tmx(2) = tsm(2)/zsm
                chk_lb = .true.
                ! to get equivalent of original Rahmstorf code, comment out next line
                if (kt > 1) chk_la = .true.
              endif
            enddo

          enddo

          ! mix all tracers from kt to kb
          do k=kt,kb
            T_prog(index_temp)%field(i,j,k,taup1) = tmx(1)
            T_prog(index_salt)%field(i,j,k,taup1) = tmx(2)
          enddo
          do n=1,size(T_prog(:))
            if (n /= index_temp .and. n /= index_salt) then
              tsm(3) = 0.0
              do k=kt,kb
                tsm(3) = tsm(3) + T_prog(n)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1)
              enddo
              tmx(3) = tsm(3)/zsm 
              do k=kt,kb
                T_prog(n)%field(i,j,k,taup1) = tmx(3)
              enddo
            endif
          enddo

          ! compute diagnostic fields
          ktot(i,j) = ktot(i,j) + kb - kt + 1
          if (kt == 1) kven(i,j) = kb

          kt = kb + 1

        else
          kt = kb
        endif

        ! continue the search for other unstable regions
        kb = kt + 1

      enddo   ! while (kt < kbo)

    enddo     ! i=isc,iec
  enddo       ! j=jsc,jec

  
  if (id_ktot > 0) used = send_data (id_ktot, ktot(isc:iec,jsc:jec), &
                          Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))
  if (id_kven > 0) used = send_data (id_kven, kven(isc:iec,jsc:jec), &
                          Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))

  do n=1,size(T_prog(:))
     if (id_convect(n) > 0) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 T_prog(n)%wrk1(i,j,k) =  dtimer*(T_prog(n)%field(i,j,k,taup1) - T_prog(n)%wrk1(i,j,k))
              enddo
           enddo
           T_prog(n)%wrk1(isc:iec,jsc:jec,k) = T_prog(n)%wrk1(isc:iec,jsc:jec,k)*&
                                               T_prog(n)%conversion*Thickness%rho_dzt(isc:iec,jsc:jec,k,taup1)
        enddo
        used = send_data(id_convect(n), T_prog(n)%wrk1(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
  enddo


end subroutine convection_full_scalar
! </SUBROUTINE> NAME="convection_full_scalar"



!#######################################################################
! <SUBROUTINE NAME="convection_full_vector">
!
! <DESCRIPTION>
! Subroutine to vertically adjust gravitationally unstable columns of ocean fluid.
! Produces updated values for all the tracers.  Code implemented on vector 
! machines at CSIRO. Has been found to be faster on these machines than 
! convection_full_scalar. Answers differ, but not significantly. 
!
! Code written by russell.fiedler@csiro.au 
! most recently modified Aug2005 
!
! internal variables:
!
!     chk_la = logical flag to check level above kt
!
!     chk_lb = logical flag to check level below kb
!
!     kb     = bottom level of (potential) instability
!
!     kbo    = bottom level of ocean
!
!     kt     = top level of (potential) instability
!
!     la     = test level above kt
!
!     lb     = test level below kb
!
!     rl     = lower level density referenced to lower level
!
!     ru     = upper level density referenced to lower level
!
!     tmx    = mixed tracer (1=temp, 2=salt, 3=other)
!
!     tsm    = sum of tracers (weighted by thickness) in the instability
!
!     zsm    = total thickness of the instability
!
! </DESCRIPTION>
!
subroutine convection_full_vector (Time, Thickness, T_prog, Dens)
  
  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens

  real, dimension(3)               :: tmx, tsm
  real, dimension(isc:iec,jsc:jec) :: ktot, kven

  integer :: i, j, k, n, kb, kbo, kt, la, lb, taup1, tau
  real    :: rl, ru, zsm
  logical :: chk_la, chk_lb

  ! additions for vectorised code
  logical :: unstable(isc:iec)
  integer :: ikven(isc:iec,jsc:jec), kt_temp(isc:iec)
  real    :: zsum(isc:iec),zsum1(isc:iec)
  real    :: rl_line(isd:ied), ru_line(isd:ied)
  real    :: ratio1(isc:iec),ratio2(isc:iec)
  integer :: ikunstable, nunstable
  integer,dimension((iec-isc+1)*nk)  :: iindex,kindex 


  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error in ocean_convect_mod (convection_full_vector): module needs to be initialized')
  endif 

  ! for diagnostics, save initial value of the tracers 

  taup1 = Time%taup1
  tau   = Time%tau
  
  do n=1, size(T_prog(:))
    T_prog(n)%wrk1(:,:,:) = T_prog(n)%field(:,:,:,taup1)    
  enddo

!
! First we treat instabilities arising at the surface.
!

  do j=jsc,jec
!
! Initialise instability flag. Mask land points.
!
    do i=isc,iec
      if(Grd%kmt(i,j) .ne. 0) then
        unstable(i)=.true.
      else
        unstable(i)=.false.
      endif
    
      ikven(i,j)=0
      zsum1(i)=Thickness%rho_dzt(i,j,1,taup1)
    enddo
    do k=1,nk-1
!
! Calculate depth of this cell and the one below it; If we hit land column has
! been stabilised.
!
      do i=isc,iec
        zsum(i)=zsum1(i)
        zsum1(i)=zsum1(i)+Thickness%rho_dzt(i,j,k+1,taup1)
        ratio1(i)=zsum(i)/zsum1(i)
        ratio2(i)=1.0-ratio1(i)
        if(k >= Grd%kmt(i,j)) then
          unstable(i) = .false.
        endif
      enddo
!
! Check if entire row has been stabilised.
!
      if(.not.any(unstable(:))) exit
!
! Main part.
! calculate densities for this row and the one below it
! If unstable, mix the layers and mark the depth to which convection has
! occured.
!
      ru_line = density(T_prog(index_salt)%field(:,j,k,taup1),&
                        T_prog(index_temp)%field(:,j,k,taup1),&
                        Dens%pressure_at_depth(:,j,k+1))
      rl_line = density(T_prog(index_salt)%field(:,j,k+1,taup1),&
                        T_prog(index_temp)%field(:,j,k+1,taup1),&
                        Dens%pressure_at_depth(:,j,k+1))

      do n=1,size(T_prog(:))
        do i=isc,iec
          if(unstable(i)) then
            if(ru_line(i) > rl_line(i)) then
              T_prog(n)%field(i,j,k+1,taup1) = T_prog(n)%field(i,j,k+1,taup1)*ratio2(i) + &
                                               T_prog(n)%field(i,j,k,taup1)*ratio1(i)
              ikven(i,j)=k+1
            else
              unstable(i)=.false.
            endif
          endif
        enddo
      enddo
      if(.not. any(unstable(:))) exit
    enddo

!
! Now adjust all the layerabove the bottom of the instablity
!
    do n=1,size(T_prog(:))
      do k=1,maxval(ikven(:,j))-1
        do i=isc,iec
          if(k <= ikven(i,j)-1) then
            T_prog(n)%field(i,j,k,taup1)=T_prog(n)%field(i,j,ikven(i,j),taup1)
          endif
        enddo
      enddo
    enddo
  enddo

!
! Now procede to treating instabilies arising within the water column.
!
  ktot=0.0
  kven=0.0

  jloop: do j=jsc,jec
!
! Preliminary sweep through density field to eliminate stable columns 
! Note that top level is necessarily stable from initial sweep.
! 
   kt_temp=nk
   nunstable=0
   iindex=0
   kindex=0     
   kloop: do k=2,maxval(Grd%kmt(isc:iec,j))-1
      ru_line = density(T_prog(index_salt)%field(:,j,k,taup1),&
                        T_prog(index_temp)%field(:,j,k,taup1),&
                        Dens%pressure_at_depth(:,j,k+1))
      rl_line = density(T_prog(index_salt)%field(:,j,k+1,taup1),&
                        T_prog(index_temp)%field(:,j,k+1,taup1),&
                        Dens%pressure_at_depth(:,j,k+1))
!
! If an instability occurs anywhere at this level exit to the main loop
! and process as normal.
!
      do i=isc,iec
        if(k < Grd%kmt(i,j) .and. ru_line(i) > rl_line(i)) then
          nunstable=nunstable+1                                    
          iindex(nunstable)=i                                      
          kindex(nunstable)=k                                      
        endif
      enddo
    enddo kloop
    if(nunstable .eq. 0) cycle jloop 

    do ikunstable=1,nunstable                                      
      i   = iindex(ikunstable)                                         
      kt  = kindex(ikunstable)                                        
      kbo = Grd%kmt(i,j)                                     
      kb  = kt+1                                                    

        ru = density (T_prog(index_salt)%field(i,j,kt,taup1),&
                      T_prog(index_temp)%field(i,j,kt,taup1),&
                      Dens%pressure_at_depth(i,j,kb))
        rl = density (T_prog(index_salt)%field(i,j,kb,taup1),&
                      T_prog(index_temp)%field(i,j,kb,taup1),&
                      Dens%pressure_at_depth(i,j,kb))

        ! sum the first pair found in an unstable region
        if (ru > rl) then
          chk_la = .true.
          chk_lb = .true.
          zsm    = Thickness%rho_dzt(i,j,kt,taup1) + Thickness%rho_dzt(i,j,kb,taup1)
          tsm(1) =   T_prog(index_temp)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1) &
                   + T_prog(index_temp)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
          tmx(1) = tsm(1)/zsm
          tsm(2) =   T_prog(index_salt)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1)  &
                   + T_prog(index_salt)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
          tmx(2) = tsm(2)/zsm

          do while (chk_lb .or. chk_la)

            ! check for an unstable level (lb) below kb
            if (kb >= kbo) chk_lb = .false.
            do while (chk_lb)
              chk_lb = .false.
              lb = kb + 1 
              ru = density (tmx(2), tmx(1), Dens%pressure_at_depth(i,j,lb))
              rl = density (T_prog(index_salt)%field(i,j,lb,taup1),&
                            T_prog(index_temp)%field(i,j,lb,taup1),&
                            Dens%pressure_at_depth(i,j,lb))
              if (ru > rl) then
                ! add new level to sums
                kb = lb
                zsm = zsm + Thickness%rho_dzt(i,j,kb,taup1)
                tsm(1) = tsm(1) + T_prog(index_temp)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
                tmx(1) = tsm(1)/zsm
                tsm(2) = tsm(2) + T_prog(index_salt)%field(i,j,kb,taup1)*Thickness%rho_dzt(i,j,kb,taup1)
                tmx(2) = tsm(2)/zsm
                chk_la = .true.
                if (kb < kbo) chk_lb = .true.
              endif
            enddo

            ! check for an unstable level (la) above kt
            ! To get equivalent of original Rahmstorf code, uncomment next line
            ! chk_la = .true.
            if (kt <= 1) chk_la = .false.
            do while (chk_la)
              chk_la = .false.
              la = kt - 1
              ru = density (T_prog(index_salt)%field(i,j,la,taup1),&
                            T_prog(index_temp)%field(i,j,la,taup1),&
                            Dens%pressure_at_depth(i,j,kt))
              rl = density (tmx(2), tmx(1), Dens%pressure_at_depth(i,j,kt))
              if (ru > rl) then
                ! add new level to sums
                kt = la
                zsm = zsm + Thickness%rho_dzt(i,j,kt,taup1) 
                tsm(1) = tsm(1) + T_prog(index_temp)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1)
                tmx(1) = tsm(1)/zsm
                tsm(2) = tsm(2) + T_prog(index_salt)%field(i,j,kt,taup1)*Thickness%rho_dzt(i,j,kt,taup1)
                tmx(2) = tsm(2)/zsm
                chk_lb = .true.
                ! to get equivalent of original Rahmstorf code, comment out next line
                if (kt > 1) chk_la = .true.
              endif
            enddo

          enddo

          ! mix all tracers from kt to kb
          do k=kt,kb
            T_prog(index_temp)%field(i,j,k,taup1) = tmx(1)
            T_prog(index_salt)%field(i,j,k,taup1) = tmx(2)
          enddo
          do n=1,size(T_prog(:))
            if (n /= index_temp .and. n /= index_salt) then
              tsm(3) = 0.0
              do k=kt,kb
                tsm(3) = tsm(3) + T_prog(n)%field(i,j,k,taup1)*Thickness%rho_dzt(i,j,k,taup1)
              enddo
              tmx(3) = tsm(3)/zsm 
              do k=kt,kb
                T_prog(n)%field(i,j,k,taup1) = tmx(3)
              enddo
            endif
          enddo

          ! compute diagnostic fields
          ktot(i,j) = ktot(i,j) + kb - kt + 1
          if (kt == 1) kven(i,j) = kb


        endif

    enddo         !ikunstable 
  enddo jloop     ! j=jsc,jec


! Merge diagnostic info.

  do j=jsc,jec
    do i=isc,iec
      if(real(ikven(i,j)) > kven(i,j)) then
        kven(i,j)=ikven(i,j)
        ktot(i,j)=ikven(i,j)
      endif
    enddo
  enddo
  
  if (id_ktot > 0) used = send_data (id_ktot, ktot(isc:iec,jsc:jec), &
                          Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))
  if (id_kven > 0) used = send_data (id_kven, kven(isc:iec,jsc:jec), &
                          Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1))

  do n=1,size(T_prog(:))
     if (id_convect(n) > 0) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 T_prog(n)%wrk1(i,j,k) =  dtimer*(T_prog(n)%field(i,j,k,taup1) - T_prog(n)%wrk1(i,j,k))
              enddo
           enddo
           T_prog(n)%wrk1(isc:iec,jsc:jec,k) = T_prog(n)%wrk1(isc:iec,jsc:jec,k)*&
            T_prog(n)%conversion*Thickness%rho_dzt(isc:iec,jsc:jec,k,taup1)
        enddo
        used = send_data(id_convect(n), T_prog(n)%wrk1(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
  enddo


end subroutine convection_full_vector
! </SUBROUTINE> NAME="convection_full_vector"


!#######################################################################
! <SUBROUTINE NAME="convection_ncon">
!
! <DESCRIPTION>
! "ncon" convection scheme
!
! Convectively adjust water column if gravitationally unstable.
! Based on algorithm from Mike Cox used in his code from 1984.  
! Algorithm has well known problems with incomplete homogenization
! and sensitivity to the ncon parameter. 
!
! Coded in mom4 for legacy purposes by Stephen.Griffies@noaa.gov
! April 2001  
!
! </DESCRIPTION>
!
subroutine convection_ncon (Time, Thickness, T_prog, Dens)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_density_type),     intent(in)    :: Dens

  integer :: i, j, k, l, ks, n, nn
  integer :: tau, taup1
  real    :: dense 

  taup1 = Time%taup1
  tau   = Time%tau

  do n=1, size(T_prog(:))
     T_Prog(n)%wrk1(:,:,:) = T_prog(n)%field(:,:,:,taup1)    
  enddo
  wrk1 = 0.0

  ! ks=1: compare levels 1 to 2; 3 to 4; etc.
  ! ks=2: compare levels 2 to 3; 4 to 5; etc.
  do nn=1,ncon
    do ks=1,2

      ! find density with reference pressure determined by vertical pairs 
      do j=jsc,jec
        do l=1,nk,2
          if (ks .eq. 1) then
            k = min(l+1,nk)
          else
            k = l
          endif
          do i=isc,iec
            wrk1(i,j,l)=density(T_prog(index_salt)%field(i,j,l,taup1),&
                                T_prog(index_temp)%field(i,j,l,taup1),&
                                Dens%pressure_at_depth(i,j,k))
          enddo
        enddo
      enddo
      do j=jsc,jec
        do l=2,nk,2
          if (ks .eq. 1) then
            k = l
          else
            k = min(l+1,nk)
          endif
          do i=isc,iec
            wrk1(i,j,l)=density(T_prog(index_salt)%field(i,j,l,taup1),&
                                T_prog(index_temp)%field(i,j,l,taup1),&
                                Dens%pressure_at_depth(i,j,k))
          enddo
        enddo
      enddo

      ! set "heavy water" in land to stop convection
      dense = 1.e30
      do j=jsc,jec
        do i=isc,iec
          k = Grd%kmt(i,j) + 1
          if (k .le. nk) wrk1(i,j,k) = dense
        enddo
      enddo

      !if unstable, then mix tracers on adjoining levels
      do n=1,size(T_prog(:))
        do j=jsc,jec
          do k=ks,nkm1,2
            do i=isc,iec
              if (wrk1(i,j,k) .gt. wrk1(i,j,k+1)) then
                T_prog(n)%field(i,j,k,taup1) = (  Thickness%rho_dzt(i,j,k,taup1)*T_prog(n)%field(i,j,k,taup1)      &
                                                 +Thickness%rho_dzt(i,j,k+1,taup1)*T_prog(n)%field(i,j,k+1,taup1)) &
                                                /(Thickness%rho_dzt(i,j,k,taup1) + Thickness%rho_dzt(i,j,k+1,taup1))
                T_prog(n)%field(i,j,k+1,taup1) = T_prog(n)%field(i,j,k,taup1)
              endif
            enddo
          enddo
        enddo
      enddo

    enddo  !ks=1,2
  enddo    !nn=1,ncon

  do n=1,size(T_prog(:))
     if (id_convect(n) > 0) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 T_prog(n)%wrk1(i,j,k) =  dtimer*(T_prog(n)%field(i,j,k,taup1) - T_prog(n)%wrk1(i,j,k))
              enddo
           enddo
           T_prog(n)%wrk1(isc:iec,jsc:jec,k) = T_prog(n)%wrk1(isc:iec,jsc:jec,k)*&
            T_prog(n)%conversion*Thickness%rho_dzt(isc:iec,jsc:jec,k,taup1)
        enddo
        used = send_data(id_convect(n), T_prog(n)%wrk1(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
  enddo


end subroutine convection_ncon
! </SUBROUTINE> NAME="convection_ncon"



end module ocean_convect_mod 
