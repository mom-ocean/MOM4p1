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
module ocean_polar_filter_mod
!
!<CONTACT EMAIL="Mike.Spelman@noaa.gov"> Mike Spelman 
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Provide polar filtering of fields for use with spherical grid. 
! Set up only for filtering tracers with 1D domain decomposition. 
! This code is no longer supported by GFDL scientists.  
!</OVERVIEW>
!
!<DESCRIPTION>
! This module provides polar filtering of tracers for use 
! with spherical coordinate grids.  Should not be used with 
! non-spherical grids. This code is provided for legacy purposes
! to allow modelers the opportunity to test their older
! spherical models in mom4 prior to moving to a generalized 
! horizontal grid, such as the tripolar. 
!
! Polar filtering has many well known problems, especially 
! when filtering the velocity and/or free surface fields. 
! Hence, this module only provides for filtering the tracer
! field.  Even so, its use is discouraged for those building 
! new ocean model configurations.  Most notably, this scheme
! DOES NOT conserve total tracer.
!
! This scheme has been implemented ONLY for cases with 1D domain
! decomposition (constant latitude rows).  Paralellization in 
! 2D is not available with this implementation. 
!
! There are two methods for polar filtering: (1) polar filter 
! the time tendencies (as in mom3) and (2) polar filtering the 
! fields themselves (as in mom1).  The mom1 method is preferred
! at GFDL when running with an ice model.  
!
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_polar_filter_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!
!  </DATA> 
!  <DATA NAME="polar_filter_field" TYPE="logical">
!  Polar filter the tracer field (default method)
!  </DATA> 
!  <DATA NAME="polar_filter_tendency" TYPE="logical">
!  Polar filter the time tendency of tracers 
!  </DATA> 
!  <DATA NAME="rjfrst" TYPE="real">
! Southern latitude below which apply no filtering  
!  </DATA> 
!  <DATA NAME="filter_reflat_s" TYPE="real">
! Southern latitude to which we reference filtering  
!  </DATA> 
!  <DATA NAME="filter_reflat_n" TYPE="real">
! Northern latitude to which we reference filtering    
!  </DATA> 
!</NAMELIST>
!
  use axis_utils_mod,   only: nearest_index 
  use constants_mod,    only: radius, radian, epsln, kelvin
  use diag_manager_mod, only: send_data, register_diag_field
  use fms_mod,          only: open_namelist_file, check_nml_error, close_file
  use fms_mod,          only: write_version_number, FATAL, NOTE 
  use mpp_domains_mod,  only: mpp_update_domains
  use mpp_domains_mod,  only: mpp_global_field
  use mpp_mod,          only: mpp_error, stdout, stdlog

  use ocean_domains_mod,    only: get_local_indices
  use ocean_parameters_mod, only: missing_value
  use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
  use ocean_types_mod,      only: ocean_time_type, ocean_prog_tracer_type

  implicit none

  private

  real, dimension(:,:), allocatable :: yt0, yu0

  integer, dimension(:), allocatable :: numflt    ! number of fir filter applications 
  integer, parameter                 :: jmtfil=23 ! max number of latitudes to be filtered
  integer :: jfrst,  jft1,  jft2,  jfu1,  jfu2
  integer :: num_prog_tracers

  real, dimension(:), allocatable :: spsin, spcos 
  real :: rjfrst=-81.0
  real :: filter_reflat_s=-70.0
  real :: filter_reflat_n=70.0
  real :: refcosn, refcoss
  real :: dtime, dtimer

! for diagnostics 
  integer, allocatable, dimension(:) :: id_filter
  logical :: used

  type(ocean_grid_type), pointer   :: Grd => NULL()
  type(ocean_domain_type), pointer :: Dom => NULL()
  integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

  character(len=128) :: version=&
       '$Id: ocean_polar_filter.F90,v 16.0.108.1 2009/10/10 00:42:45 nnz Exp $'
  character (len=128) :: tagname = &
       '$Name: mom4p1_pubrel_dec2009_nnz $'

  public ocean_polar_filter_init
  public polar_filter_tracers
  private fast_fir
  private set_polar_filtering_indices

  logical :: use_this_module       = .false.
  logical :: polar_filter_field    = .true.
  logical :: polar_filter_tendency = .false.
  logical :: module_is_initialized = .FALSE.

  namelist /ocean_polar_filter_nml/ use_this_module, polar_filter_tendency, polar_filter_field, &
                                    rjfrst, filter_reflat_s, filter_reflat_n

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_polar_filter_init">
!
! <OVERVIEW>
!  Initialize polar filtering module. 
! </OVERVIEW>
!
  subroutine ocean_polar_filter_init(Grid, Domain, Time, T_prog, d_time)

    type(ocean_grid_type),        intent(in), target  :: Grid
    type(ocean_domain_type),      intent(in), target  :: Domain
    type(ocean_time_type),        intent(in)          :: Time
    type(ocean_prog_tracer_type), intent(in)          :: T_prog(:)
    real,                         intent(in)          :: d_time

    real    :: fx, fxa, fxb
    real    :: phi, phit, arg
    real    :: refcos 
    integer :: ioun, ierr, io_status
    integer :: jskpt, jskpu, njtbft, njtbfu
    integer :: numfmx
    integer :: i, j, n

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

    if ( module_is_initialized ) then 
        call mpp_error(FATAL, &
        '==>Error from ocean_polar_filter_mod (ocean_polar_filter_init): module initialized')
    endif

    module_is_initialized = .TRUE.

    call write_version_number( version, tagname )

    Grd => Grid
    Dom => Domain

    call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grd%nk

    num_prog_tracers = size(T_prog(:))

    allocate( id_filter(num_prog_tracers) )
    id_filter(:) = -1

    dtime = d_time
    dtimer = 1.0/dtime

    ! provide for namelist over-ride of default values
    ioun =  open_namelist_file()
    read (ioun,ocean_polar_filter_nml,IOSTAT=io_status)
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_polar_filter_nml)  
    write (stdlogunit,ocean_polar_filter_nml)
    ierr = check_nml_error(io_status,'ocean_polar_filter_nml')
    call close_file (ioun)

    if(use_this_module) then 
        call mpp_error(NOTE, '==>Note from ocean_polar_filter_mod: USING ocean_polar_filter_mod.')
        write (stdoutunit,'(a)') '==>Tracer polar filtering DOES NOT conserve total tracer.'
        write (stdoutunit,'(a)') '   This scheme NOT supported by GFDL scientists. Use at your own risk. '
        allocate (spsin(1-Dom%xhalo:Grd%ni+Dom%xhalo))
        allocate (spcos(1-Dom%xhalo:Grd%ni+Dom%xhalo))
        allocate (numflt(Grd%nj))
    else 
        call mpp_error(NOTE, &
        '==>NOTE from ocean_polar_filter_mod: NOT using ocean_polar_filter_mod.')
        return 
    endif

    if(isc .ne. 1 .or. iec .ne. Grd%ni) then 
       call mpp_error(FATAL,&
       '==>Error from ocean_polar_filter_mod: only implemented for 1D domain decomposition with ni_local=ni_global.')
    endif
    if(Grd%tripolar) then
        call mpp_error(FATAL,&
             '==>Error from ocean_polar_filter_mod: module only useful for spherical grids. ')
    endif

    write (stdoutunit,'(a)') 'Polar filtering details for latitudes'

!   compute sin and cos values for vector correction before filtering

    fx =  1.0e-10
    fxa = (Grd%xt(isc+1,jsc)-Grd%xt(isc,jsc))/radius

    do i=1,Grd%ni
       fxb = fxa*float(i-1)
       spsin(i) = sin(fxb)
       spcos(i) = cos(fxb)
       if (abs(spsin(i)) .lt. fx) spsin(i) = 0.0
       if (abs(spcos(i)) .lt. fx) spcos(i) = 0.0
    enddo

    do i=1,Dom%xhalo
       spsin(Grd%ni+i)=spsin(i)
       spcos(Grd%ni+i)= spcos(i)
       spsin(1-i)=spsin(Grd%ni+1-i)
       spcos(1-i)=spcos(Grd%ni+1-i)
    enddo


!   set up model indices for filtering polar latitudes

    call set_polar_filtering_indices

    jskpt  = jft2-jft1
    jskpu  = jfu2-jfu1
    njtbft = (jft1-jfrst+1)+(Grd%nj-jft2+1)
    njtbfu = (jfu1-jfrst+1)+(Grd%nj-jfu2+1)

    if (njtbft .gt. jmtfil .or. njtbfu .gt. jmtfil) then
        write (stdoutunit,9599) njtbft, njtbfu
        call mpp_error(FATAL,'==>Error from ocean_polar_filter_mod: jmtfil must be >= max(njtbft,njtbfu')
    endif
    if (jfrst .gt. jft1) then
        call mpp_error(FATAL,'==>Error from ocean_polar_filter_mod: jfrst > jft1 is not allowed')
    endif

    write(stdoutunit,9501) filter_reflat_s, filter_reflat_n

    allocate (yu0(Grd%ni, Grd%nj))
    call mpp_global_field(Dom%domain2d, Grd%yu, yu0)

    numfmx = Grd%ni*0.25

    do j=1,Grd%nj
       if ((j .le. jft1 .or. j .ge. jft2) .and. j .ge. jfrst) then
           phi  = yu0(1,j)/radian
           phit = yt0(1,j)/radian
           if (phi .gt. 0.0) then
               refcos = refcosn
           else
               refcos = refcoss
           endif
           if (abs(yt0(1,j)) .ge. 90.0) then
               numflt(j) = 0
           else
               arg = refcos/cos(phit)
               numflt(j) = min(max(1,int(arg)),numfmx)
           endif
           write(stdoutunit,9502) j, numflt(j), yt0(1,j)
       else
           numflt(j) = 0
       endif
    enddo

    do n=1,num_prog_tracers
       id_filter(n) = -1
       if (T_prog(n)%name == 'temp') then
           id_filter(n) = register_diag_field ('ocean_model', 'filter_heat',             &
                Grid%tracer_axes(1:3), Time%model_time, 'heat flux due to polar filter', &
                'Watts/m^2', missing_value=missing_value, range=(/-1.e10,1.e10/))
       else
           id_filter(n) = register_diag_field ('ocean_model', 'filter_'//trim(T_prog(n)%name), &
                Grid%tracer_axes(1:3), Time%model_time, 'rho_dzt * polar filter tendency',     & 
                'kg/(m^2*sec)', missing_value=missing_value, range=(/-1.e10,1.e10/))
       endif
    enddo

9501 format(/' southern reference latitude =',f7.3,' deg.  ' &
         , ' northern reference latitude =',f7.3,' deg.'     &
         , /' T-cell j  numflt  latitude')
9502 format(2i7,6x,f7.3)

9599 format (/,' error => jmtfil must be >= max(njtbft,njtbfu)', &
         /,'          njtbft=',i8,' njtbfu=',i8)


  end subroutine ocean_polar_filter_init
! </SUBROUTINE> NAME="ocean_polar_filter_init"


!#######################################################################
! <SUBROUTINE NAME="polar_filter_tracers">
!
! <OVERVIEW>
!  set up input needed for symmetric finite impulse response 
!  filtering at the specified high latitude row.
! </OVERVIEW>
!
    subroutine polar_filter_tracers(Time, Thickness, T_prog, itemp)

    type(ocean_time_type),        intent(in)    :: Time
    type(ocean_thickness_type),   intent(in)    :: Thickness
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
    integer,                      intent(in)    :: itemp

    integer :: i, j, n, k
    integer :: taum1, tau, taup1
    integer :: numfilter
    real, dimension(isd:ied,nk) :: field_ik  

    if(.not. use_this_module) return 

    if (.not. module_is_initialized ) then 
      call mpp_error(FATAL, &
       '==>Error from ocean_polar_filter_mod (polar_filter_tracers): module not initialized')
    endif 

    taum1 = Time%taum1
    tau   = Time%tau
    taup1 = Time%taup1

    do n = 1, num_prog_tracers
       call mpp_update_domains (T_prog(n)%field(:,:,:,taup1), Dom%domain2d)
    enddo

    ! save pre-filtered field for diagnostics 
    do n=1,num_prog_tracers
       T_prog(n)%wrk1(:,:,:) = T_prog(n)%field(:,:,:,taup1)    
    enddo

    ! filter the tracer field 
    if(polar_filter_field) then 
        do j=jsc,jec
           if (numflt(j+Dom%joff) .ne. 0) then
              numfilter = numflt(j+Dom%joff)
              do n = 1, num_prog_tracers
                  field_ik(:,:) = T_prog(n)%field(:,j,:,taup1)
                  if (n == itemp) field_ik(:,:) = (field_ik(:,:) + kelvin) * Grd%tmask(:,j,:)
                  call fast_fir (field_ik, Grd%tmask(:,j,:), numfilter)
                  if (n == itemp) field_ik(:,:) = (field_ik(:,:) - kelvin) * Grd%tmask(:,j,:)
                  T_prog(n)%field(:,j,:,taup1) = field_ik(:,:)
               enddo
           endif
        enddo
    endif

    ! filter the tracer tendency 
    if(polar_filter_tendency) then 
        do j=jsc,jec
           if (numflt(j+Dom%joff) .ne. 0) then
               numfilter = numflt(j+Dom%joff)
               do n = 1, num_prog_tracers
                  field_ik(:,:) = dtimer*(T_prog(n)%field(:,j,:,taup1)-T_prog(n)%field(:,j,:,taum1))
                  call fast_fir (field_ik, Grd%tmask(:,j,:), numfilter)
                  T_prog(n)%field(:,j,:,taup1) = T_prog(n)%field(:,j,:,taum1) + dtime*field_ik(:,:)
               enddo
           endif
        enddo
    endif

    ! compute tendency due to polar filtering to be used in diagnostics 
    do n = 1, num_prog_tracers
       if (id_filter(n) > 0) then 
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    T_prog(n)%wrk1(i,j,k) =  dtimer*(T_prog(n)%field(i,j,k,taup1) - T_prog(n)%wrk1(i,j,k))
                 enddo
              enddo
               T_prog(n)%wrk1(isc:iec,jsc:jec,k) = T_prog(n)%wrk1(isc:iec,jsc:jec,k)*&
                   T_prog(n)%conversion*Thickness%dzt(isc:iec,jsc:jec,k)
           enddo
           used = send_data  (id_filter(n), T_prog(n)%wrk1(:,:,:), &
                  Time%model_time, rmask=Grd%tmask(:,:,:), &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
       endif
    enddo

  end subroutine polar_filter_tracers
! </SUBROUTINE> NAME="polar_filter_tracers"


!#######################################################################
! <SUBROUTINE NAME="fast_fir">
!
! <OVERVIEW>
! Finite impulse response filter with [.25, .5, .25] weights
! using built in symmetric boundary conditions at land
! </OVERVIEW>
!
! <DESCRIPTION>
!
!     input:
!
!             f     = functions to be filtered
!
!             rmask = mask. must be (1.0,0.0) on (ocean,land) points
!
!             num   = number of filter applications
!
!     output:
!
!             f    = filtered quantities
!
!     author:  r.c.pacanowski   e-mail  Ronald.Pacanowski@noaa.gov
!
! </DESCRIPTION>
!
  subroutine fast_fir (f, rmask, num)

    integer, intent(in)                        :: num
    real, intent(in), dimension(isd:ied,nk)    :: rmask    
    real, intent(inout), dimension(isd:ied,nk) :: f    

    real, dimension(isd:ied,nk) :: s, cw, cc, ce    
    integer                     :: k, i, npass

    if (.not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error from ocean_polar_filter_mod (fast_fir): module not initialized')
    endif 

    if (num .ne. 0) then

!     build weighting functions
        do k=1,Grd%nk
           do i=isc,iec
              cw(i,k) = 0.25*rmask(i-1,k)
              ce(i,k) = 0.25*rmask(i+1,k)
              cc(i,k) = 0.5*(1.0+rmask(i,k)) - cw(i,k) - ce(i,k)
           enddo
        enddo

!     filter "f"
        do npass=1,num
           do k=1,Grd%nk
              do i=isc,iec
                 s(i,k) = cw(i,k)*f(i-1,k) + cc(i,k)*f(i,k)+ce(i,k)*f(i+1,k)
              enddo
              s(isd,k) = s(iec,k)
              s(ied,k) = s(isc,k)
           enddo
           do k=1,nk
              do i=isd,ied
                 f(i,k) = s(i,k)*rmask(i,k)
              enddo
           enddo
        enddo

    endif

  end subroutine fast_fir
! </SUBROUTINE> NAME="fast_fir"


!#######################################################################
! <SUBROUTINE NAME="set_polar_filtering_indices">
!
! <OVERVIEW>
! Set up model indices for filtering polar latitudes
! </OVERVIEW>
!
    subroutine set_polar_filtering_indices

    allocate (yt0(Grd%ni,Grd%nj))
    call mpp_global_field(Dom%domain2d, Grd%yt, yt0)

    if (.not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error ocean_polar_filter_mod (set_polar_filtering_indices): module not initialized')
    endif 

    jfrst  = nearest_index (rjfrst, yt0(1,:))
    refcosn = cos(filter_reflat_n/radian)
    refcoss = cos(filter_reflat_s/radian)
    jft1   = nearest_index (filter_reflat_s, yt0(1,:))
    if (yt0(1,jft1) .gt. filter_reflat_s) jft1 = jft1 - 1
    jft2   = nearest_index (filter_reflat_n, yt0(1,:))
    if (yt0(1,jft2) .lt. filter_reflat_n) jft2 = jft2 + 1
    jfu1   = jft1 - 1
    jfu2   = jft2

  end subroutine set_polar_filtering_indices
! </SUBROUTINE> NAME="set_polar_filtering_indices"

end module ocean_polar_filter_mod
