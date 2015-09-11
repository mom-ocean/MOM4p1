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
module ocean_riverspread_mod
!  
!<CONTACT EMAIL="Mike.Spelman@noaa.gov"> Mike Spelman
!</CONTACT>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! Spread runoff or calving horizontally over a region determined 
! by a table. 
!</OVERVIEW>
!
!<DESCRIPTION>
! At some coastal ocean gridpoints, the runoff or calving flux contribution to 
! (p-e+r) may be very large because of local insertion to ocean. Therefore, 
! we may choose to spread the large river runoff over neighboring pairs
! of gridpoints. Annual mean river runoff greater than 0.05 m/day is 
! considered to be very large.
!</DESCRIPTION>
!
! <INFO>
!
! <NOTE>
! Spreading in a 2D lat-long field is implemented in the following manner:
!
! If the spreading region lives within the halo region
! (i.e., within same local_domain), 
! then no added mpp communication required. However, more generally 
! the spreading region can extend beyond the existing halo region.
! In this case, spread_domain
! is defined so that its halos incorporate the maximum separation 
! of spreading points.  New tracer and grid arrays
! are defined over this extended spread_domain.  This added domain
! size will come at some computational cost, so it is advantageous
! to choose the spreading region carefully.
! </NOTE>
!
! <NOTE>
! The current implementation of spreading has not been tested 
! for a spreading region that is separated by either 
! a zonal cyclic condition or across the tripolar fold.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_riverspread_nml">
!
!  <DATA NAME="use_this_module"  TYPE="logical">
!  Must be true to enable this module.  Default=.false.
!  </DATA>
!  <DATA NAME="debug_this_module"  TYPE="logical">
!  For debugging. Default=.false.
!  </DATA>
!
!</NAMELIST>
!
use field_manager_mod,   only: MODEL_OCEAN, find_field_index, get_field_methods
use field_manager_mod,   only: method_type, get_field_info, parse
use fms_mod,             only: write_version_number, open_namelist_file, check_nml_error, close_file
use fms_mod,             only: stdout, stdlog, FATAL, NOTE
use mpp_domains_mod,     only: domain2D, mpp_update_domains, mpp_define_domains
use mpp_domains_mod,     only: cyclic_global_domain, global_data_domain
use mpp_mod,             only: mpp_error

use ocean_domains_mod,   only: get_local_indices
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type, ocean_options_type

implicit none

private 

public ocean_riverspread_init
public spread_river_horz
private on_comp_domain
private on_data_domain

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()

! Variables set inside subroutine ocean_riverspread_init
! See comments in ocean_riverspread_init for description
integer, dimension(:),   allocatable :: ispread, jspread
integer, dimension(:),   allocatable :: is, ie
integer, dimension(:),   allocatable :: js, je
integer, dimension(:),   allocatable :: halon     ! halo size for each spreading region
integer                              :: nspread=0 ! number of spreading regions

! Variables calculated from user input
integer                              :: spread_halo     ! halo size needed to perform spreading
type(ocean_domain_type), save        :: Spread_domain   ! domain for spreading 

! variables defined with halos given by spread_halo
real, dimension(:,:), allocatable    :: field_spreadx  ! 
real, dimension(:,:), allocatable    :: xtx            ! zonal position of tracer points (degrees)
real, dimension(:,:), allocatable    :: ytx            ! meridional position of tracer points (degrees)
real, dimension(:,:), allocatable    :: datx           ! area of tracer cell (m^2)
real, dimension(:,:), allocatable    :: tmaskx         !

integer                              :: unit=6         ! processor zero writes to unit 6

logical, dimension(:), allocatable   :: error_spread   ! for checking that all spread points are OK.  

character(len=256) :: version=&
       '$Id: ocean_riverspread.F90,v 16.0.70.2.18.1 2009/10/10 00:42:47 nnz Exp $)'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

logical :: use_this_module       = .false.
logical :: debug_this_module     = .false.
logical :: module_is_initialized = .FALSE.

namelist /ocean_riverspread_nml/ use_this_module, debug_this_module 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_riverspread_init">
!
! <DESCRIPTION>
!    Initial set up for spreading of tracers 
!
!    (i,j) locations of points to be spread are set in data
!    statements.
!
!    Checks are performed to ensure that the spreading
!    grid locations are valid according to model configuration.
!
!    A summary of the locations of spreading points is written out.
!
!    User specified inputs in "USER INPUT" section:
!
!    ispread and jspread = user specified i,j grid locations of data for spreading.
!
!    is, ie, js, je = user specified i,j grid locations of 
!                     the corners of the spreading regions.
!                     is and ie are the east and west coord of each region,
!                     js and je are the south and north coord of each region. 
!
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_riverspread_nml">
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging 
!  </DATA> 
!</NAMELIST>

subroutine ocean_riverspread_init(Grid, Domain, Ocean_options)

  type(ocean_grid_type),    intent(in), target :: Grid
  type(ocean_domain_type),  intent(in), target :: Domain
  type(ocean_options_type), intent(inout)      :: Ocean_options

  type(method_type), allocatable, dimension(:) :: spread_methods

  character(len=32) :: fld_type, fld_name
  integer           :: i, n, model
  integer           :: io_status, ioun, ierr
  integer           :: isp, jsp
  integer           :: nsp, imax, imin, spread_halo_test
  integer           :: parse_ok

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then  
    call mpp_error(FATAL, &
    '==>Error in ocean_riverspread_mod (ocean_riverspread_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk

  Dom => Domain
  Grd => Grid

  n = find_field_index(MODEL_OCEAN,'riverspread')

  if (n < 1) then 
      write(stdoutunit,*)' '
      write(stdoutunit,*) &
     '==>Warning: ocean_riverspread_init: n<1 for riverspread table. Will NOT use ocean_riverspread.'  
      write(stdoutunit,*)' '
      Ocean_options%riverspread = 'Did NOT use riverspread module.'
      return
  endif

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_riverspread_nml,IOSTAT=io_status)
  write (stdlogunit,ocean_riverspread_nml)
  write (stdoutunit,'(/)') 
  write (stdoutunit,ocean_riverspread_nml)
  ierr = check_nml_error(io_status,'ocean_riverspread_nml')
  call close_file (ioun)

  call get_field_info(n,fld_type,fld_name,model,nspread)

    if(use_this_module) then 
      call mpp_error(NOTE,&
      '==>From ocean_riverspread_mod: Using riverspread module to mix river water into the ocean.')
      Ocean_options%riverspread = 'Used riverspread module to horizontally spread river water away from land.'
    else
      call mpp_error(NOTE,&
      '==>From ocean_riverspread_mod: NOT using riverspread module, yet there are n > 0 riverspread points.')
      Ocean_options%riverspread = 'Did NOT use riverspread module.'
      return
    endif

  if(nspread == 0) then 
     call mpp_error(NOTE,&
     '==>Warning in ocean_riverspread_mod: No river spreading points chosen.')
     Ocean_options%riverspread = 'Did NOT use riverspread module.'
     return
  endif 

  if(Grd%tripolar) then 
    call mpp_error(NOTE,&
    '==>Warning in ocean_riverspread_mod: river spreading not implemented for points across bipolar fold.')
  endif

  allocate(spread_methods(nspread))

  allocate (ispread(nspread))
  ispread(:) = 0
  allocate (jspread(nspread))
  jspread(:) = 0
  allocate (is(nspread), ie(nspread))
  is(:) = 0; ie(:) = 0
  allocate (js(nspread), je(nspread))
  js(:) = 0; je = 0
  allocate (error_spread(nspread))
  error_spread(:) = .false.
  allocate (halon(nspread))
  halon(:) = 0

  call get_field_methods(n,spread_methods)
  do i=1, nspread
     parse_ok = parse(spread_methods(i)%method_control,'iloc',ispread(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error in ocean_riverspread_mod: spread table entry "iloc" error')
     parse_ok = parse(spread_methods(i)%method_control,'jloc',jspread(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error in ocean_riverspread_mod: spread table entry "jloc" error')
     parse_ok = parse(spread_methods(i)%method_control,'is',is(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error in ocean_riverspread_mod: spread table entry "is" error')
     parse_ok = parse(spread_methods(i)%method_control,'ie',ie(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error in ocean_riverspread_mod: spread table entry "ie" error')
     parse_ok = parse(spread_methods(i)%method_control,'js',js(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error in ocean_riverspread_mod: spread table entry "js" error')
     parse_ok = parse(spread_methods(i)%method_control,'je',je(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error in ocean_riverspread_mod: spread table entry "je" error')
  enddo

  spread_halo = 1
  do nsp=1,nspread

    ! check for invalid spreading grid locations
    if (ispread(nsp) < 1 .or. ispread(nsp) > Grd%ni) then
        write(unit,991) nsp, ispread(nsp)
        error_spread(nsp) = .true.
     endif
     if (jspread(nsp) < 1 .or. jspread(nsp) > Grd%nj ) then
        write(unit,992) nsp, jspread(nsp)
        error_spread(nsp) = .true.
     endif

   ! determine size for spread_halo
   if(Grd%cyclic_x) then 
     imax = max(is(nsp),ie(nsp))
     imin = min(is(nsp),ie(nsp))
     spread_halo_test = min( abs(imax-imin), abs(imax-Grd%ni-imin)) 
   else 
     spread_halo_test = abs(is(nsp)-ie(nsp))
   endif 
   if (spread_halo_test > spread_halo ) spread_halo = spread_halo_test

   spread_halo_test = abs(js(nsp)-je(nsp))
   if (spread_halo_test > spread_halo ) spread_halo = spread_halo_test
   halon(nsp) = spread_halo_test

  enddo  ! end of do-loop nsp=1,nspread

  if(spread_halo > Dom%xhalo .or. spread_halo > Dom%yhalo ) then 
    write(stdoutunit,*)'Defining extra tracer arrays for "riverspread_domain". Halo for spreading = ', spread_halo
    write(stdoutunit,*)'Halos for the spreading regions = ', halon
  endif 

  ! define arrays and spread_domain
  if(spread_halo <=Dom%xhalo .and. spread_halo <= Dom%yhalo) then 

     write(stdoutunit,*)'Local computational domain is large enough for chosen points in river spreading,'

     allocate (datx(isd:ied,jsd:jed))
     allocate (xtx(isd:ied,jsd:jed))
     allocate (ytx(isd:ied,jsd:jed))
     allocate (tmaskx(isd:ied,jsd:jed))
     xtx(:,:)    = Grd%xt(:,:)
     ytx(:,:)    = Grd%yt(:,:)
     datx(:,:)   = Grd%dat(:,:)
     tmaskx(:,:) = Grd%tmask(:,:,1)

     allocate (field_spreadx(isd:ied,jsd:jed))
     field_spreadx(:,:) = 0.0
  endif 

  if(spread_halo > Dom%xhalo .or. spread_halo > Dom%yhalo) then 
     write(stdoutunit,*)'The model local computational domain has a halo = ',Dom%xhalo, Dom%yhalo
     write(stdoutunit,*)'This is smaller than the halo required for river spreading.'
     write(stdoutunit,*)'For spreading, will define a new domain type with halo = ', spread_halo

     call mpp_define_domains((/1,Grid%ni,1,Grid%nj/),Domain%layout,Spread_domain%domain2d, maskmap=Domain%maskmap,&
           xflags = CYCLIC_GLOBAL_DOMAIN, xhalo = spread_halo, yhalo = spread_halo, name='spread',&
           x_cyclic_offset = Domain%x_cyclic_offset, y_cyclic_offset = Domain%y_cyclic_offset)

     ! time independent arrays defined over the spread_domain region 
     allocate (datx(isc-spread_halo:iec+spread_halo,jsc-spread_halo:jec+spread_halo))
     allocate (xtx(isc-spread_halo:iec+spread_halo,jsc-spread_halo:jec+spread_halo))
     allocate (ytx(isc-spread_halo:iec+spread_halo,jsc-spread_halo:jec+spread_halo))
     allocate (tmaskx(isc-spread_halo:iec+spread_halo,jsc-spread_halo:jec+spread_halo))
     datx(:,:)   = 0.0
     xtx(:,:)    = 0.0
     ytx(:,:)    = 0.0
     tmaskx(:,:) = 0.0

     datx(isc:iec,jsc:jec)    = Grd%dat(isc:iec,jsc:jec)
     xtx(isc:iec,jsc:jec)     = Grd%xt(isc:iec,jsc:jec)
     ytx(isc:iec,jsc:jec)     = Grd%yt(isc:iec,jsc:jec)
     tmaskx(isc:iec,jsc:jec)  = grd%tmask(isc:iec,jsc:jec,1)
     call mpp_update_domains (datx(:,:), Spread_domain%domain2d)
     call mpp_update_domains (xtx(:,:), Spread_domain%domain2d)
     call mpp_update_domains (ytx(:,:), Spread_domain%domain2d)
     call mpp_update_domains (tmaskx(:,:), Spread_domain%domain2d)

     ! time dependent arrays defined over the spread_domain region 
     allocate (field_spreadx(isc-spread_halo:iec+spread_halo,jsc-spread_halo:jec+spread_halo))
     field_spreadx(:,:) = 0.0

  endif


  do nsp=1,nspread

    isp = ispread(nsp)-Dom%ioff
    jsp = jspread(nsp)-Dom%joff

    ! check for attempts to spread land rather than sea
     if(on_comp_domain(nsp)) then 
        if (Grd%kmt(isp,jsp) == 0) then
           write (unit,192) ispread(nsp),jspread(nsp),Grd%kmt(isp,jsp)
           error_spread(nsp)=.true.
        endif
     endif

     ! write out summary information for this spreading point
     if(on_comp_domain(nsp)) then
       write(unit,191) nsp, ispread(nsp), jspread(nsp) &
          ,xtx(isp,jsp), ytx(isp,jsp) &
          ,is(nsp),ie(nsp),js(nsp),je(nsp)
     endif

  enddo

  ! bring model down if problems with spreading points  
  do nsp=1,nspread
    if(error_spread(nsp)) then
      call mpp_error(FATAL,&
      '==>Error in ocean_riverspread_mod: problem in river spread initialization.')
    endif
  enddo

191 format(/' ===== from ocean_riverspread_init ====='/                            &
       ' for spreading location number',i3,/                                       &
       ' spread (i,j) gridpoint (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/    &
       ' over (is, ie, js, je)= ', 4i6)
192 format(/'=>Error: problem with spreading:',/,                                  &
       ' land point requested'/                                                    &
       '      at ispread(',i4,'), jspread(',i4,'),', ' kmt = ',i4,/)
991 format(/'=>Error: problem with spreading:',/,                                  &
       ' out of bounds i grid location requested'/                                 &
       '      ispread(',i4,') was set incorrectly set to = ',i8,/)
992 format(/'=>Error: problem with spreading:',/,                                  &
       ' out of bounds j-row grid location requested'/                             &
       '      jspread(',i4,') was set incorrectly set to = ',i8,/)


end subroutine ocean_riverspread_init
! </SUBROUTINE> NAME="ocean_riverspread_init"


!#######################################################################
! <SUBROUTINE NAME="spread_river_horz">
!
! <DESCRIPTION>
! Provide conservative spreading of river runoff field. 
! </DESCRIPTION>
!
subroutine spread_river_horz (field_spread)

  real, intent(inout), dimension(isd:,jsd:) :: field_spread

  real    :: spread, wtsum
  integer :: i, j, nsp

  if(nspread == 0)  return

  if(spread_halo > Dom%xhalo .or. spread_halo > Dom%yhalo) then 

    field_spreadx(isc:iec,jsc:jec) = field_spread(isc:iec,jsc:jec)
    call mpp_update_domains (field_spreadx(:,:), Spread_domain%domain2d)

    do nsp=1,nspread 

      if(at_least_one_in_spread_domain(nsp)) then  

           if(debug_this_module) then 
             write(unit,'(a,i3,a,i3,a,8i5)') '=>in spread_river_horz: spreading region= ', nsp, &
                                         ' halo= ', halon(nsp), ' region (is, ie, js, je)= ',   &
                                         is(nsp),ie(nsp),js(nsp),je(nsp),                       &
                                         is(nsp)-Dom%ioff, ie(nsp)-Dom%ioff, js(nsp)-Dom%joff, je(nsp)-Dom%joff
           endif 
 
           spread = 0.0; wtsum  = 0.0
           do j=js(nsp)-Dom%joff,je(nsp)-Dom%joff
             do i=is(nsp)-Dom%ioff,ie(nsp)-Dom%ioff
               spread = spread + field_spreadx(i,j) * datx(i,j) * tmaskx(i,j)
               wtsum = wtsum + datx(i,j) * tmaskx(i,j)
             enddo
           enddo

           if (wtsum > 0.0) then
             spread = spread / wtsum
           else
             write(unit,'(a,2i6)') '=>ERROR: spreading = 0 at i,j=', ispread(nsp), jspread(nsp)
             call mpp_error(FATAL,'==>Error: Problem in spread_river_horz.')
           endif

           do j=js(nsp)-Dom%joff,je(nsp)-Dom%joff
             do i=is(nsp)-Dom%ioff,ie(nsp)-Dom%ioff
               field_spreadx(i,j) = spread * tmaskx(i,j)
             enddo
           enddo
           do j=jsc,jec
             do i=isc,iec
               field_spread(i,j) = field_spreadx(i,j)
             enddo
           enddo

      endif ! end of at_least_one_in_spread_domain(nsp) if-test 

    enddo  ! end of the nspread do-loop 

  endif ! end of spread_halo > halo if-test 

  if(spread_halo <= Dom%xhalo .and. spread_halo <= Dom%yhalo) then 

    do nsp=1,nspread 

      if(at_least_one_in_spread_domain(nsp)) then  
  
           if(debug_this_module) then 
             write(unit,'(a,i3,a,i3,a,8i5)') '=>in spread_river_horz, compute spreading for region= ', nsp, &
                                         ' halo= ', halon(nsp), ' region (is, ie, js, je)= ',               &
                                           is(nsp),ie(nsp),js(nsp),je(nsp),                                 &
                                           is(nsp)-Dom%ioff, ie(nsp)-Dom%ioff, js(nsp)-Dom%joff, je(nsp)-Dom%joff
           endif 

           spread = 0.0; wtsum  = 0.0
           do j=js(nsp)-Dom%joff,je(nsp)-Dom%joff
             do i=is(nsp)-Dom%ioff,ie(nsp)-Dom%ioff
               spread = spread + field_spread(i,j) * Grd%dat(i,j) * Grd%tmask(i,j,1)
               wtsum = wtsum + Grd%dat(i,j) * Grd%tmask(i,j,1)
             enddo
           enddo

           if (wtsum > 0.0) then
             spread = spread / wtsum
           else
             write(unit,'(a,2i6)') '=>ERROR: spreading = 0 at i,j=', ispread(nsp), jspread(nsp)
             call mpp_error(FATAL,'==>Error: Problem in spread_river_horz.')
           endif

           do j=js(nsp)-Dom%joff,je(nsp)-Dom%joff
             do i=is(nsp)-Dom%ioff,ie(nsp)-Dom%ioff
               field_spread(i,j) = spread *Grd%tmask(i,j,1)
             enddo
           enddo

      endif ! end of at_least_one_in_spread_domain(nsp) if-test 

    enddo  ! end of the nspread do-loop 

  endif ! end of spread_halo <= halo if-test 

end subroutine spread_river_horz
! </SUBROUTINE> NAME="spread_river_horz"


!#######################################################################
! <FUNCTION NAME="at_least_one_in_spread_domain">
!
! <DESCRIPTION>
! Function to see if at least one of the points in the spreading region
! is within the computational domain for the processor. 
! </DESCRIPTION>
!
! <IN NAME="nsp" TYPE="integer">
! Integer labeling the particular spreading region
! </IN>
!
function at_least_one_in_spread_domain(nsp)

  integer, intent(in) :: nsp
  logical :: in_domain, at_least_one_in_spread_domain

  if(isc-halon(nsp)+Dom%ioff <= is(nsp) .and. ie(nsp) <= iec+halon(nsp)+Dom%ioff .and. &
     jsc-halon(nsp)+Dom%joff <= js(nsp) .and. je(nsp) <= jec+halon(nsp)+Dom%joff) then
     in_domain = .true.
  else 
     in_domain = .false.  
  endif

  at_least_one_in_spread_domain = .false.
  if(in_domain) at_least_one_in_spread_domain = .true.


end function at_least_one_in_spread_domain 
! </FUNCTION> NAME="at_least_one_in_spread_domain"


!#######################################################################
! <FUNCTION NAME="on_comp_domain">
!
! <DESCRIPTION>
! Determine if the point is in comp-domain for the processor
! </DESCRIPTION>
!
function on_comp_domain(nsp)

  integer, intent(in)  :: nsp
  logical :: on_comp_domain

  if(isc+Dom%ioff <= ispread(nsp) .and. ispread(nsp) <= iec+Dom%ioff .and. &
     jsc+Dom%joff <= jspread(nsp) .and. jspread(nsp) <= jec+Dom%joff) then
     on_comp_domain = .true.
  else 
     on_comp_domain = .false.  
  endif

end function on_comp_domain
! </FUNCTION> NAME="on_comp_domain"


!#######################################################################
! <FUNCTION NAME="on_data_domain">
!
! <DESCRIPTION>
! Determine if the point is in data-domain for the processor
! </DESCRIPTION>
!
function on_data_domain(nsp)

  integer, intent(in)  :: nsp
  logical :: on_data_domain

  if(isd+Dom%ioff <= ispread(nsp) .and. ispread(nsp) <= ied+Dom%ioff .and. &
     jsd+Dom%joff <= jspread(nsp) .and. jspread(nsp) <= jed+Dom%joff) then
     on_data_domain = .true.
  else 
     on_data_domain = .false.  
  endif

end function on_data_domain
! </FUNCTION> NAME="on_data_domain"


end module ocean_riverspread_mod
