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
module ocean_workspace_mod
!
! <CONTACT EMAIL="Matthew.Harrison@noaa.gov">
! M.J. Harrison 
! </CONTACT>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov">
! S.M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! This module allocates some workspace for use in mom4. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module allocates some workspace for use in mom4. 
!</DESCRIPTION>
!
use ocean_types_mod, only: ocean_domain_type, ocean_grid_type

implicit none

private

#include <ocean_memory.h>

#ifdef MOM4_STATIC_ARRAYS
real, public, dimension(isd:ied,jsd:jed)      :: wrk1_2d
real, public, dimension(isd:ied,jsd:jed)      :: wrk2_2d
real, public, dimension(isd:ied,jsd:jed)      :: wrk3_2d
real, public, dimension(isd:ied,jsd:jed)      :: wrk4_2d

real, public, dimension(isd:ied,jsd:jed,2)    :: wrk1_v2d
real, public, dimension(isd:ied,jsd:jed,2)    :: wrk2_v2d

real, public, dimension(isd:ied,jsd:jed,nk)   :: wrk1
real, public, dimension(isd:ied,jsd:jed,nk)   :: wrk2
real, public, dimension(isd:ied,jsd:jed,nk)   :: wrk3
real, public, dimension(isd:ied,jsd:jed,nk)   :: wrk4

real, public, dimension(isd:ied,jsd:jed,0:nk) :: wrk1_zw
real, public, dimension(isd:ied,jsd:jed,0:nk) :: wrk2_zw

real, public, dimension(isd:ied,jsd:jed,nk,2) :: wrk1_v
real, public, dimension(isd:ied,jsd:jed,nk,2) :: wrk2_v
real, public, dimension(isd:ied,jsd:jed,nk,2) :: wrk3_v
real, public, dimension(isd:ied,jsd:jed,nk,2) :: wrk4_v

#else

real, public, allocatable, dimension(:,:)     :: wrk1_2d
real, public, allocatable, dimension(:,:)     :: wrk2_2d
real, public, allocatable, dimension(:,:)     :: wrk3_2d
real, public, allocatable, dimension(:,:)     :: wrk4_2d

real, public, allocatable, dimension(:,:,:)   :: wrk1
real, public, allocatable, dimension(:,:,:)   :: wrk2
real, public, allocatable, dimension(:,:,:)   :: wrk3
real, public, allocatable, dimension(:,:,:)   :: wrk4

real, public, allocatable, dimension(:,:,:)   :: wrk1_zw
real, public, allocatable, dimension(:,:,:)   :: wrk2_zw

real, public, allocatable, dimension(:,:,:,:) :: wrk1_v
real, public, allocatable, dimension(:,:,:,:) :: wrk2_v
real, public, allocatable, dimension(:,:,:,:) :: wrk3_v
real, public, allocatable, dimension(:,:,:,:) :: wrk4_v

real, public, allocatable, dimension(:,:,:)   :: wrk1_v2d
real, public, allocatable, dimension(:,:,:)   :: wrk2_v2d

#endif

public ocean_workspace_init, ocean_workspace_end



logical, private :: module_is_initialized = .FALSE.

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_workspace_init">
!
! <DESCRIPTION>
! Initialize mom4 workspace module.
! </DESCRIPTION>
!
subroutine ocean_workspace_init(Domain, Grid)

type(ocean_domain_type), intent(in) :: Domain
type(ocean_grid_type), intent(in)   :: Grid

module_is_initialized = .TRUE.

#ifndef MOM4_STATIC_ARRAYS
allocate(wrk1_2d(Domain%isd:Domain%ied,Domain%jsd:Domain%jed))
allocate(wrk2_2d(Domain%isd:Domain%ied,Domain%jsd:Domain%jed))
allocate(wrk3_2d(Domain%isd:Domain%ied,Domain%jsd:Domain%jed))
allocate(wrk4_2d(Domain%isd:Domain%ied,Domain%jsd:Domain%jed))

allocate(wrk1(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk))
allocate(wrk2(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk))
allocate(wrk3(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk))
allocate(wrk4(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk))

allocate(wrk1_zw(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,0:Grid%nk))
allocate(wrk2_zw(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,0:Grid%nk))

allocate(wrk1_v(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk,2))
allocate(wrk2_v(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk,2))
allocate(wrk3_v(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk,2))
allocate(wrk4_v(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,Grid%nk,2))

allocate(wrk1_v2d(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,2))
allocate(wrk2_v2d(Domain%isd:Domain%ied,Domain%jsd:Domain%jed,2))
#endif

end subroutine ocean_workspace_init
! </SUBROUTINE> NAME="ocean_workspace_init">


!#######################################################################
! <SUBROUTINE NAME="ocean_workspace_end">
!
! <DESCRIPTION>
! End mom4 workspace.
! </DESCRIPTION>
!
subroutine ocean_workspace_end()

module_is_initialized = .FALSE.

#ifndef MOM4_STATIC_ARRAYS
deallocate(wrk1,wrk2,wrk3,wrk4,wrk1_zw,wrk2_zw,wrk1_v,wrk2_v,wrk3_v,wrk4_v)
deallocate(wrk1_2d,wrk2_2d,wrk3_2d,wrk4_2d,wrk1_v2d,wrk2_v2d)
#endif

end subroutine ocean_workspace_end
! </SUBROUTINE> NAME="ocean_workspace_end">

end module ocean_workspace_mod

