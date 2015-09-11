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

module fs_profile_mod

use fms_mod, only : mpp_pe, mpp_root_pe, write_version_number, &
                    error_mesg, FATAL

implicit none
private

!-----------------------------------------------------------------------
! **      THIS PROGRAM CALCULATES TEMPERATURES ,H2O MIXING RATIOS     **
! **      AND O3 MIXING RATIOS BY USING AN ANALYTICAL                 **
! **      FUNCTION WHICH APPROXIMATES                                 **
! **      THE US STANDARD (1976).  THIS IS                            **
! **      CALCULATED IN FUNCTION 'ANTEMP', WHICH IS CALLED BY THE     **
! **      MAIN PROGRAM.  THE FORM OF THE ANALYTICAL FUNCTION WAS      **
! **      SUGGESTED TO ME IN 1971 BY RICHARD S. LINDZEN.              **
!
!*****THIS VERSION IS ONLY USABLE FOR 1976 US STD ATM AND OBTAINS
!     QUANTITIES FOR CO2 INTERPOLATION AND INSERTION INTO OPERA-
!     TIONAL RADIATION CODES
!
!    definitions:       
!    -----------
!      pd,pd8: pressures (mb) for data levels. pd is for the case where
!              p(sfc)=1013.25 mb; pd8 applies when p(sfc)=810.6 mb.
!              in either case, index (nlev+1) is at the sfc.
!      press:  same as pd, but with indices reversed,index 1 at the
!              surface, and index (nlev+1) at the top (nonzero) data
!              level.
!-----------------------------------------------------------------------

 Public   fs_profile, fs_profile_init, fs_profile_end

 Real, Public, Allocatable, Dimension(:) :: pd1013,plm1013,pd810,plm810

!------------ VERSION NUMBER ----------------

 character(len=128) :: version = '$Id: fs_profile.F90,v 10.0 2003/10/24 22:00:32 fms Exp $'
 character(len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'
 logical            :: module_is_initialized = .false.

CONTAINS

!#######################################################################

  subroutine fs_profile (Pref,stemp,gtemp)

!-----------------------------------------------------------------------
  Real, Intent(IN) , Dimension(:,:) :: Pref
  Real, Intent(OUT), Dimension(:)   :: stemp,gtemp
!-----------------------------------------------------------------------

      call error_mesg('fs_profile', &
      'This module is not supported as part of the public release', FATAL)

  end subroutine fs_profile

!#######################################################################

      Subroutine fs_profile_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('fs_profile_init', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine fs_profile_init

!#######################################################################

      Subroutine fs_profile_end

      module_is_initialized = .false.
!---------------------------------------------------------------------

      call error_mesg('fs_profile_end', &
      'This module is not supported as part of the public release', FATAL)

      End Subroutine fs_profile_end

!#######################################################################


end module fs_profile_mod

