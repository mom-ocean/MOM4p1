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
      subroutine upper(string,length)

!***********************************************************************
!
!     upper.f - change lower case letter to upper case letter          *
!                                                                      *
!     George Lai Tue Jun 28 16:37:00 1994                              *
!                                                                      *
!***********************************************************************

      implicit         none

!      character string(length)
!      character*20 string
!      character, dimension(length) :: string
!      character (len=*), intent(inout) ::  string
!      character (len=*) ::  string
!      character (len=1), intent(inout) ::  string(20)
!ok      character (len=20), intent(inout) ::  string
      character (len=*), intent(inout) ::  string
      character char1
      integer,   intent(in)    ::  length
      integer i
      integer a, z, dist
      a = ichar('a')
      z = ichar('z')
      dist = ichar('A') - a

      do i = 1,length
        char1=string(i:i)
        if (ichar(char1) .ge. a .and.       &
            ichar(char1) .le. z) then
          string(i:i) = char(ichar(char1)+dist)
        endif
      end do

      return
      end

