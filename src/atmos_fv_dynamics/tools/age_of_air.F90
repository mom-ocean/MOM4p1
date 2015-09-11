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
      subroutine age_of_air(im, jm, km, jfirst, jlast, ng, time, pfull, q)

      implicit none
      integer im
      integer jm
      integer km
      integer jfirst
      integer jlast  
      integer ng

! q is the age tracer
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.

      real, intent(in):: pfull(im,jfirst:jlast,km)   ! model full level
      real, intent(in):: time        ! accumulated time since init
      real, intent(inout):: q(im,jfirst-ng:jlast+ng,km)

! Local
      integer i, j, k
      real p_source      ! source level (pa)
      real ascale
      real tiny
      parameter ( tiny = 1.e-6 )
      parameter ( p_source = 75000. )
      parameter ( ascale = 5.e-6 / 60. )

!$omp parallel do private(i, j, k)
      do k=1,km
        do j=jfirst, jlast
          do i=1,im
            if( time < tiny ) then
                q(i,j,k) = 0.
            elseif( pfull(i,j,k) >= p_source ) then
                q(i,j,k) = ascale * time
            endif
          enddo
        enddo          ! j-loop
      enddo             ! k-loop

      return
      end
