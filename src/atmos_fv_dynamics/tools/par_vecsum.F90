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
 subroutine par_vecsum(jm, jfirst, jlast, InVector, te0)

! !ROUTINE: par_vecsum --- Calculate vector sum bit-wise consistently

#ifdef SPMD
      use mod_comm, only : mp_sum1d
#ifdef use_shared_pointers
      use fv_arrays_mod, only: fv_thread_bcast
#endif
#endif
      implicit none
! Input:
      integer jm                   ! global latitudes
      integer jfirst               ! first latitude on this PE
      integer jlast                ! last latitude on this PE
      real, intent(in):: InVector(jfirst:jlast)  ! input vector to be summed

! OUTPUT:
      real, intent(out):: te0   ! sum of all vector entries
! Local
      integer j

#ifdef SPMD 
      call mp_sum1d(jm, jfirst, jlast, InVector, te0)
#ifdef use_shared_pointers
      call fv_thread_bcast(te0)
#endif
#else
      te0 = 0.0
      do j=1,jm
        te0 = te0 + InVector(j)
      enddo
#endif

 end subroutine par_vecsum
