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
      real function gmean(im, jm, jfirst, jlast, q)
      use fv_pack, only: grid_weight
      implicit none

! !INPUT PARAMETERS:
      integer  im, jm                        ! Horizontal dimensions
      integer  jfirst, jlast                 ! Latitude strip
      real q(im,jfirst:jlast)            ! 2D field 

      integer i, j
      real xsum(jfirst:jlast)

          xsum = 0.
      do j=jfirst,jlast
        do i=1,im
          xsum(j) = xsum(j) + q(i,j)
        enddo
          xsum(j) = xsum(j)*grid_weight(j)
      enddo

      call par_vecsum( jm, jfirst, jlast, xsum, gmean)
      end
