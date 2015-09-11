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
      module MO_JPL_MOD

implicit none
character(len=128), parameter :: version     = '$Id: mo_jpl.F90,v 13.0 2006/03/28 21:16:17 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: mom4p1_pubrel_dec2009_nnz $'
logical                       :: module_is_initialized = .false.

      CONTAINS

      subroutine JPL( rate, m, factor, ko, kinf,plnplv )
!-----------------------------------------------------------------
!        ... Calculate JPL troe rate
!-----------------------------------------------------------------
      
      implicit none

!-----------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------
      integer,intent(in)::    plnplv
      real, intent(in)  ::   factor
      real, intent(in)  ::   ko(plnplv)
      real, intent(in)  ::   kinf(plnplv)
      real, intent(in)  ::   m(plnplv)
      real, intent(out) ::   rate(plnplv)

!-----------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------
      real  ::     xpo( SIZE(rate) )
      

      xpo(:)  = ko(:) * m(:) / kinf(:)
      rate(:) = ko(:) / (1. + xpo(:))
      xpo(:)  = LOG10( xpo(:) )
      xpo(:)  = 1. / (1. + xpo(:)*xpo(:))
      rate(:) = rate(:) * factor**xpo(:)

      end subroutine JPL

      end module MO_JPL_MOD
