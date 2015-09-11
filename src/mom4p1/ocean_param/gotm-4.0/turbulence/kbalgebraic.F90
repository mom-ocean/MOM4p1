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
!$Id: kbalgebraic.F90,v 1.1.2.1 2007/09/11 13:01:20 smg Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic kb-equation\label{sec:kbalgebraic}
!
! !INTERFACE:
   subroutine kbalgebraic(nlev)
!
! !DESCRIPTION:
! The algebraic equation for $k_b$ simply assumes equilibrium in \eq{kbeq},
! \begin{equation}
!   \label{kbEquilibrium}
!   P_b = \epsilon_b
!   \point
! \end{equation}
! This equation can be re-written as
! \begin{equation}
!   \label{kbAgebraic}
!   k_b = \dfrac{k_b \epsilon}{k \epsilon_b} \dfrac{k}{\epsilon} P_b
!       = r \dfrac{k}{\epsilon} P_b = c_b \dfrac{k}{\epsilon} P_b
!   \comma
! \end{equation}
! where we used the definition of the time scale ratio $r$ in
! \eq{DefR}, and assumed that $r=c_b$ is a constant.
!

!
! !USES:
   use turbulence,  only:     tke,eps,kb,Pb
   use turbulence,  only:     ctt,kb_min

  IMPLICIT NONE
!
! !INPUT PARAMETERS:

! number of vertical layers
   integer,  intent(in)                 :: nlev

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: kbalgebraic.F90,v $
!  Revision 1.1.2.1  2007/09/11 13:01:20  smg
!  Add these files to the mom4p1 branch.
!  AUTHOR:Griffies
!  REVIEWERS:
!  TEST STATUS:
!  CHANGES PUBLIC INTERFACES?
!  CHANGES ANSWERS?
!
!  Revision 1.1  2005-06-27 10:54:33  kbk
!  new files needed
!
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:

   integer                             :: i

!-----------------------------------------------------------------------
!BOC

   do i=0,nlev
      kb(i) = ctt*tke(i)/eps(i)*Pb(i)

      !  clip at kb_min
      kb(i) = max(kb(i),kb_min)
   enddo

   return
   end subroutine kbalgebraic

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
