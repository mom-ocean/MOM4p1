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
!$Id: fk_craig.F90,v 1.1.2.1 2007/09/11 13:01:19 smg Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: TKE flux from wave-breaking\label{sec:fkCraig}
!
! !INTERFACE:
   REALTYPE  function fk_craig(u_tau)
!
! !DESCRIPTION:
! This functions returns the flux of $k$ caused by breaking surface waves
! according to
! \begin{equation}
!  \label{craig}
!   F_k = \eta u_*^3
!  \point
! \end{equation}
! This form has also been used by \cite{CraigBanner94}, who suggested
! $\eta \approx 100$.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: u_tau
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter                 :: eta=100.
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: fk_craig.F90,v $
!  Revision 1.1.2.1  2007/09/11 13:01:19  smg
!  Add these files to the mom4p1 branch.
!  AUTHOR:Griffies
!  REVIEWERS:
!  TEST STATUS:
!  CHANGES PUBLIC INTERFACES?
!  CHANGES ANSWERS?
!
!  Revision 1.4  2005-06-27 13:44:07  kbk
!  modified + removed traling blanks
!
!  Revision 1.3  2004/08/18 12:50:57  lars
!  updated documentation
!
!  Revision 1.2  2003/03/28 09:20:35  kbk
!  added new copyright to files
!
!  Revision 1.1  2003/03/10 09:00:36  gotm
!  Part of new generic turbulence model
!
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fk_craig = eta*u_tau**3.

   end function fk_craig
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
