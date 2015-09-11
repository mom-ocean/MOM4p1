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
!$Id: r_ratio.F90,v 1.1.2.1 2007/09/11 13:01:20 smg Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update time scale ratio
!
! !INTERFACE:
   subroutine r_ratio(nlev)
!
! !DESCRIPTION:
! This routine updates the ratio $r$ of the dissipation
! time scales as defined in \eq{DefR}.
!
! !USES:
  use turbulence,  only:     tke,eps,kb,epsb
  use turbulence,  only:     r

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer, intent(in)        :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: r_ratio.F90,v $
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
!BOC

!KBK   r = kb/epsb*eps/tke
   r = kb*eps/(epsb*tke)

   return
end subroutine r_ratio

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
