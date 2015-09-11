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
!$Id: variances.F90,v 1.1.2.1 2007/09/11 13:01:21 smg Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: The algebraic velocity variances \label{sec:variances}
!
! !INTERFACE:
   subroutine variances(nlev,SSU,SSV)
!
! !DESCRIPTION:

!  Using \eq{bijVertical} and the solution shown in \eq{b13} and
!  the variances of the turbulent velocity fluctations can be
!  evaluated according to
!  \begin{equation}
!    \label{variances}
!    \begin{array}{rcl}
!      \dfrac{\langle u'^2 \rangle}{k}
!      &=& \dfrac23 + \dfrac{1}{\mathcal{N}\eps} \left(
!        \left(\dfrac{a_2}{3}+a_3\right) \nu_t \left( \partder{U}{z} \right)^2
!                          -\dfrac23 a_2 \nu_t \left( \partder{V}{z} \right)^2
!                          -\dfrac43 a_5 G \right)
!      \comma
!      \\[7mm]
!      \dfrac{\langle v'^2 \rangle}{k}
!      &=& \dfrac23
!      +\dfrac{1}{\mathcal{N}\eps} \left(
!        \left(\dfrac{a_2}{3}+a_3\right) \nu_t \left(\partder{V}{z}\right)^2
!                          -\dfrac23 a_2 \nu_t \left(\partder{U}{z}\right)^2
!                          -\dfrac43 a_5 G \right)
!      \comma
!      \\[7mm]
!      \dfrac{\langle w'^2 \rangle}{k}
!      &=& \dfrac23
!      +\dfrac{1}{\mathcal{N}\eps} \left(
!        \left(\dfrac{a_2}{3}-a_3\right) P
!        +\dfrac83 a_5 G
!      \right)
!      \comma
!    \end{array}
!  \end{equation}
!  where the diffusivities are computed according to \eq{nu}
!  (also see \sect{sec:cmueC} and \sect{sec:cmueD}),
!  and the buoyancy production, $G$, follows from \eq{computeG}.
!
! !USES:
  use turbulence,  only:     uu,vv,ww
  use turbulence,  only:     tke,eps,P,B,num
  use turbulence,  only:     cc1,ct1,a2,a3,a5
  IMPLICIT NONE
!
! !INPUT PARAMETERS:

! number of vertical layers
  integer,  intent(in)                 :: nlev

! square of shear frequency (1/s^2)
! (from u- and v-component)
  REALTYPE, intent(in)                :: SSU(0:nlev),SSV(0:nlev)

! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  $Log: variances.F90,v $
!  Revision 1.1.2.1  2007/09/11 13:01:21  smg
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
   REALTYPE                            :: N,Nt
   REALTYPE                            :: fac1,fac2,fac3,fac4,fac5

!
!-----------------------------------------------------------------------
!BOC

   N    =   0.5*cc1
   Nt   =   ct1

   do i=0,nlev

      fac1 = 2./3.
      fac2 = 1.0/( N*eps(i) )
      fac3 = a2/3.0 + a3
      fac4 = a2/3.0 - a3
      fac5 = 2./3.*a2

      uu(i) = tke(i)*( fac1 + fac2*( fac3*num(i)*SSU(i)                &
                          - fac5*num(i)*SSV(i) - 4./3.*a5*B(i) ) )

      vv(i) = tke(i)*( fac1 + fac2*( fac3*num(i)*SSV(i)                &
                          - fac5*num(i)*SSU(i) - 4./3.*a5*B(i) ) )

      ww(i) = tke(i)*( fac1 + fac2*( fac4*P(i) + 8./3.*a5*B(i) ) )

   enddo

   return
   end subroutine variances

!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
