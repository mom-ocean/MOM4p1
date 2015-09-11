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

!---------------------------------------------------------------------
!------------ FMS version number and tagname for this file -----------
        
! $Id: load_hydrometeor_classes.f90,v 1.1.2.1.2.1 2009/08/10 10:48:14 rsh Exp $
! $Name: mom4p1_pubrel_dec2009_nnz $

  subroutine load_hydrometeor_classes(Nprmts_max,dist_prmts_hydro,hp,nhclass)
  use radar_simulator_types
  implicit none
  
! Purpose:
!   Loads the hydrometeor classes to be used in calculations
!   Part of QuickBeam v1.03 by John Haynes
!   http://reef.atmos.colostate.edu/haynes/radarsim
!
! Inputs:  
!   [dist_prmts_hydro]   from data in hydrometeor class input 
!
! Outputs:
!   [hp]            structure that define hydrometeor types
!
! Modified:
!   08/23/2006  placed into subroutine form (Roger Marchand)
   
! ----- INPUT -----
  integer, intent(in) :: nhclass,Nprmts_max
  real,dimension(Nprmts_max,nhclass), intent(in) :: dist_prmts_hydro
! ----- OUTPUTS -----  
  type(class_param), intent(out) :: hp
  
! ----- INTERNAL -----  
  integer :: i
	
    hp%rho(:) = -1

    do i = 1,nhclass,1
    hp%dtype(i) = dist_prmts_hydro(1,i)
    hp%col(i) = dist_prmts_hydro(2,i)
    hp%phase(i) = dist_prmts_hydro(3,i)
    hp%cp(i) = dist_prmts_hydro(4,i)
    hp%dmin(i) = dist_prmts_hydro(5,i)
    hp%dmax(i) = dist_prmts_hydro(6,i)
    hp%apm(i) = dist_prmts_hydro(7,i)
    hp%bpm(i) = dist_prmts_hydro(8,i)
    hp%rho(i) = dist_prmts_hydro(9,i)
    hp%p1(i) = dist_prmts_hydro(10,i)
    hp%p2(i) = dist_prmts_hydro(11,i)
    hp%p3(i) = dist_prmts_hydro(12,i)
    enddo
        
!   // setup scaling arrays
    hp%fc = -999.
    hp%scaled = .false.	
    hp%z_flag = .false.
    hp%rho_eff = -999.
    hp%ifc = -9
    hp%idd = -9
   
  
  end subroutine load_hydrometeor_classes
