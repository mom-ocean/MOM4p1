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

      module m_tracname_mod
!-----------------------------------------------------------
! 	... List of advected and non-advected trace species, and
!           surface fluxes for the advected species.
!-----------------------------------------------------------

      use mo_grid_mod,   only : pcnst
      use chem_mods_mod, only : grpcnt

      implicit none

character(len=128), parameter :: version     = '$Id: m_tracname.F90,v 13.0 2006/03/28 21:16:39 fms Exp $'
character(len=128), parameter :: tagname     = '$Name: mom4p1_pubrel_dec2009_nnz $'
logical                       :: module_is_initialized = .false.

      save

      character(len=8) :: tracnam(pcnst)          ! species names
      character(len=8) :: natsnam(max(1,grpcnt))  ! names of non-advected trace species

      end module m_tracname_mod
