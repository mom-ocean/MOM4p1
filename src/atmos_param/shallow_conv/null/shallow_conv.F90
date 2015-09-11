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
  MODULE SHALLOW_CONV_MOD

!=======================================================================
! --- SHALLOW CONVECTION MODULE - GFDL SPECTRAL MODEL VERSION
!=======================================================================

 use  Sat_Vapor_Pres_Mod, ONLY: ESCOMP, DESCOMP
 use       Fms_Mod,       ONLY: FILE_EXIST, ERROR_MESG, FATAL,   &
                                CHECK_NML_ERROR, OPEN_NAMELIST_FILE,      &
                                CLOSE_FILE, mpp_pe, mpp_root_pe, &
                                write_version_number, stdlog

 use constants_mod, only: Hlv, Cp_Air, RDgas, RVgas, Kappa, grav

!---------------------------------------------------------------------
 implicit none
 private
!---------------------------------------------------------------------

 public  :: SHALLOW_CONV, SHALLOW_CONV_INIT, SHALLOW_CONV_END
 public  :: MYLCL

!---------------------------------------------------------------------

 character(len=128) :: version = '$Id: shallow_conv.F90,v 10.0 2003/10/24 22:00:49 fms Exp $'
 character(len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'

 logical :: module_is_initialized = .false.

!---------------------------------------------------------------------

 contains

!#######################################################################
!#######################################################################

 SUBROUTINE SHALLOW_CONV_INIT( kx )

!=======================================================================
! ***** INITIALIZE SHALLOW CONVECTION
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!     kx     - Number of levels in vertical
!---------------------------------------------------------------------
 integer, intent(in) :: kx
 
!------- write version number and namelist ---------

  if ( mpp_pe() == mpp_root_pe() ) then
       call write_version_number(version, tagname)
  endif

!-------------------------------------------------------------------
  module_is_initialized = .true.
!---------------------------------------------------------------------

      call error_mesg('SHALLOW_CONV_INIT', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE SHALLOW_CONV_INIT

!#######################################################################

  SUBROUTINE SHALLOW_CONV_END
!-------------------------------------------------------------------
  module_is_initialized = .false.
!---------------------------------------------------------------------

      call error_mesg('SHALLOW_CONV_END', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE SHALLOW_CONV_END

!#######################################################################

  SUBROUTINE SHALLOW_CONV( Temp, qmix0, pfull, phalf, akhsc, kbot )

!=======================================================================
! --- SHALLOW CONVECTION
!=======================================================================
!----------------------------------------------------------------------
! Arguments (Intent in)
!       Temp    -  Temperature
!       qmix0   -  Specific humidity
!       pfull   -  Pressure at full levels
!       phalf   -  Pressure at half levels
!       kbot    -  OPTIONAL; lowest model level index (integer)
!----------------------------------------------------------------------
  real, intent(in), dimension(:,:,:) :: Temp, qmix0, pfull, phalf

  integer, intent(in), OPTIONAL, dimension(:,:) :: kbot

!----------------------------------------------------------------------
! Arguments (Intent out)
!       akhsc  -  mixing coefficient for heat and moisture
!                 due to shallow convection
!----------------------------------------------------------------------
  real, intent(out), dimension(:,:,:) :: akhsc

!---------------------------------------------------------------------

      call error_mesg('SHALLOW_CONV', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE SHALLOW_CONV

!#######################################################################

  SUBROUTINE MYLCL ( tlparc, qlparc, plparc, phalf, plcl, kbase )

!=======================================================================
! ***** COMPUTE LCL ( CLOUD BASE )
!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!       tlparc   Initial parcel temperature
!       qlparc   Initial parcel mixing ratio
!       plparc   Initial parcel pressure
!       phalf    Pressure at half levels
! Arguments (Intent out)
!       plcl     Pressure at LCL
!       kbase    Index of LCL in column
!---------------------------------------------------------------------
  real,    intent(in),  dimension(:,:)   :: tlparc, qlparc, plparc
  real,    intent(in),  dimension(:,:,:) :: phalf
  real,    intent(out), dimension(:,:)   :: plcl
  integer, intent(out), dimension(:,:)   :: kbase

!---------------------------------------------------------------------

      call error_mesg('MYLCL', &
      'This module is not supported as part of the public release', FATAL)

  end SUBROUTINE MYLCL

!#######################################################################
!#######################################################################
  end MODULE SHALLOW_CONV_MOD

