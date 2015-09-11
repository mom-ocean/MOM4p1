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

                     MODULE RAD_DIAG_MOD

!-----------------------------------------------------------------------

      USE   RDPARM_MOD, ONLY:  LMAX, LP1, NBLW, NBLY, NBLM

      USE   HCONST_MOD, ONLY:  RADCON, RADCON1

      USE LONGWAVE_MOD, ONLY: OSOUR, CSOUR, SS1
      USE LONGWAVE_MOD, ONLY: FLX1E1, GXCTS, FCTSG
      USE LONGWAVE_MOD, ONLY: CLDFAC
      USE LONGWAVE_MOD, ONLY: DELP2, DELP
      USE LONGWAVE_MOD, ONLY: TO3, CO21, EMISS, EMISS2, CTS, EXCTS,  &
                              EXCTSN, E1FLX, CO2SP
      USE LONGWAVE_MOD, ONLY: IBAND, BANDLO, BANDHI

!     -------------------------------------------------------------
      Use       Fms_Mod, ONLY: write_version_number, mpp_pe, mpp_root_pe, &
                               error_mesg, FATAL


implicit none
private

!-----------------------------------------------------------------------
      character(len=128) :: version = '$Id: rad_diag.F90,v 10.0 2003/10/24 22:00:33 fms Exp $'
      character(len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'
      logical            :: module_is_initialized = .false.

public RADIAG, RAD_DIAG_init, RAD_DIAG_end

      CONTAINS

!#######################################################################
!#######################################################################

      SUBROUTINE RADIAG  &
          ( PRESS,TEMP,RH2O,RRVCO2,QO3,CAMT,KTOP,KBTM,NCLDS,  &
            HEATRA,GRNFLX,  &
            FSW,DFSW,UFSW,HSW,  &
            KTOPSW,KBTMSW,EMCLD,CUVRF,CIRRF,CIRAB,  &
            SALB,COSZRO,SSOLAR,   ip,jp)

!-----------------------------------------------------------------------
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRESS,TEMP,RH2O
      REAL,    INTENT(IN)                   :: RRVCO2
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: QO3,CAMT
      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOP,KBTM
      INTEGER, INTENT(IN), DIMENSION(:,:)   :: NCLDS

      REAL,    INTENT(IN), DIMENSION(:,:,:) :: HEATRA
      REAL,    INTENT(IN), DIMENSION(:,:)   :: GRNFLX

      REAL,    INTENT(IN), DIMENSION(:,:,:) :: FSW,DFSW,UFSW,HSW

      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOPSW,KBTMSW
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: EMCLD,CUVRF,CIRRF,CIRAB

      REAL,    INTENT(IN), DIMENSION(:,:)   :: SALB,COSZRO,SSOLAR
      INTEGER, INTENT(IN)                 :: ip,jp
!-----------------------------------------------------------------------

      call error_mesg('RADIAG', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE RADIAG

!#######################################################################
!#######################################################################

      subroutine RAD_DIAG_init
!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('RAD_DIAG_init', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine RAD_DIAG_init

!#######################################################################
!#######################################################################

      subroutine RAD_DIAG_end

      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('RAD_DIAG_end', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine RAD_DIAG_end


                 END MODULE RAD_DIAG_MOD

