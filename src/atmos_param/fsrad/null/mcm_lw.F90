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
      MODULE MCM_LW_MOD

!   Added interface routine (lw_rad_ss) which is called by
!     fsrad and which calls lwcool in this
!     module after constructing the appropriate inputs.

      USE Constants_Mod, ONLY: grav, tfreeze

      Use       Fms_Mod, ONLY: Error_Mesg, FATAL, &
                               write_version_number, mpp_pe, mpp_root_pe

implicit none
      private

!------------ VERSION NUMBER ----------------

      character(len=128) :: version = '$Id: mcm_lw.F90,v 10.0 2003/10/24 22:00:33 fms Exp $'
      character(len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'
      logical            :: module_is_initialized = .false.

      public  MCM_LW_RAD, mcm_lw_init, mcm_lw_end

!     -------------------------------------------------

!-----------------------------------------------------------------------
!--------------------- G L O B A L   D A T A ---------------------------
!-----------------------------------------------------------------------


      contains

!#######################################################################
      subroutine mcm_lw_init(ix_in, jx_in, kx_in)
      integer, intent(in) :: ix_in, jx_in, kx_in

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('mcm_lw_init', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine mcm_lw_init
!#######################################################################

      subroutine MCM_LW_END

      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('MCM_LW_END', &
      'This module is not supported as part of the public release', FATAL)

      end subroutine MCM_LW_END

!#######################################################################
      SUBROUTINE MCM_LW_RAD (KTOP,KBTM,NCLDS,EMCLD, &
                      PRES,TEMP,RH2O,QO3,CAMT, &
                      RRVCO2,  HEATRA,GRNFLX,TOPFLX, phalf)

      INTEGER, INTENT(IN), DIMENSION(:,:,:) :: KTOP,KBTM
      INTEGER, INTENT(IN), DIMENSION(:,:)    :: NCLDS
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: EMCLD
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: PRES,TEMP
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: RH2O,QO3
      REAL,    INTENT(IN), DIMENSION(:,:,:) :: CAMT
      REAL,    INTENT(IN)                      :: RRVCO2
 
      REAL,   INTENT(OUT), DIMENSION(:,:,:) :: HEATRA
      REAL,   INTENT(OUT), DIMENSION(:,:)    :: GRNFLX,TOPFLX

      REAL,    INTENT(IN), DIMENSION(:,:,:) :: phalf



!---------------------------------------------------------------------

      call error_mesg('MCM_LW_RAD', &
      'This module is not supported as part of the public release', FATAL)

      END SUBROUTINE MCM_LW_RAD

      end module mcm_lw_mod
