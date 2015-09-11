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

! $Id: cosp_simulator.f90,v 1.1.2.1.2.1 2009/08/10 10:44:31 rsh Exp $
! $Name: mom4p1_pubrel_dec2009_nnz $

! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
!
!

MODULE MOD_COSP_SIMULATOR
  USE MOD_COSP_TYPES
  USE MOD_COSP_RADAR
  USE MOD_COSP_LIDAR
  USE MOD_COSP_ISCCP_SIMULATOR
  USE MOD_COSP_MISR_SIMULATOR
  USE MOD_COSP_STATS
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_SIMULATOR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_SIMULATOR(me, gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,stradar,stlidar)

  ! Arguments
  integer, intent(in) :: me
  type(cosp_gridbox),intent(in) :: gbx      ! Grid-box inputs
  type(cosp_subgrid),intent(in) :: sgx      ! Subgrid inputs
  type(cosp_sghydro),intent(in) :: sghydro  ! Subgrid info for hydrometeors
  type(cosp_config),intent(in) :: cfg       ! Configuration options
  type(cosp_vgrid),intent(in)   :: vgrid    ! Information on vertical grid of stats
  type(cosp_sgradar),intent(inout) :: sgradar ! Output from radar simulator
  type(cosp_sglidar),intent(inout) :: sglidar ! Output from lidar simulator
  type(cosp_isccp),intent(inout)   :: isccp   ! Output from ISCCP simulator
  type(cosp_misr),intent(inout)    :: misr    ! Output from MISR simulator
  type(cosp_radarstats),intent(inout) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics from lidar simulator
  ! Local variables
  ! ***Timing variables (to be deleted in final version)
  integer :: t0,t1,count_rate,count_max

  !+++++++++ Radar model ++++++++++  
  if (cfg%Lradar_sim) then
!   call system_clock(t0,count_rate,count_max)
    call cosp_radar(me, gbx,sgx,sghydro,sgradar)
!   call system_clock(t1,count_rate,count_max)
!   print *, '%%%%%%  Radar:', (t1-t0)*1.0/count_rate, ' s'
! else 
!   print *, '%%%%%%  Radar not used'
  endif
  
  !+++++++++ Lidar model ++++++++++
  if (cfg%Llidar_sim) then
!   call system_clock(t0,count_rate,count_max)
    call cosp_lidar(gbx,sgx,sghydro,sglidar)
!   call system_clock(t1,count_rate,count_max)
!   print *, '%%%%%%  Lidar:', (t1-t0)*1.0/count_rate, ' s'
! else 
!   print *, '%%%%%%  Lidar not used'
  endif

  
  !+++++++++ ISCCP simulator ++++++++++
  if (cfg%Lisccp_sim) then
!   call system_clock(t0,count_rate,count_max)
    call cosp_isccp_simulator(gbx,sgx,isccp)
!   call system_clock(t1,count_rate,count_max)
!   print *, '%%%%%%  ISCCP:', (t1-t0)*1.0/count_rate, ' s'
! else 
!   print *, '%%%%%%  ISCCP not used'
  endif
  
  !+++++++++ MISR simulator ++++++++++
  if (cfg%Lmisr_sim) then
!   call system_clock(t0,count_rate,count_max)
    call cosp_misr_simulator(gbx,sgx,misr)
!   call system_clock(t1,count_rate,count_max)
!   print *, '%%%%%%  MISR:', (t1-t0)*1.0/count_rate, ' s'
! else 
!   print *, '%%%%%%  MISR not used'
  endif
  

  !+++++++++++ Summary statistics +++++++++++
!   write(*,*) 'Stats:'
!   read(*,*) c 
  if (cfg%Lstats) then
!   call system_clock(t0,count_rate,count_max)
    call cosp_stats(gbx,sgx,cfg,sgradar,sglidar,vgrid,stradar,stlidar)
!   call system_clock(t1,count_rate,count_max)
!   print *, '%%%%%%  Stats:', (t1-t0)*1.0/count_rate, ' s'
  endif
  !+++++++++++ change of units after computation of statistics +++++++++++
  if (cfg%Llidar_sim) then
    where((sglidar%beta_tot > 0.0) .and. (sglidar%beta_tot /= R_UNDEF)) 
        sglidar%beta_tot = log10(sglidar%beta_tot)
    elsewhere
        sglidar%beta_tot = R_UNDEF
    end where
  endif

END SUBROUTINE COSP_SIMULATOR

END MODULE MOD_COSP_SIMULATOR
