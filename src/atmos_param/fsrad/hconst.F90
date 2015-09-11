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

                   MODULE HCONST_MOD

!-----------------------------------------------------------------------
!     ----- The following are physical constants -----

           REAL,PARAMETER :: AMOLWT   = 28.9644
           REAL,PARAMETER :: CSUBP    = 1.00484E7 
           REAL,PARAMETER :: DIFFCTR  = 1.66
           REAL,PARAMETER :: GRAV     = 980.665 
           REAL,PARAMETER :: GINV     = 1./GRAV 
           REAL,PARAMETER :: GRAVDR   = 980.0
           REAL,PARAMETER :: O3DIFCTR = 1.90 
           REAL,PARAMETER :: P0       = 1013250. 
           REAL,PARAMETER :: P0INV    = 1./P0 
           REAL,PARAMETER :: GP0INV   = GINV*P0INV 
           REAL,PARAMETER :: P0XZP2   = 202649.902 
           REAL,PARAMETER :: P0XZP8   = 810600.098 
           REAL,PARAMETER :: P0X2     = 2.*1013250.
           REAL,PARAMETER :: RADCON   = 8.427
           REAL,PARAMETER :: RADCON1  = 1./8.427
           REAL,PARAMETER :: RATCO2MW = 1.519449738
           REAL,PARAMETER :: RATH2OMW = 0.622 
           REAL,PARAMETER :: RGAS     = 8.3142E7 
           REAL,PARAMETER :: RGASSP   = 8.31432E7
           REAL,PARAMETER :: SECPDA   = 8.64E4 

!-----------------------------------------------------------------------

                  END MODULE HCONST_MOD

