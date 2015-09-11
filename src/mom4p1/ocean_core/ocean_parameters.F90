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
module ocean_parameters_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov">
! S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module contains some parameters used in mom4.
!</OVERVIEW>
!
!<DESCRIPTION>
! The parameters include physical constants and 
! settings for numerical and/or physical schemes.  
!</DESCRIPTION>
!

  implicit none

  private

  logical :: module_is_initialized=.false.
  character(len=128) :: version = &
     '$Id: ocean_parameters.F90,v 16.0.2.1.16.1.16.1 2009/10/14 14:37:24 smg Exp $'
  character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

  ! Boussinesq density (kg/m3) and its inverse.  Ideally, this 
  ! density should be set equal to the domain averaged density 
  ! in the model.   
  real,    parameter, public :: rho0  = 1035.0
  real,    parameter, public :: rho0r = 1.0/1035.0 

  ! specific heat capacity J/(kg degC) for seawater 
  ! from Jackett etal (2006)
  ! This value differs from that in shared/constants since 
  ! the present value is more updated.  
  real,    parameter, public :: cp_ocean = 3992.10322329649d0  

  ! specific heat capacity J/(kg degC) for calving land ice. 
  ! this value is consistent with that used in the GFDL land model. 
  real,    parameter, public :: cp_solid_runoff = 2106.

  ! specific heat capacity J/(kg degC) for liquid runoff from land.
  ! this value is consistent with that used in the GFDL land model. 
  real,    parameter, public :: cp_liquid_runoff = 4218.

  ! product of rho0*cp_ocean
  ! (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C)
  real,    parameter, public :: rho_cp = rho0*cp_ocean

  ! freezing point of fresh water at standard atmos pressure 
  real,    parameter, public :: tfreeze  = 273.15

  ! rotation of earth in radians per second 
  ! from equation (4.1) in Griffies (2004)
  real,    parameter,  public :: omega_earth = 7.2921e-5 

  ! for setting some missing values 
  real,    parameter, public :: missing_value=-1.e20

  ! some numerical constants 
  real,    parameter, public  :: onehalf      = 1.0/2.0 
  real,    parameter, public  :: onefourth    = 1.0/4.0 
  real,    parameter, public  :: onesixth     = 1.0/6.0 
  real,    parameter, public  :: oneeigth     = 1.0/8.0 
  real,    parameter, public  :: onesixteenth = 1.0/16.0 
  real,    parameter, public  :: sec_in_day   = 86400.0 
  real,    parameter, public  :: sec_in_day_r = 1.0/86400.0 
  real,    parameter, public  :: sec_in_yr    = 86400.0*365.25 
  real,    parameter, public  :: sec_in_yr_r  = 1.0/(86400.0*365.25)
  real,    parameter, public  :: m3toSv       = 1.e-6   ! volume based Sv (10^6 m^3/s)
  real,    parameter, public  :: kgtoSv       = 1.e-9   ! mass   based Sv (10^9 kg/s)

 
  ! von Karman constant (dimensionless)
  real,    parameter, public :: von_karman=0.4

  ! parameters for choosing the time tendency calculation
  integer, parameter, public :: TWO_LEVEL   = 2
  integer, parameter, public :: THREE_LEVEL = 3

  ! parameters for vertical coordinates 
  integer, parameter, public :: GEOPOTENTIAL = 1
  integer, parameter, public :: ZSTAR        = 2
  integer, parameter, public :: ZSIGMA       = 3
  integer, parameter, public :: PRESSURE     = 4
  integer, parameter, public :: PSTAR        = 5 
  integer, parameter, public :: PSIGMA       = 6 

  integer, parameter, public :: DEPTH_BASED    = 1
  integer, parameter, public :: PRESSURE_BASED = 2

  integer, parameter, public :: QUASI_HORIZONTAL  = 1
  integer, parameter, public :: TERRAIN_FOLLOWING = 2

  ! parameters for choosing method of computing vertical cell thickness 
  integer, parameter, public :: ENERGETIC   =1
  integer, parameter, public :: FINITEVOLUME=2

  ! parameters for choosing the temperature variable
  integer, parameter, public :: CONSERVATIVE_TEMP = 1
  integer, parameter, public :: POTENTIAL_TEMP    = 2

  ! parameters for choosing the vertical mixing scheme 
  integer, parameter, public :: VERTMIX_CONST      = 1
  integer, parameter, public :: VERTMIX_KPP        = 2
  integer, parameter, public :: VERTMIX_KPP_MOM4P0 = 3
  integer, parameter, public :: VERTMIX_PP         = 4
  integer, parameter, public :: VERTMIX_CHEN       = 5
  integer, parameter, public :: VERTMIX_GOTM       = 6

  ! parameters for grid type
  integer, parameter, public :: TGRID = 1
  integer, parameter, public :: UGRID = 2
  
  integer, parameter, public :: TEMP_ID = 1
  integer, parameter, public :: SALT_ID = 2

  ! parameters for tracer advection 
  integer, parameter, public :: ADVECT_UPWIND     = 1
  integer, parameter, public :: ADVECT_2ND_ORDER  = 2
  integer, parameter, public :: ADVECT_4TH_ORDER  = 3
  integer, parameter, public :: ADVECT_6TH_ORDER  = 4 
  integer, parameter, public :: ADVECT_QUICKER    = 5
  integer, parameter, public :: ADVECT_QUICKMOM3  = 6
  integer, parameter, public :: ADVECT_MDFL_SUP_B = 7
  integer, parameter, public :: ADVECT_MDFL_SWEBY = 8
  integer, parameter, public :: ADVECT_PSOM       = 9
  integer, parameter, public :: ADVECT_MDPPM      = 10
  integer, parameter, public :: ADVECT_DST_LINEAR = 11

  ! parameters for linear mometum  advection 
  integer, parameter, public :: VEL_ADVECT_UPWIND     = 1
  integer, parameter, public :: VEL_ADVECT_2ND_ORDER  = 2

end module ocean_parameters_mod


