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

                 module specified_clouds_W_mod

use time_manager_mod,   only:  time_type
use       fms_mod,      only:  error_mesg, FATAL, &
                               mpp_pe, mpp_root_pe, &
                               write_version_number
use rad_utilities_mod,  only:  cld_specification_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!             specified clouds radiative properties module;
!             used with cloud_obs_mod and cloud_zonal_mod
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128)  :: version =  '$Id: specified_clouds_W.F90,v 15.0 2007/08/14 03:55:16 fms Exp $'
  character(len=128)  :: tagname =  '$Name: mom4p1_pubrel_dec2009_nnz $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          specified_clouds_W_init, specified_clouds_amt, &
          specified_clouds_W_end

logical  :: module_is_initialized = .false.

!------------------------------------------------------------------
!------------------------------------------------------------------



contains 

subroutine specified_clouds_W_init (lonb, latb)


real, dimension(:,:), intent(in) :: lonb, latb


      integer          :: unit, ierr, io

!------- write version number and namelist ---------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

!---------------------------------------------------------------------

      call error_mesg('specified_clouds_W_init', &
      'This module is not supported as part of the public release', FATAL)

end subroutine specified_clouds_W_init

subroutine specified_clouds_W_end
        
!----------------------------------------------------------------------
!    specified_clouds_end is the destructor for specified_clouds_W_mod.
!----------------------------------------------------------------------
        
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------

      call error_mesg('specified_clouds_W_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine specified_clouds_W_end


!######################################################################

subroutine specified_clouds_amt (is, ie, js, je, Rad_time, lat, pflux, &
                                 Cld_spec)

!----------------------------------------------------------------------
!    specified_clouds_amt defines the location, amount (cloud fraction),
!    number and type (hi, mid, low) of clouds present on the model grid.
!----------------------------------------------------------------------
 
!--------------------------------------------------------------------
integer,                      intent(in)    :: is, ie, js, je
type(time_type),              intent(in)    :: Rad_time
real, dimension(:,:),         intent(in)    :: lat
real, dimension(:,:,:),       intent(in)    :: pflux
type(cld_specification_type), intent(inout) :: Cld_spec
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in 
!                   the physics_window being integrated
!      Rad_time     time at which the climatologically-determined, 
!                   time-varying specified cloud fields should apply
!                   [ time_type, days and seconds]
!      lat          latitude of model points  [ radians ]
!      pflux        average of pressure at adjacent model levels
!                   [ (kg /( m s^2) ] 
!
!   intent(inout) variables:
!
!      Cld_spec     cld_specification_type variable containing the 
!                   cloud specification input fields needed by the 
!                   radiation package
!
!               the following elements of Cld_spec are defined here:
!
!                  %cmxolw  fraction of maximally overlapped clouds
!                           seen by the longwave radiation 
!                           [ dimensionless ]
!                  %crndlw  fraction of randomly overlapped clouds
!                           seen by the longwave radiation 
!                           [ dimensionless ]
!                  %camtsw  cloud fraction seen by the shortwave
!                           radiation; the sum of the maximally
!                           overlapped and randomly overlapped 
!                           longwave cloud fractions  [ dimensionless ]
!                  %nmxolw  number of maximally overlapped longwave 
!                           clouds in each grid column.
!                  %nrndlw  number of randomly overlapped longwave 
!                           clouds in each grid column.
!                  %ncldsw  number of clouds seen by he shortwave
!                           radiation in each grid column.
!                  %hi_cld  logical flag indicating the presence of
!                           high clouds in a grid box
!                 %mid_cld  logical flag indicating the presence of 
!                           middle clouds in a grid box
!                 %low_cld  logical flag indicating the presence of 
!                           low clouds in a grid box
!                                                                  
!---------------------------------------------------------------------

      call error_mesg('specified_clouds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine specified_clouds_amt 

!######################################################################

                 end module specified_clouds_W_mod



