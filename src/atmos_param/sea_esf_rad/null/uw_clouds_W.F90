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
!FDOC_TAG_GFDL

                 module uw_clouds_W_mod
! <CONTACT EMAIL="fei.liu@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!          uw shallow convection cloud radiative properties module
!   
! </OVERVIEW>
! <DESCRIPTION>
!   
! </DESCRIPTION>
!

use time_manager_mod,       only: time_type
use       fms_mod,      only:  error_mesg, FATAL, &
                               mpp_pe, mpp_root_pe, &
                               write_version_number
use rad_utilities_mod,      only: microphysics_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!          uw shallow convection cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id: uw_clouds_W.F90,v 15.0 2007/08/14 03:55:19 fms Exp $'
   character(len=128)  :: tagname =  '$Name: mom4p1_pubrel_dec2009_nnz $'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          uw_clouds_W_init,   &
          uw_clouds_W_end , uw_clouds_amt

!---------------------------------------------------------------------
!-------- namelist  ---------

logical   :: dummy = .true.


namelist /uw_clouds_W_nml /     &
                                     dummy                          


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

  logical :: module_is_initialized = .false.
!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 





! <SUBROUTINE NAME="uw_clouds_W_init">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call uw_clouds_W_init  (pref, lonb, latb, axes, Time)
!		
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
! 
!  </IN>
!  <IN NAME="lonb" TYPE="real">
! 
!  </IN>
!  <IN NAME="latb" TYPE="real">
! 
!  </IN>
!  <IN NAME="axes" TYPE="integer">
! 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
! 
!  </IN>
! </SUBROUTINE>
!
subroutine uw_clouds_W_init  (pref, lonb, latb, axes, Time)

real, dimension(:,:), intent(in) :: pref
real, dimension(:,:), intent(in) :: lonb, latb
integer, dimension(4), intent(in)      :: axes
type(time_type),       intent(in)      :: Time

      integer            :: unit, ierr, io


      if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
      endif

      module_is_initialized = .true.

      call error_mesg('uw_clouds_W_init', &
      'This module is not supported as part of the public release', FATAL)


end subroutine uw_clouds_W_init

! <SUBROUTINE NAME="uw_clouds_W_end">
!  <OVERVIEW>
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call uw_clouds_W_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine uw_clouds_W_end
       
!----------------------------------------------------------------------
!    uw_clouds_W_end is the destructor for uw_clouds_W_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!--------------------------------------------------------------------

      call error_mesg('uw_clouds_W_end', &
      'This module is not supported as part of the public release', FATAL)

end subroutine uw_clouds_W_end


!#################################################################


!---------------------------------------------------------------------

! <SUBROUTINE NAME="uw_clouds_amt">
!  <OVERVIEW>
!    uw_clouds_amt defines the distribution of cloud water and cloud ice 
!    amounts [ g / m**3 ] and liquid and ice particle sizes and total cloud 
!    fraction for the clouds associated with uw shallow convection. these 
!    values will later be combined with other cloud fields to produce the 
!    cloud radiative properties that will be seen by the radiation package.
!   
!  </OVERVIEW>
!  <DESCRIPTION>
!    uw_clouds_amt defines the distribution of cloud water and cloud ice 
!    amounts [ g / m**3 ] and liquid and ice particle sizes and total cloud 
!    fraction for the clouds associated with uw shallow convection. these 
!    values will later be combined with other cloud fields to produce the 
!    cloud radiative properties that will be seen by the radiation package.
!   
!  </DESCRIPTION>
!  <TEMPLATE>
!   call uw_clouds_amt (is, ie, js, je, Shallow_microphys)   
!		
!  </TEMPLATE>
!  <IN NAME="is" TYPE="integer">
! 
!  </IN>
!  <IN NAME="ie" TYPE="integer">
! 
!  </IN>
!  <IN NAME="js" TYPE="integer">
! 
!  </IN>
!  <IN NAME="je" TYPE="integer">
! 
!  </IN>
!  <INOUT NAME="Shallow_microphys" TYPE="microphysics_type">
! 
!  </INOUT>
! </SUBROUTINE>
!
subroutine uw_clouds_amt (is, ie, js, je,   &
                   shallow_cloud_area, shallow_liquid, shallow_ice, &
                   shallow_droplet_number, land, pfull, tkel, &
                   Shallow_microphys)

!---------------------------------------------------------------------
!    uw_clouds_amt defines the distribution of cloud water and cloud ice 
!    amounts [ g / m**3 ] and liquid and ice particle sizes and total cloud 
!    fraction for the clouds associated with uw shallow convection. these 
!    values will later be combined with other cloud fields to produce the 
!    cloud radiative properties that will be seen by the radiation package.
!----------------------------------------------------------------------

integer,                 intent(in)    :: is,ie,js,je
real, dimension(:,:,:),  intent(in)    :: shallow_cloud_area,  &
                                          shallow_liquid, shallow_ice, &
                                          shallow_droplet_number
real, dimension(:,:),    intent(in)    :: land
real, dimension(:,:,:),  intent(in)    :: pfull, tkel
type(microphysics_type), intent(inout) :: Shallow_microphys


      call error_mesg('uw_clouds_amt', &
      'This module is not supported as part of the public release', FATAL)

end subroutine uw_clouds_amt  



!####################################################################


                     end module uw_clouds_W_mod

