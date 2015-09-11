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

! $Id: lmd_ipsl_stats.f90,v 1.1.2.1.2.1 2009/08/10 10:45:27 rsh Exp $
! $Name: mom4p1_pubrel_dec2009_nnz $

! Copyright (c) 2009, Centre National de la Recherche Scientifique
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
!     * Neither the name of the LMD/IPSL/CNRS/UPMC nor the names of its
!       contributors may be used to endorse or promote products derived from this software without 
!       specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


!------------------------------------------------------------------------------------
! Authors: Sandrine Bony and Helene Chepfer (LMD/IPSL, CNRS, UPMC, France).
!------------------------------------------------------------------------------------
MODULE MOD_LMD_IPSL_STATS
  USE MOD_LLNL_STATS
  use mod_cosp_constants
  IMPLICIT NONE

!RSH made module variables so can be accessed from cosp_driver_init for
! use in defining netcdf axes.
! c threshold for cloud detection :
      real S_clr 
      parameter (S_clr = 1.2) 
      real S_cld
!      parameter (S_cld = 3.0)  ! Previous thresold for cloud detection
      parameter (S_cld = 5.0)  ! New (dec 2008) thresold for cloud detection
      real S_att
      parameter (S_att = 0.01)

CONTAINS
      SUBROUTINE diag_lidar(npoints,ncol,llm,max_bin,nrefl &
                  ,pnorm,pmol,refl,land,pplay,undef,ok_lidar_cfad &
                  ,cfad2,srbval &
                  ,ncat,lidarcld,cldlayer,parasolrefl)
!
! -----------------------------------------------------------------------------------
! Lidar outputs :
! 
! Diagnose cloud fraction (3D cloud fraction + low/middle/high/total cloud fraction
! from the lidar signals (ATB and molecular ATB) computed from model outputs
!      +
! Compute CFADs of lidar scattering ratio SR and of depolarization index
! 
! Authors: Sandrine Bony and Helene Chepfer (LMD/IPSL, CNRS, UPMC, France).
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - change of the cloud detection threshold S_cld from 3 to 5, for better
! with both day and night observations. The optical thinest clouds are missed.
! - remove of the detection of the first fully attenuated layer encountered from above.
! December 2008, A. Bodas-Salcedo:
! - Dimensions of pmol reduced to (npoints,llm)
!
! Version 1.0 (June 2007)
! Version 1.1 (May 2008)
! Version 1.2 (June 2008)
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
! c------------------------------------------------------------------------------------

! c inputs :
      integer npoints
      integer ncol
      integer llm
      integer max_bin               ! nb of bins for SR CFADs
      integer ncat                  ! nb of cloud layer types (low,mid,high,total)
      integer nrefl                 ! nb of solar zenith angles for parasol reflectances

      real undef                    ! undefined value
      real pnorm(npoints,ncol,llm)  ! lidar ATB 
      real pmol(npoints,llm)        ! molecular ATB
      real land(npoints)            ! Landmask [0 - Ocean, 1 - Land]    
      real pplay(npoints,llm)       ! pressure on model levels (Pa)
      logical ok_lidar_cfad         ! true if lidar CFAD diagnostics need to be computed
      real refl(npoints,ncol,nrefl) ! subgrid parasol reflectance ! parasol

! c outputs :
      real lidarcld(npoints,llm)     ! 3D "lidar" cloud fraction 
      real cldlayer(npoints,ncat)    ! "lidar" cloud fraction (low, mid, high, total)
      real cfad2(npoints,max_bin,llm) ! CFADs of SR  
      real srbval(max_bin)           ! SR bins in CFADs  
      real parasolrefl(npoints,nrefl)! grid-averaged parasol reflectance


! c local variables :
      integer ic,k
      real x3d(npoints,ncol,llm)
      real x3d_c(npoints,llm),pnorm_c(npoints,llm)
      real xmax
!
! c -------------------------------------------------------
! c 0- Initializations
! c -------------------------------------------------------
!
! Parasol reflectance algorithm is not valid over land. Write
! a warning if there is no land. Landmask [0 - Ocean, 1 - Land] 
!RSH  IF ( MAXVAL(land(:)) .EQ. 0.0) THEN
!RSH      WRITE (*,*) 'WARNING. PARASOL reflectance is not valid over land' &
!RSH         & ,' and there is only land'
!RSH    END IF

!  Should be modified in future version
      xmax=undef-1.0

! c -------------------------------------------------------
! c 1- Lidar scattering ratio :
! c -------------------------------------------------------
!
!       where ((pnorm.lt.xmax) .and. (pmol.lt.xmax) .and. (pmol.gt. 0.0 )) 
!          x3d = pnorm/pmol
!       elsewhere
!           x3d = undef
!       end where
! A.B-S: pmol reduced to 2D (npoints,llm) (Dec 08)
      do ic = 1, ncol
        pnorm_c = pnorm(:,ic,:)
        where ((pnorm_c.lt.xmax) .and. (pmol.lt.xmax) .and. (pmol.gt. 0.0 )) 
            x3d_c = pnorm_c/pmol
        elsewhere
            x3d_c = undef
        end where
        x3d(:,ic,:) = x3d_c
      enddo

! c -------------------------------------------------------
! c 2- Diagnose cloud fractions (3D, low, middle, high, total)
! c from subgrid-scale lidar scattering ratios :
! c -------------------------------------------------------

      CALL COSP_CLDFRAC(npoints,ncol,llm,ncat,  &
              x3d,pplay, S_att,S_cld,undef,lidarcld, &
              cldlayer)

! c -------------------------------------------------------
! c 3- CFADs 
! c -------------------------------------------------------
      if (ok_lidar_cfad) then
!
! c CFADs of subgrid-scale lidar scattering ratios :
! c -------------------------------------------------------
      CALL COSP_CFAD_SR(npoints,ncol,llm,max_bin, &
                 x3d, &
                 S_att,S_clr,xmax,cfad2,srbval)

      endif   ! ok_lidar_cfad
! c -------------------------------------------------------

! c -------------------------------------------------------
! c 4- Compute grid-box averaged Parasol reflectances
! c -------------------------------------------------------

      parasolrefl(:,:) = 0.0

      do k = 1, nrefl
       do ic = 1, ncol
         parasolrefl(:,k) = parasolrefl(:,k) + refl(:,ic,k)
       enddo
      enddo

      do k = 1, nrefl
        parasolrefl(:,k) = parasolrefl(:,k) / float(ncol)
! if land=1 -> parasolrefl=undef
! if land=0 -> parasolrefl=parasolrefl
        parasolrefl(:,k) = parasolrefl(:,k) * MAX(1.0-land(:),0.0) &
                           + (1.0 - MAX(1.0-land(:),0.0))*undef 
      enddo

      RETURN
      END SUBROUTINE diag_lidar
	  
	  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- FUNCTION COSP_CFAD_SR ------------------------
! Author: Sandrine Bony (LMD/IPSL, CNRS, Paris)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE COSP_CFAD_SR(Npoints,Ncolumns,Nlevels,Nbins, &
                      x,S_att,S_clr,xmax,cfad,srbval)
      IMPLICIT NONE

!--- Input arguments
! Npoints: Number of horizontal points
! Ncolumns: Number of subcolumns
! Nlevels: Number of levels
! Nbins: Number of x axis bins
! xmax: maximum value allowed for x
! S_att: Threshold for full attenuation
! S_clr: Threshold for clear-sky layer
!
!--- Input-Outout arguments
! x: variable to process (Npoints,Ncolumns,Nlevels), mofified where saturation occurs
!
! -- Output arguments
! srbval : values of the histogram bins
! cfad: 2D histogram on each horizontal point

! Input arguments
      integer Npoints,Ncolumns,Nlevels,Nbins
      real xmax,S_att,S_clr
! Input-output arguments
      real x(Npoints,Ncolumns,Nlevels)
! Output :
      real cfad(Npoints,Nbins,Nlevels)
      real srbval(Nbins)
! Local variables
      integer i, j, k, ib

! c -------------------------------------------------------
! c 0- Initializations
! c -------------------------------------------------------
      if ( Nbins .lt. 6) return

      call define_srbval (srbval)



      cfad(:,:,:) = 0.0


! c -------------------------------------------------------
! c c- Compute CFAD
! c -------------------------------------------------------

        do j = Nlevels, 1, -1 
          do k = 1, Ncolumns
              where ( x(:,k,j).le.srbval(1) ) &
                        cfad(:,1,j) = cfad(:,1,j) + 1.0
          enddo  !k
        enddo  !j

      do ib = 2, Nbins
        do j = Nlevels, 1, -1 
          do k = 1, Ncolumns
              where ( ( x(:,k,j).gt.srbval(ib-1) ) .and. ( x(:,k,j).le.srbval(ib) ) ) &
                        cfad(:,ib,j) = cfad(:,ib,j) + 1.0
          enddo  !k
        enddo  !j
      enddo  !ib

      cfad(:,:,:) = cfad(:,:,:) / float(Ncolumns)

! c -------------------------------------------------------
      RETURN
      END SUBROUTINE COSP_CFAD_SR


subroutine define_srbval (srbval)

real, dimension(:), intent(out) :: srbval

     integer :: i, Nbins

     Nbins = SR_BINS

     srbval(1) =  S_att
     srbval(2) =  S_clr
     srbval(3) =  3.0
     srbval(4) =  5.0
     srbval(5) =  7.0
     srbval(6) = 10.0
     do i = 7, MIN(10,Nbins)
      srbval(i) = srbval(i-1) + 5.0
     enddo
     DO i = 11, MIN(13,Nbins)
      srbval(i) = srbval(i-1) + 10.0
     enddo
     srbval(MIN(14,Nbins)) = 80.0
     srbval(Nbins) = LIDAR_UNDEF - 1.0

end subroutine define_srbval



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- SUBROUTINE COSP_CLDFRAC -------------------
! c Purpose: Cloud fraction diagnosed from lidar measurements 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE COSP_CLDFRAC(Npoints,Ncolumns,Nlevels,Ncat, &
                  x,pplay,S_att,S_cld,undef,lidarcld, &
                  cldlayer)
      IMPLICIT NONE
! Input arguments
      integer Npoints,Ncolumns,Nlevels,Ncat
      real x(Npoints,Ncolumns,Nlevels)
      real pplay(Npoints,Nlevels)
      real S_att,S_cld
      real undef
! Output :
      real lidarcld(Npoints,Nlevels) ! 3D cloud fraction
      real cldlayer(Npoints,Ncat)    ! low, middle, high, total cloud fractions
! Local variables
      integer ip, k, iz, ic
      real p1
      real cldy(Npoints,Ncolumns,Nlevels)
      real srok(Npoints,Ncolumns,Nlevels)
      real cldlay(Npoints,Ncolumns,Ncat)
      real nsublay(Npoints,Ncolumns,Ncat), nsublayer(Npoints,Ncat)
      real nsub(Npoints,Nlevels)


! ---------------------------------------------------------------
! 1- initialization 
! ---------------------------------------------------------------

      if ( Ncat .ne. 4 ) then
         print *,'Error in lmd_ipsl_stats.cosp_cldfrac, Ncat must be 4, not',Ncat
         stop
      endif

      lidarcld = 0.0
      nsub = 0.0
      cldlay = 0.0
      nsublay = 0.0

! ---------------------------------------------------------------
! 2- Cloud detection
! ---------------------------------------------------------------

      do k = 1, Nlevels

! cloud detection at subgrid-scale:
         where ( (x(:,:,k).gt.S_cld) .and. (x(:,:,k).ne. undef) )
           cldy(:,:,k)=1.0
         elsewhere
           cldy(:,:,k)=0.0
         endwhere

! number of usefull sub-columns:
         where ( (x(:,:,k).gt.S_att) .and. (x(:,:,k).ne. undef)  ) 
           srok(:,:,k)=1.0
         elsewhere
           srok(:,:,k)=0.0
         endwhere

      enddo ! k

! ---------------------------------------------------------------
! 3- grid-box 3D cloud fraction and layered cloud fractions (ISCCP pressure
! categories) :
! ---------------------------------------------------------------

      do k = Nlevels, 1, -1
       do ic = 1, Ncolumns
        do ip = 1, Npoints

          iz=1
          p1 = pplay(ip,k)
          if ( p1.gt.0. .and. p1.lt.(440.*100.)) then ! high clouds
            iz=3
          else if(p1.ge.(440.*100.) .and. p1.lt.(680.*100.)) then  ! mid clouds
            iz=2
         endif

         cldlay(ip,ic,iz) = MAX(cldlay(ip,ic,iz),cldy(ip,ic,k))
         cldlay(ip,ic,4) = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))
         lidarcld(ip,k)=lidarcld(ip,k) + cldy(ip,ic,k)

         nsublay(ip,ic,iz) = MAX(nsublay(ip,ic,iz),srok(ip,ic,k))
         nsublay(ip,ic,4) = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
         nsub(ip,k)=nsub(ip,k) + srok(ip,ic,k)

        enddo
       enddo
      enddo

! -- grid-box 3D cloud fraction

      where ( nsub(:,:).gt.0.0 )
         lidarcld(:,:) = lidarcld(:,:)/nsub(:,:)
      elsewhere
         lidarcld(:,:) = undef
      endwhere

! -- layered cloud fractions

      cldlayer = 0.0
      nsublayer = 0.0

      do iz = 1, Ncat
       do ic = 1, Ncolumns

          cldlayer(:,iz)=cldlayer(:,iz) + cldlay(:,ic,iz)    
          nsublayer(:,iz)=nsublayer(:,iz) + nsublay(:,ic,iz) 

       enddo
      enddo
      where ( nsublayer(:,:).gt.0.0 )
         cldlayer(:,:) = cldlayer(:,:)/nsublayer(:,:)
      elsewhere
         cldlayer(:,:) = undef
      endwhere

      RETURN
      END SUBROUTINE COSP_CLDFRAC
! ---------------------------------------------------------------
	  
END MODULE MOD_LMD_IPSL_STATS
