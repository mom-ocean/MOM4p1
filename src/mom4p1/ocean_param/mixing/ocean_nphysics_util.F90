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
module ocean_nphysics_util_mod
! 
!<CONTACT EMAIL="stephen.griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Utilities for neutral physics, including the code to compute 
! space-time dependent diffusivities.  
!</OVERVIEW>
!
!<DESCRIPTION>
! Utilities for neutral physics, including the code to compute 
! space-time dependent diffusivities.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! D.B. Chelton,  R.A. deSzoeke, M.G. Schlax, K.E. Naggar, N. Siwertz
! Geographical Variability of the First Baroclinic Rossby Radius of Deformation
! Journal of Physical Oceanography (1998) vol 28 pages 433-460 
! </REFERENCE>
!
! <REFERENCE>
! G. Danabasoglu and J. C. McWilliams
! Sensitivity of the global ocean circulation to 
! parameterizations of mesoscale tracer transports
! Journal of Climate (1995) vol 8 pages 2967--2987 
! </REFERENCE>
!
! <REFERENCE>
! Held and Larichev
! A scaling theory for horizontally homogeneous baroclinically 
! unstable flow on a beta plane
! Journal of Atmospheric Sciences (1996) vol 53 pages 946-952
! </REFERENCE>
!
! <REFERENCE>
! M. Visbeck, J.C. Marshall, T. Haine and M. Spall
! Specification of eddy transfer coefficients in coarse resolution ocean
! circulation models
! Journal of Physical Oceanography (1997) vol 27 pages 381--402
! </REFERENCE>
!
! <REFERENCE>
! D. Ferreira, J. Marshall, and P. Heimbach, 
! Estimating eddy stresses by fitting dynamics to observations 
! using a residual-mean ocean circulation omdel and its adjoint. 
! Journal of Physical Oceanography (2005) vol 35 pages 1891-1910.
! </REFERENCE>
!
! <REFERENCE>
! K. Eden, Eddy length scales in the North Atlantic, 2007,
! Preprint. 
! </REFERENCE>
!
! <REFERENCE>
! K. Eden and R. Greatbatch, 2008: Towards a mesoscale eddy closure,
! Ocean Modelling, vol. 20, pages 223-239
! </REFERENCE>
!
! <NOTE>
! Diffusivities can be determined in a number of manners
!
! TIME INDEPENDENT
!
! Various methods are available for specifying a time 
! independent diffusivity, either globally uniform or 
! with selections of spatial dependence.  
!
! TIME DEPENDENT (as a function of the flow)
!
! Various methods are available for determining the 
! diffusivity that changes in time according to the 
! properties of the fluid.  There are various means 
! for specifying the length and time scales needed
! to set the diffusivity. 
!
! LENGTH SCALES 
!
! 1. First baroclinic Rossby radius (estimated as per Chelton etal).  
! Equatorial Rossby radius is used within 5deg of the equator.
! 
! 2. Width of the baroclinic zone, as done in the Hadley Centre 
! model and documented in the MOM3 Manual.
!
! 3. Specified length scale set independent of the flow. 
!
! TIME SCALE 
!
! When using either of the above for the length scale, 
! the time scale is determined by the Eady growth rate.    
!
! COMBINED LENGTH/TIME SCALE 
! 
! Another option, used in the GFDL CM2.X coupled climate models,
! is to set the diffusivity proportional to the depth averaged
! absolute value of the horizontal density gradient.  
! </NOTE> 
!
! </INFO>
!
!
!<NAMELIST NAME="ocean_nphysics_util_nml">
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  </DATA> 
!
!  <DATA NAME="nphysics_util_zero_init" TYPE="logical">
!  For Time%init=.true. and wishing to ensure starting with a clean 
!  suite of nphysics_util fields, even if ocean_neutral.res.nc exists.  
!  </DATA> 
!
!  <DATA NAME="horz_z_derivative" TYPE="logical">
!  When computing "horizontal" derivatives of tracers, we have the choice of 
!  doing so on surfaces of constant z or surfaces of constant s.  Choosing to 
!  perform derivatives on constant z-surfaces can lead to spurious extrema in
!  in certain regions, such as next to topography, where z and s surfaces deviate
!  a lot.  This problem has been seen in the horz_lap_diffusion module. 
!  It is for this reason that horz_z_diffusion=.false. is the default.
!  </DATA> 
!
!  <DATA NAME="horz_s_derivative" TYPE="logical">
!  This is the preferred setting for computing horizontal derivatives of tracers.
!  This approach ensures that when neutral physics defaults to "horizontal" physics
!  next to boundaries, it will do so as horizontal, defined along surfaces of constant 
!  s-surfaces, and so will not generate spurious extrema.  
!
!  Additionally, when using generalized vertical coordinates, the neutral diffusion
!  slope should be computed relative to the s-surfaces.  The skew diffusion slope 
!  should ideally be computed with respect to z-surfaces, as z-surfaces define
!  available potential energy. However, when s and z surfaces are reasonably close, 
!  as they are in the interior for zstar and pstar vertical coordinates, then we 
!  choose to to dissipate thickness as defined relative to the zstar or pstar surfaces. 
!  This should not be such a big deal, and it is certainly easier computationally than
!  worrying about computing two separate sets of slopes.  More on this detail is 
!  discussed in "Elements of mom4p1".
!  </DATA> 
!
!  <DATA NAME="epsln_drhodz" UNITS="kg/m^3"  TYPE="real">
!  For computing drhodz used in slope calculation.
!  We must keep drhodz < 0 in order to maintain integrity of the 
!  quasi-Stokes streamfunction as well as computation of buoyancy frequency.  
!  Default epsln_drhodz=1e-30.
!  </DATA> 
!  <DATA NAME="drhodz_mom4p1" TYPE="logical">
!  For computing the vertical neutral density deriviative 
!  as in the preferred mom4p1 algorithm rather than the 
!  mom4p0 approach. Default drhodz_mom4p1=.true.
!  </DATA> 
!  <DATA NAME="drhodz_smooth_horz" TYPE="logical">
!  For horizontal laplacian smoothing the vertical derivative 
!  of density prior to its use in computing the neutral slopes. 
!  This smoothing helps to produce regularized slopes.  
!  Note that this option breaks the integrity of the triads
!  and is thus NOT generally recommended.  
!  Default drhodz_smooth_horz=.false.
!  </DATA> 
!  <DATA NAME="drhodz_smooth_vert" TYPE="logical">
!  For vertical 1-2-1 smoothing the vertical derivative of 
!  density prior to its use in computing the neutral slopes.
!  This smoothing helps to produce regularized slopes.  
!  Note that this option breaks the integrity of the triads
!  and is thus NOT generally recommended.  
!  Default drhodz_smooth_vert=.false.
!  </DATA> 
!
!  <DATA NAME="vel_micom_smooth" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in smoothing of drhodzb. 
!  Default vel_micom_smooth=0.2. 
!  </DATA>
!  <DATA NAME="num_121_passes" TYPE="integer">
!  For number of 1-2-1 passes through to smooth drhodz or 
!  eady_rate in vertical. Default num_121_passes=1. 
!  </DATA>
!
!  <DATA NAME="aredi" TYPE="real">
!  Neutral diffusivity used for experiments using a constant diffusivity. 
!  </DATA> 
!  <DATA NAME="agm" TYPE="real">
!  GM-skew diffusivity used for experiments using a constant diffusivity. 
!  </DATA> 
!  <DATA NAME="aredi_equal_agm" TYPE="logical">
!  Will set aredi_array=agm_array, over-riding any other specification 
!  of aredi_array. Default aredi_equal_agm=.true. 
!  </DATA> 
!
!  <DATA NAME="smax" TYPE="real">
!  Value of the maximum neutral direction slope above which the neutral fluxes are 
!  either tapered to zero or saturated.  Typical value is smax=0.01 or smaller. 
!  </DATA> 
!  <DATA NAME="swidth" TYPE="real">
!  Width in slope over which use tanh with dm_taper scheme to taper fluxes in 
!  steep sloped regions. Typical value swidth=0.1*smax
!  </DATA> 
!
!  <DATA NAME="neutral_horz_mix_bdy" TYPE="logical">
!  If .true., then use a horizontal diffusivity in the neutral boundary layer. 
!  </DATA> 
!  <DATA NAME="vel_micom_bdy" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM horizontal diffusivity 
!  within the neutral boundary layer. 
!  </DATA> 
!  <DATA NAME="ah_bdy" UNITS="m^2/sec" TYPE="real">
!  Constant horizontal diffusivity for the boundary layer.  
!  Default ah_bdy=0.0.
!  </DATA> 
!
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the GM-skew diffusivity is set according to a velocity scale 
!  times the grid spacing. 
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!
!  <DATA NAME="agm_lat_zones" TYPE="logical">
!  If true, will set agm_array as constant within two latitudinal zones.  
!  The idea is that one may wish to use a larger agm in the ACC than 
!  elsewhere. 
!  </DATA> 
!  <DATA NAME="agm_lat_zones_boundary" TYPE="real">
!  Boundary between agm in the south and north zones. 
!  </DATA> 
!  <DATA NAME="agm_lat_zones_ratio" TYPE="real">
!  Ratio between the large agm used in the southern latitudinal zone
!  to that used in the north.  
!  agm_array(north) = agm
!  agm_array(south) = agm*agm_lat_zones_ratio
!  </DATA> 
!
!  <DATA NAME="bryan_lewis_aredi" TYPE="logical">
!  Set bryan_lewis_aredi=.true. when want to have aredi a function of depth
!  according to the Bryan and Lewis (1979) profile. Maintained for legacy 
!  purposes, and not recommended for new models.                                 
!  </DATA>
!  <DATA NAME="ahs" TYPE="real">
!  ahs = adjustable parameter at the surface for bryan_lewis_aredi   
!  </DATA>
!  <DATA NAME="ahb" TYPE="real">
!  ahb = adjustable parameter at the bottom for bryan_lewis_aredi   
!  </DATA>
!
!  <DATA NAME="agm_read_restart" TYPE="logical">
!  For those cases with agm_closure=.false. where we wish to read in 
!  the agm_array from restart files and keep the value from the restart.
!  This approach allows us to read in a spatially dependent agm_array 
!  that may have been computed from another integration, but to leave
!  the coefficient static in time.  
!  Default agm_read_restart=.false.
!  </DATA> 
!
!  <DATA NAME="agm_closure" TYPE="logical">
!  If .true. then will compute the GM-skew diffusivity as a function of the flow.
!  The length scale is determined by the Rossby radius and the time scale is 
!  determined by the Eady growth rate.  Diffusivities are depth independent.  
!  </DATA> 
!  <DATA NAME="agm_closure_max" UNITS="m^2/sec" TYPE="real">
!  Maximum GM diffusivity allowed when using agm_closure=.true. 
!  </DATA> 
!  <DATA NAME="agm_closure_min" UNITS="m^2/sec" TYPE="real">
!  Minimum GM diffusivity allowed when using agm_closure=.true. 
!  </DATA> 
!  <DATA NAME="agm_closure_scaling" UNITS="dimensionless" TYPE="logical">
!  Dimensionless tuning parameter for computing flow dependent diffusivities. 
!  </DATA> 
!  <DATA NAME="agm_closure_upper_depth" UNITS="m" TYPE="real">
!  Upper depth where start the depth integration to compute the Eady 
!  growth rate and/or baroclinicity.
!  </DATA> 
!  <DATA NAME="agm_closure_lower_depth" UNITS="m" TYPE="real">
!  Deeper depth where finish the depth integration to compute the Eady 
!  growth rate and/or baroclinicity.
!  </DATA> 
!
!  <DATA NAME="agm_closure_n2_scale" TYPE="logical">
!  For computing the agm coefficient using a 3-dimensional 
!  scaling by (N/Nref)^2, with N=buoyancy frequency and 
!  Nref the buoyancy frequency at the base of the neutral 
!  blayer. Default agm_closure_n2_scale=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_n2_scale_coeff" UNITS="m^2/s" TYPE="real">
!  Coefficient setting the scale for the diffusivity computed from 
!  agm_closure_n2_scale. 
!  Default agm_closure_n2_scale_coeff=1e3. 
!  </DATA> 
!  <DATA NAME="agm_closure_n2_scale_nref_cst" TYPE="logical">
!  For taking the reference buoyancy frequency as agm_closure_buoy_freq
!  for the (N/Nref)^2 scaling.  
!  Default agm_closure_n2_scale_nref_cst=.false.
!  </DATA> 
!
!  <DATA NAME="agm_closure_baroclinic" TYPE="logical">
!  For computing the agm coefficient using only the vertically
!  averaged magnitude of the horizontal density gradient 
!  (i.e., the "baroclinicity").
!  </DATA> 
!  <DATA NAME="agm_closure_buoy_freq" UNITS="sec^-1" TYPE="real">
!  For computing the agm coefficient using only the vertically
!  averaged horizontal density gradient, we need to specify a 
!  buoyancy frequency,  which is taken to be fixed over all space-time.
!  </DATA> 
!
!  <DATA NAME="agm_closure_length_cap" TYPE="logical">
!  For setting a maximum length scale for the agm_closure calculation.
!  Default agm_closure_length_cap=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_length_max" UNITS="metre" TYPE="real">
!  Maximum length scale used for computing agm_closure.  
!  Default agm_closure_length_max=50e3.
!  </DATA> 
!
!  <DATA NAME="agm_closure_length_fixed" TYPE="logical">
!  Use fixed length scale for computing agm_closure diffusivity 
!  </DATA> 
!  <DATA NAME="agm_closure_length" UNITS="meter" TYPE="real">
!  Fixed length scale for use with agm_closure_fixed_length
!  </DATA>
!  
!  <DATA NAME="agm_closure_length_rossby" TYPE="logical">
!  For computing the agm_closure length scale according to Rossby radius. 
!  </DATA> 
!  <DATA NAME="rossby_radius_max" UNITS="meter" TYPE="real">
!  Maximum Rossby radius used for agm_closure_length_rossby and 
!  the neutral_sine_taper. Default = 100e3 m.
!  </DATA> 
!  <DATA NAME="rossby_radius_min" UNITS="meter" TYPE="real">
!  Minimum Rossby Radius used for agm_closure_length_rossby and 
!  the neutral_sine_taper. Default = 15e3 m.
!  </DATA> 
!
!  <DATA NAME="agm_closure_length_bczone" TYPE="logical">
!  For computing the agm_closure length scale according to radius of baroclinic zone. 
!  </DATA> 
!  <DATA NAME="bczone_max_pts" TYPE="integer">
!  Max number of horizontal grid points for use in computing the baroclinic zone radius.  
!  </DATA> 
!  <DATA NAME="agm_closure_bczone_crit_rate" UNITS="sec^-1" TYPE="real">
!  Critical growth rate for determining width of the baroclinic zone. 
!  </DATA> 
!  <DATA NAME="agm_closure_growth_scale" UNITS="dimensionless" TYPE="real">
!  Dimensionless scaling used to set a maximum for agm_growth. 
!  </DATA> 
!
!  <DATA NAME="agm_closure_eden_greatbatch" TYPE="logical">
!  For computing the agm_closure length scale according to minimum 
!  of the Rhines scale and the Rossby radius, and using 3d Eady 
!  growth rate. 
!  </DATA> 
!  <DATA NAME="agm_closure_eden_gamma" TYPE="real" UNITS="dimensionless">
!  For use in regularizing the growth rate used in the eden/greatbatch approach.
!  Default agm_closure_eden_gamma=200. Setting to zero removes the regularization.  
!  </DATA> 
!  <DATA NAME="agm_closure_eden_length_const" TYPE="logical">
!  To set the length scale for agm_closure_eden_greatbatch to constant. 
!  Default agm_closure_eden_length_const=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eden_length" TYPE="real" UNITS="metre">
!  Length scale for use with agm_closure_eden_length_const=.true.
!  Default agm_closure_eden_length=10e3. 
!  </DATA> 
!
!  <DATA NAME="agm_closure_eady_smooth_vert" TYPE="logical">
!  For vertical 1-2-1 smoothing the eady_rate 
!  Default agm_closure_eady_smooth=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eady_smooth_horz" TYPE="logical">
!  For horizontal Laplacian smoothing of eady growth rate. 
!  Default agm_closure_eady_smooth_horz=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eady_cap" TYPE="logical">
!  For capping the eady growth rate to avoid huge values.  
!  Default agm_closure_eady_cap=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_eady_ave_mixed" TYPE="logical">
!  To set the Eady growth rate to its average within mixed layer region.  
!  This is used to avoid spuriously large values which often appear just 
!  in the upper regions of the ocean mixed layer.  
!  Default agm_closure_eady_ave_mixed=.false.
!  </DATA> 
!
!  <DATA NAME="agm_smooth_space" TYPE="logical">
!  For smoothing the agm diffusivity in space when nonconstant diffusivity used. 
!  Default is agm_smooth_space=.false.
!  </DATA> 
!  <DATA NAME="agm_smooth_time" TYPE="logical">
!  For smoothing the agm diffusivity in time when nonconstant diffusivity used. 
!  Default is agm_smooth_time=.false.
!  </DATA> 
!  <DATA NAME="agm_damping_time" UNITS="days" TYPE="real">
!  The damping time used for time smoothing agm_array.
!  Default agm_damping_time=10days.  
!  </DATA> 
!
!  <DATA NAME="agm_closure_grid_scaling" TYPE="logical">
!  For an overall scaling of the agm coefficient, according to
!  the relative resolution of the grid and deformation radius. 
!  Default is agm_closure_grid_scaling=.false.
!  </DATA> 
!  <DATA NAME="agm_closure_grid_scaling_power" TYPE="real">
!  Power used to scale the agm_closure diffusivity. 
!  Default is agm_closure_grid_scaling_power=2.0
!  </DATA> 
!  <DATA NAME="aredi_diffusivity_grid_scaling" TYPE="logical">
!  For an overall scaling of the aredi coefficient, according to
!  the relative resolution of the grid and deformation radius.
!  This option is used only when aredi_equal_agm=.false. 
!  Default is aredi_diffusivity_grid_scaling=.false.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln, pi, grav
use diag_manager_mod, only: register_diag_field, register_static_field, send_data, need_data
use fms_mod,          only: FATAL, WARNING, NOTE, file_exist
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_io_mod,       only: register_restart_field, save_restart, restore_state
use fms_io_mod,       only: restart_file_type
use mpp_mod,          only: mpp_pe, mpp_min, mpp_error, mpp_chksum, stdout, stdlog
use mpp_mod,          only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use mpp_domains_mod,  only: mpp_update_domains, EUPDATE, NUPDATE
use time_manager_mod, only: time_type, increment_time

use ocean_domains_mod,    only: get_local_indices, set_ocean_domain
use ocean_operators_mod,  only: FDX_T, FDY_T, FDX_ZT, FDY_ZT, FMX, FMY, LAP_T, S2D
use ocean_parameters_mod, only: missing_value, onehalf, onefourth
use ocean_parameters_mod, only: rho0, rho0r
use ocean_tracer_diag_mod,only: calc_mixed_layer_depth
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_density_type
use ocean_types_mod,      only: ocean_time_type, ocean_thickness_type, ocean_time_steps_type
use ocean_types_mod,      only: tracer_3d_0_nk_type, tracer_3d_1_nk_type 
use ocean_util_mod,       only: write_timestamp
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk1_2d, wrk2_2d, wrk3_2d


implicit none

public ocean_nphysics_util_init
public ocean_nphysics_coeff_init
public ocean_nphysics_coeff_end
public tracer_derivs 
public neutral_slopes
public compute_eady_rate
public compute_baroclinicity
public compute_rossby_radius
public compute_bczone_radius
public compute_diffusivity
public ocean_nphysics_util_restart
public transport_on_rho_gm
public transport_on_nrho_gm
public transport_on_theta_gm

private 

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: BCzone_domain   

! to write some output to all processors 
integer :: unit=6

! clock ids
integer :: id_clock_tracer_derivs
integer :: id_clock_neutral_slopes
integer :: id_clock_compute_eady_rate
integer :: id_clock_compute_baroclinicity
integer :: id_clock_compute_rossby_radius
integer :: id_clock_compute_bczone_radius
integer :: id_clock_compute_diffusivity 
integer :: id_transport_on_nrho_gm
integer :: id_transport_on_rho_gm
integer :: id_transport_on_theta_gm

! for diagnostics manager 
integer :: id_slope31               =-1
integer :: id_slope32               =-1
integer :: id_aredi                 =-1
integer :: id_agm                   =-1
integer :: id_aredi_3d              =-1
integer :: id_agm_3d                =-1
integer :: id_eady_mld              =-1
integer :: id_eady_rate             =-1
integer :: id_eady_rate_zave        =-1
integer :: id_rossby                =-1
integer :: id_rossby_radius         =-1
integer :: id_rossby_equator        =-1
integer :: id_rossby_nonequator     =-1
integer :: id_bczone                =-1
integer :: id_gw_speed              =-1
integer :: id_sqrt2betaCr           =-1
integer :: id_growth_rate_baroclinic=-1
integer :: id_baroclinicity         =-1
integer :: id_agm_growth_rate       =-1
integer :: id_agm_length            =-1
integer :: id_rhines_length         =-1
integer :: id_agm_qg                =-1
integer :: id_N2slope               =-1
integer :: id_N2_nblayer_base       =-1
integer :: id_N2_for_agm            =-1
integer :: id_coriolis_param        =-1
integer :: id_beta_param            =-1
integer :: id_grid_length           =-1
integer :: id_agm_grid_scaling      =-1
integer :: id_aredi_grid_scaling    =-1
integer :: id_ksurf_blayer          =-1


! for transport_on_nrho_gm
integer :: id_tx_trans_nrho_gm=-1
integer :: id_ty_trans_nrho_gm=-1

! for transport_on_rho_gm
integer :: id_tx_trans_rho_gm=-1
integer :: id_ty_trans_rho_gm=-1

! for transport_on_theta_gm
integer :: id_tx_trans_theta_gm=-1
integer :: id_ty_trans_theta_gm=-1

! for restart
type(restart_file_type), save :: Nphysics_util_restart

integer :: num_prog_tracers=0
integer :: index_temp=-1
integer :: index_salt=-1

logical :: used

! 2D array of micom gm diffusivities (m^2/sec)
real, dimension(:,:), allocatable :: agm_micom 

! 2D array of deformation radius scaling
real, dimension(:,:), allocatable :: agm_grid_scaling
real, dimension(:,:), allocatable :: aredi_grid_scaling

! 3D array of gm length scale (metre)
real, dimension(:,:,:), allocatable :: agm_length 

! 3D array of agm_array for local purposes (m^2/sec)
real, dimension(:,:,:), allocatable :: agm_array_local

! 3D array of gm growth rate scale (sec^-1)
real, dimension(:,:,:), allocatable :: agm_growth_rate

! 2D array of gm growth rate max (sec^-1)
real, dimension(:,:), allocatable :: agm_growth_rate_max

! grid length array (metre)
real, dimension(:,:),  allocatable :: grid_length

! absolute value of the Coriolis parameter (sec^-1)
real, dimension(:,:),  allocatable :: coriolis_param 

! beta = d(Coriolis)/dy (m^-1 sec^-1)
real, dimension(:,:),  allocatable :: beta_param  

! for determining baroclinic zone radius (metre) 
real, dimension(:,:), allocatable :: bczone_dxt  
real, dimension(:,:), allocatable :: bczone_dyt  
real, dimension(:,:), allocatable :: bczone_rate 

! for Eady growth rate and baroclinicity calculation 
real, dimension(:,:),  allocatable :: count_x 
real, dimension(:,:),  allocatable :: count_y 

! Eady growth rate (sec^-1)
real, dimension(:,:,:),  allocatable :: eady_rate
real, dimension(:,:),    allocatable :: eady_rate_zave

! 1st baroclinic gravity wave speed (m/s)
real, dimension(:,:),  allocatable :: gravity_wave_speed

! for Eden and Greatbatch near equator (m)
real, dimension(:,:),  allocatable :: sqrt2betaCr

!vertically ave abs(horiz density gradient) (kg/m^4)
real, dimension(:,:),  allocatable :: baroclinicity

! (m^2/sec) for smoothing
real, dimension(:,:),  allocatable :: smooth_lap

! some useful constants 
real :: grav_r
real :: gravrho0r
real :: gravrho0r_buoyr
real :: gamma_damp  
real :: dtime


!*************************************
! nml settings 

! for computing horizontal dervatives on constant z-surfaces or constant s-surfaces
logical :: horz_z_derivative = .false. 
logical :: horz_s_derivative = .true.  

! for how to compute drhodz prior to computing slopes 
logical :: drhodz_mom4p1      =.true.
logical :: drhodz_smooth_horz =.false.
logical :: drhodz_smooth_vert =.false.
integer :: num_121_passes     = 1
real    :: vel_micom_smooth   = 0.2   ! m/sec
real    :: epsln_drhodz       = 1e-30    

! globally constant diffusivities 
real    :: aredi = 1.0e3 ! constant neutral diffusion tracer diffusivity (m^2/sec)
real    :: agm   = 1.0e3 ! constant gent-mcwilliams skew-diffusion diffusivity (m^2/sec)

! will set aredi_array=agm_array 
logical :: aredi_equal_agm =.true.

! maximum neutral direction slope allowed before tapering fluxes
real    :: smax=0.01

! slope width over which fluxes tapered using tanh function using dm_taper scheme 
real    :: swidth = 0.05*.01

! if true, apply horizontal diffusivity in neutral boundary layer for nphysicsA module
logical :: neutral_horz_mix_bdy =.false. 
real    :: vel_micom_bdy        = 0.0  ! velocity scale (m/s) for horizontal diffusivity w/i boundary 
real    :: ah_bdy               = 0.0  ! constant horiz diffusivity in neutral bdy

! time independent, vertically dependent neutral diffusivity according to Bryan-Lewis. 
! maintained in mom4 only for legacy purposes.  not recommended for new models. 
logical :: bryan_lewis_aredi = .false.  
real    :: ahs               = 0.0      
real    :: ahb               = 0.0      

! for setting agm according to latitude zones 
logical :: agm_lat_bands          = .false. 
real    :: agm_lat_bands_boundary = -999.   ! boundary between agm in the south and north zones
real    :: agm_lat_bands_ratio    = 1.0     ! ratio agm(south)/agm(north)

! set diffusivity according to local horizontal area 
! of grid and a specified velocity scale (m/s)
logical :: tracer_mix_micom =.false. 
real    :: vel_micom        = 0.0    

! determining how to read in restart diffusivity
logical :: agm_read_restart=.false.

! for setting diffusivity according to flow properties 
logical :: agm_closure         =.false. 
real    :: agm_closure_scaling = 2.0    ! dimensionless tuning parameter 
real    :: agm_closure_max     = 2.e3   ! maximum diffusivity allowed when agm_closure=.true.
real    :: agm_closure_min     = 2.e2   ! minimum diffusivity allowed when agm_closure=.true.

! for fixed length scale (metres) set by agm_closure_length
logical :: agm_closure_length_fixed =.false. 
real    :: agm_closure_length       = 50.e3  ! metres

! for length scale set according to estimate of first baroclinic Rossby radius
logical :: agm_closure_length_rossby =.false. 
real    :: rossby_radius_max=100e3  ! metres 
real    :: rossby_radius_min=15e3   ! metres 

! for length scale set according to radius of baroclinic zone
logical :: agm_closure_length_bczone    =.false. 
integer :: bczone_max_pts               =10      ! max # points searched for determining baroclinic zone width 
real    :: agm_closure_bczone_crit_rate =1.4e-6  ! critical growth rate for determining baroclinic zone (sec^-1)

! for capping the agm_length scale 
logical :: agm_closure_length_cap = .false.
real    :: agm_closure_length_max = 50e3

! for scaling diffusivity by (N/Nref)^2
logical :: agm_closure_n2_scale=.false. 
logical :: agm_closure_n2_scale_nref_cst=.false.
real    :: agm_closure_n2_scale_coeff=1e3

! for diffusity according to Eden and Greatbatch (2008) and Eden(2007)
logical :: agm_closure_eden_greatbatch  =.false. 
logical :: agm_closure_eden_length_const=.false.
real    :: agm_closure_eden_gamma  = 200.0
real    :: agm_closure_eden_length = 10e3

! for diffusivity set proportional to vertically averaged horiz density gradient  
logical :: agm_closure_baroclinic = .true. 
! sec^-1 buoyancy frequency for use with agm_closure_baroclinic
real    :: agm_closure_buoy_freq  = 0.004  

! depths over which compute the Eady growth and/or horizontal density gradient 
real    :: agm_closure_upper_depth  =100.0  ! metre
real    :: agm_closure_lower_depth  =2000.0 ! metre 

! for smoothing and capping eady growth rate
logical :: agm_closure_eady_smooth_vert=.false.
logical :: agm_closure_eady_smooth_horz=.false.
logical :: agm_closure_eady_cap        =.false.
logical :: agm_closure_eady_ave_mixed  =.false.

! dimensionless number to set maximum value of agm_growth. 
! agm_closure_growth_scale =0.5 yields agm_growth max=coriolis parameter.
real    :: agm_closure_growth_scale =0.5    

! for smoothing agm_array in space and time
logical :: agm_smooth_space=.false.
logical :: agm_smooth_time =.false.
real    :: agm_damping_time=10.0    

! for overall scaling of the agm coefficient 
logical :: agm_closure_grid_scaling=.false.
real    :: agm_closure_grid_scaling_power=2.0

! for scaling aredi_array according to size of grid and Rossby radius  
logical :: aredi_diffusivity_grid_scaling=.false.

! for Tim%init=.true. and ocean_neutral.res.nc exists, but still wish
! to start from initial fields.
logical :: nphysics_util_zero_init=.false.

! for debugging 
logical :: debug_this_module = .false.

!*************************************

character(len=128) :: version=&
     '$Id: ocean_nphysics_util.F90,v 1.1.2.28.14.5.18.2.18.2.30.1.4.1 2009/10/21 20:07:28 smg Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: module_is_initialized = .FALSE.

namelist /ocean_nphysics_util_nml/ debug_this_module, nphysics_util_zero_init,      &
          horz_z_derivative, horz_s_derivative, smax, swidth,                       &
          epsln_drhodz, drhodz_mom4p1,                                              &
          drhodz_smooth_horz, drhodz_smooth_vert, num_121_passes,                   &
          aredi, agm, aredi_equal_agm, tracer_mix_micom, vel_micom,                 &
          bryan_lewis_aredi, ahs, ahb,                                              &
          neutral_horz_mix_bdy, vel_micom_bdy, ah_bdy,                              &
          agm_lat_bands, agm_lat_bands_boundary, agm_lat_bands_ratio,               &  
          rossby_radius_max, rossby_radius_min,  agm_read_restart,                  &
          agm_closure, agm_closure_scaling, agm_closure_max, agm_closure_min,       &
          agm_closure_growth_scale, agm_closure_length_fixed, agm_closure_length,   &
          agm_closure_length_rossby, agm_closure_length_bczone,                     &
          bczone_max_pts, agm_closure_bczone_crit_rate,                             &
          agm_closure_eden_greatbatch, agm_closure_eden_gamma,                      &
          agm_closure_eden_length_const, agm_closure_eden_length,                   &
          agm_closure_eady_smooth_vert, agm_closure_eady_smooth_horz,               &
          agm_closure_eady_ave_mixed, agm_closure_eady_cap,                         &
          agm_closure_baroclinic, agm_closure_buoy_freq,                            &   
          agm_closure_upper_depth, agm_closure_lower_depth,                         &
          agm_closure_length_cap, agm_closure_length_max,                           &
          agm_smooth_space, vel_micom_smooth, agm_smooth_time, agm_damping_time,    &
          agm_closure_grid_scaling, agm_closure_grid_scaling_power,                 &
          aredi_diffusivity_grid_scaling,                                           &
          agm_closure_n2_scale, agm_closure_n2_scale_coeff,                         &
          agm_closure_n2_scale_nref_cst
 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_util_init">
!
! <DESCRIPTION>
! Initialize the utility module for neutral physics.
! </DESCRIPTION>
!
subroutine ocean_nphysics_util_init(Grid, Domain, Time, Time_steps, Dens, T_prog, &
           agm_closure_lower_dept, agm_closure_upper_dept, agm_closure_buoy_frq,  &
           smx, swidt, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  real,                         intent(inout)        :: agm_closure_lower_dept
  real,                         intent(inout)        :: agm_closure_upper_dept
  real,                         intent(inout)        :: agm_closure_buoy_frq
  real,                         intent(inout)        :: smx 
  real,                         intent(inout)        :: swidt
  logical,                      intent(in), optional :: debug

  integer :: ioun, io_status, ierr  
  integer :: i,j,n
  integer :: num_schemes 
  real    :: param 
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysics_util_mod: already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  Dom => Domain
  Grd => Grid

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_nphysics_util_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_nphysics_util_nml)  
  write (stdlogunit,ocean_nphysics_util_nml)
  ierr = check_nml_error(io_status,'ocean_nphysics_util_nml')
  call close_file (ioun)

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_nphysics_util_mod with debug_this_module=.true.'  
  endif

#ifndef MOM4_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  ! for setting up temp and salt tracers
  num_prog_tracers = size(T_prog(:))
  index_temp=-1;index_salt=-1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp == -1 .or. index_salt == -1) then 
     call mpp_error(FATAL, &
     '==>Error: temp and/or salt not identified in call to ocean_nphysics_util_init')
  endif 

  ! some useful constants 
  dtime           = Time_steps%dtime_t 
  grav_r          = 1.0/grav 
  gravrho0r       = grav*rho0r                       !for buoyancy frequency calculation  
  gravrho0r_buoyr = gravrho0r/agm_closure_buoy_freq  !for agm_closure_baroclinic
  gamma_damp      = dtime/(86400.0*agm_damping_time) !for damping time dependent agm_array 

  ! for sending back to ocean_nphysics_mod 
  agm_closure_lower_dept = agm_closure_lower_depth
  agm_closure_upper_dept = agm_closure_upper_depth
  agm_closure_buoy_frq   = agm_closure_buoy_freq
  smx                    = smax
  swidt                  = swidth 

  ! Coriolis parameter and beta parameter 
  allocate (coriolis_param(isd:ied,jsd:jed))
  allocate (beta_param(isd:ied,jsd:jed))
  coriolis_param(:,:) = 2.0*7.292e-5*abs(sin(Grd%phit(:,:)))
  beta_param(:,:)     = 2.28e-11*abs(cos(Grd%phit(:,:)))

  allocate (gravity_wave_speed(isd:ied,jsd:jed))
  allocate (sqrt2betaCr(isd:ied,jsd:jed))
  allocate (eady_rate(isd:ied,jsd:jed,nk))
  allocate (eady_rate_zave(isd:ied,jsd:jed))
  allocate (baroclinicity(isd:ied,jsd:jed))
  sqrt2betaCr(:,:)        = 0.0
  gravity_wave_speed(:,:) = 0.0
  eady_rate(:,:,:)        = 0.0
  eady_rate_zave(:,:)     = 0.0
  baroclinicity(:,:)      = 0.0

  ! for computing the baroclinic zone radius using Hadley Centre search algorithm 
  call set_ocean_domain(BCzone_domain,Grd,xhalo=bczone_max_pts,yhalo=bczone_max_pts,name='bczone',maskmap=Dom%maskmap)
  allocate (bczone_rate(isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  allocate (bczone_dxt (isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  allocate (bczone_dyt (isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  bczone_dxt=0.0 ; bczone_dyt=0.0 ; bczone_rate=0.0
  bczone_dxt(isc:iec,jsc:jec) = Grd%dxt(isc:iec,jsc:jec)
  bczone_dyt(isc:iec,jsc:jec) = Grd%dyt(isc:iec,jsc:jec)
  call mpp_update_domains (bczone_dxt(:,:), BCzone_domain%domain2d)
  call mpp_update_domains (bczone_dyt(:,:), BCzone_domain%domain2d)

  allocate (grid_length(isd:ied,jsd:jed))
  grid_length(:,:) = 2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))

  allocate (agm_micom(isd:ied,jsd:jed))
  if(tracer_mix_micom) then 
       agm_micom(:,:) = vel_micom*grid_length(:,:)
  else 
       agm_micom(:,:) = 0.0
  endif 

  allocate (agm_length(isd:ied,jsd:jed,nk))
  agm_length(:,:,:) = agm_closure_length*Grd%tmask(:,:,:)

  allocate (agm_grid_scaling(isd:ied,jsd:jed))
  agm_grid_scaling(:,:) = Grd%tmask(:,:,1)

  allocate (aredi_grid_scaling(isd:ied,jsd:jed))
  aredi_grid_scaling(:,:) = Grd%tmask(:,:,1)

  allocate (agm_growth_rate(isd:ied,jsd:jed,nk))
  agm_growth_rate(:,:,:) = 0.0

  allocate (agm_array_local(isd:ied,jsd:jed,nk))
  agm_array_local(:,:,:) = agm*Grd%tmask(:,:,:)


  !max agm_growth_rate is 2.0*agm_closure_growth_scale*param
  !if agm_closure_growth_scale=0.5, then agm_growth_rate is <= param
  allocate (agm_growth_rate_max(isd:ied,jsd:jed))
  do j=jsd,jed
     do i=isd,ied
        param = max(coriolis_param(i,j), sqrt2betaCr(i,j))
        agm_growth_rate_max(i,j) = agm_closure_growth_scale*param
     enddo
  enddo

  if(Time_steps%aidif /= 1.0 .and. aredi /= 0) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_nphysics_util_mod: stability requires aidif=1.0 to handle K33 implicitly.')
  endif

  if(agm_closure) then

      num_schemes = 0 
      if(agm_closure_length_fixed) then 
          num_schemes = num_schemes + 1
          write(stdoutunit,'(1x,a)') &
          '==>Note: Computing 2d flow-dependent tracer diffusivity with agm_closure_length_fixed.'
      endif
      if(agm_closure_length_rossby) then 
          num_schemes = num_schemes + 1
          write(stdoutunit,'(1x,a)') &
          '==>Note: Computing 2d flow-dependent tracer diffusivity with agm_closure_length_rossby.'
      endif
      if(agm_closure_length_bczone) then 
          num_schemes = num_schemes + 1
          write(stdoutunit,'(1x,a)') &
          '==>Note: Computing 2d flow-dependent tracer diffusivity with agm_closure_length_bczone.'
      endif
      if(agm_closure_eden_greatbatch) then 
          num_schemes = num_schemes + 1
          write(stdoutunit,'(1x,a)') &
          '==>Note: Computing 3d flow-dependent tracer diffusivity with agm_closure_eden_greatbatch.'
      endif
      if(agm_closure_baroclinic) then 
          num_schemes = num_schemes + 1
          write(stdoutunit,'(1x,a)') &
          '==>Note: Computing 2d flow-dependent tracer diffusivity with agm_closure_baroclinic.'
      endif
      if(agm_closure_n2_scale) then 
          num_schemes = num_schemes + 1
          write(stdoutunit,'(1x,a)') &
          '==>Note: Computing 3d flow-dependent tracer diffusivity with agm_closure_n2_scale.'
      endif

      if(num_schemes == 0) then 
          call mpp_error(FATAL, &
          '==>Error: with agm_closure=.true., must choose one of the "agm_closure" methods')
      endif
      if(num_schemes > 1) then 
          call mpp_error(FATAL, &
          '==>Error: with agm_closure=.true., can choose only one of the available agm_closure methods')
      endif

      write(stdoutunit,'(a,e10.5)') &
      '        The maximum allowable diffusivity (m^2/s) is given by ',agm_closure_max
      write(stdoutunit,'(a,e10.5)') &
      '        The minimum allowable diffusivity (m^2/s) is given by ',agm_closure_min

      write(stdoutunit,'(a,e10.5,a,e10.5)') &
      '  Depths (m) between which compute eady growth and baroclinicity = ', &
         agm_closure_upper_depth, ' ',agm_closure_lower_depth


      if(agm_closure_n2_scale) then 
          write(stdoutunit,'(a)')  '  Diffusivity will be computed as agm = coeff*(N/Nref)^2.'
      endif 

      if(agm_closure_baroclinic) then 
          write(stdoutunit,'(a)')       &
               '        Length and time scales set by vertically averaged baroclinicity |grad(rho)|,'
          write(stdoutunit,'(a,e10.5)') &
               '        as well as the constant buoyancy freq(sec^-1) = ',agm_closure_buoy_freq
          write(stdoutunit,'(a,e10.5)') &
               '        and the constant length scale (m) = ',agm_closure_length
      else 
          write(stdoutunit,'(a)') &
               '        Eady growth rate gives inverse time scale.'

          if(agm_closure_length_fixed) then 
              write(stdoutunit,'(a,e12.5,a)') &
                   '        Length scale set by nml parameter agm_closure_length =',agm_closure_length,' metre'
          elseif(agm_closure_length_rossby) then 
              write(stdoutunit,'(a)') &
                   '        First baroclinic Rossby radius gives the length scale.'
          elseif(agm_closure_length_bczone) then 
              write(stdoutunit,'(a)') &
                   '        Radius of baroclinic zone gives the length scale.'
          elseif(agm_closure_eden_greatbatch) then
              if(agm_closure_eden_length_const) then 
                  write(stdoutunit,'(a,e12.5,a)') &
                       '        Take constant length scale of ',agm_closure_eden_length, 'metre'
              else 
                  write(stdoutunit,'(a)') &
                       '        Min(rossby,rhines) gives the length scale.'
              endif
          endif
      endif

  endif   ! endif for agm_closure 


  ! for smoothing
  allocate (smooth_lap(isd:ied,jsd:jed))
  smooth_lap(:,:) = vel_micom_smooth*grid_length(:,:)

  ! for Eady growth rate and baroclinicity calculation 
  allocate (count_x(isd:ied,jsd:jed))
  allocate (count_y(isd:ied,jsd:jed))
  count_x(:,:) = 0.0
  count_y(:,:) = 0.0
  do i=isc,iec
    do j=jsc,jec
      count_x(i,j) =  &
      min(Grd%tmask(i+1,j,1) + Grd%tmask(i-1,j,1), 1.0/(2.0*(Grd%tmask(i+1,j,1) + Grd%tmask(i-1,j,1) + epsln)))
      count_y(i,j) =  &
      min(Grd%tmask(i,j+1,1) + Grd%tmask(i,j-1,1), 1.0/(2.0*(Grd%tmask(i,j+1,1) + Grd%tmask(i,j-1,1) + epsln)))
    enddo
  enddo

  ! initialize clock ids 
  id_clock_tracer_derivs         = mpp_clock_id('(Ocean neutral: tracer derivs)'  ,grain=CLOCK_ROUTINE)
  id_clock_neutral_slopes        = mpp_clock_id('(Ocean neutral: slopes)'         ,grain=CLOCK_ROUTINE)
  id_clock_compute_eady_rate     = mpp_clock_id('(Ocean neutral: eady rate)'      ,grain=CLOCK_ROUTINE)
  id_clock_compute_baroclinicity = mpp_clock_id('(Ocean neutral: baroclinic)'     ,grain=CLOCK_ROUTINE)
  id_clock_compute_rossby_radius = mpp_clock_id('(Ocean neutral: rossby radius)'  ,grain=CLOCK_ROUTINE)
  id_clock_compute_bczone_radius = mpp_clock_id('(Ocean neutral: bc zone radius)' ,grain=CLOCK_ROUTINE)
  id_clock_compute_diffusivity   = mpp_clock_id('(Ocean neutral: diffusivity)'    ,grain=CLOCK_ROUTINE)
  id_transport_on_nrho_gm        = mpp_clock_id('(Ocean neutral: nrho-trans)'     ,grain=CLOCK_ROUTINE)
  id_transport_on_rho_gm         = mpp_clock_id('(Ocean neutral: rho-trans)'      ,grain=CLOCK_ROUTINE)
  id_transport_on_theta_gm       = mpp_clock_id('(Ocean neutral: theta-trans)'    ,grain=CLOCK_ROUTINE)


  ! for diagnostics manager 

  id_ksurf_blayer= -1            
  id_ksurf_blayer= register_diag_field ('ocean_model', 'ksurf_blayer',     &
                Grd%tracer_axes(1:2), Time%model_time,                     &
                'k-value at base of surf nblayer region', 'dimensionless', &
                missing_value=missing_value, range=(/-1.0,1.e1/))

  id_N2slope = -1
  id_N2slope = register_diag_field ('ocean_model', 'N2slope',               & 
            Grd%tracer_axes_wt(1:3), Time%model_time,                       &
            'Squared buoyancy frequency used in neutral slope calculation', &
            '1/s^2', missing_value=missing_value, range=(/-1e3,1.e10/))

  id_N2_for_agm = -1
  id_N2_for_agm = register_diag_field ('ocean_model', 'N2_for_agm',           & 
            Grd%tracer_axes(1:3), Time%model_time,                            &
            'Squared buoyancy frequency used for (N/Nref)^2 agm calculation', &
            '1/s^2', missing_value=missing_value, range=(/-1e3,1.e10/))

  id_N2_nblayer_base= -1
  id_N2_nblayer_base= register_diag_field ('ocean_model', 'N2_nblayer_base',& 
            Grd%tracer_axes(1:2), Time%model_time,                          &
            'Squared buoyancy frequency at base of nblayer',                &
            '1/s^2', missing_value=missing_value, range=(/-1e3,1.e10/))

  id_slope31 = -1
  id_slope31 = register_diag_field ('ocean_model', 'slope31',   &
               Grd%tracer_axes(1:3), Time%model_time,           &
               'neutral slope -(rho_x/rho_z)', 'dimensionless', &
               missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_slope32 = -1
  id_slope32 = register_diag_field ('ocean_model', 'slope32',   &
               Grd%tracer_axes(1:3), Time%model_time,           &
               'neutral slope -(rho_y/rho_z)', 'dimensionless', &
               missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_eady_mld = -1          
  id_eady_mld = register_diag_field ('ocean_model', 'eady_mld',        &
            Grd%tracer_axes(1:2), Time%model_time,                     &
            'Mixed layer depth inside of which compute ave Eady rate', &
            's^-1', missing_value=missing_value, range=(/-10.0,1.e10/))

  id_eady_rate_zave = -1          
  id_eady_rate_zave = register_diag_field ('ocean_model', 'eady_rate_zave',&
            Grd%tracer_axes(1:2), Time%model_time,                         &
            'Eady growth rate averaged over depth range for nphysics ',    &
            's^-1', missing_value=missing_value, range=(/-10.0,1.e10/))

  id_eady_rate = -1          
  id_eady_rate = register_diag_field ('ocean_model', 'eady_rate',&
            Grd%tracer_axes(1:3), Time%model_time,               &
            'Eady growth rate used in neutral physics ', 's^-1', &
            missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby = -1      
  id_rossby = register_diag_field ('ocean_model', 'rossby', &
              Grd%tracer_axes(1:2), Time%model_time,        &
              'Rossby radius used in neutral physics', 'm', &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby_radius = -1      
  id_rossby_radius = register_diag_field ('ocean_model', 'rossby_radius', &
              Grd%tracer_axes(1:2), Time%model_time,                      &
              'Rossby radius computed without min/max bounds', 'm',       &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby_equator = -1      
  id_rossby_equator = register_diag_field ('ocean_model', 'rossby_equator', &
              Grd%tracer_axes(1:2), Time%model_time,                        &
              'Equatorial Rossby radius', 'm',                              &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rossby_nonequator = -1      
  id_rossby_nonequator = register_diag_field ('ocean_model', 'rossby_nonequator', &
              Grd%tracer_axes(1:2), Time%model_time,                              &
              'Rossby radius outside equatorial region', 'm',                     &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_bczone = -1      
  id_bczone = register_diag_field ('ocean_model', 'bczone', &
              Grd%tracer_axes(1:2), Time%model_time,        &
              'radius of baroclinic zone', 'm',             &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_sqrt2betaCr = -1          
  id_sqrt2betaCr = register_diag_field ('ocean_model', 'sqrt2betaCr',         &
            Grd%tracer_axes(1:2), Time%model_time,                            &
            'sqrt(2 * 1st b/c wave speed *beta) in neutral physics ', '1/sec',&
            missing_value=missing_value, range=(/-1e0,1e10/))

  id_gw_speed = -1          
  id_gw_speed = register_diag_field ('ocean_model', 'gw_speed',  &
            Grd%tracer_axes(1:2), Time%model_time,               &
            'Gravity wave speed used in neutral physics ', 'm/s',&
            missing_value=missing_value, range=(/-1e2,1e2/))

  id_baroclinicity = -1          
  id_baroclinicity = register_diag_field ('ocean_model', 'baroclinicity',   &
                     Grd%tracer_axes(1:2), Time%model_time,                 &
                     'vertically averaged horz density gradient', 'kg/m^4', &
                     missing_value=missing_value, range=(/-1e10,1.e10/))

  id_growth_rate_baroclinic = -1          
  id_growth_rate_baroclinic = register_diag_field ('ocean_model', 'growth_rate_baroclinic', &
                              Grd%tracer_axes(1:3), Time%model_time,                        &
                              'growth rate using baroclinicity', 's^-1',                    &
                              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm_growth_rate = -1          
  id_agm_growth_rate = register_diag_field ('ocean_model', 'agm_growth_rate',&
                  Grd%tracer_axes(1:3), Time%model_time,                     &
                  'effective growth rate for agm', 's^-1',                   &
                  missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm_length = -1          
  id_agm_length = register_diag_field ('ocean_model', 'agm_length',  &
                  Grd%tracer_axes(1:3), Time%model_time,             &
                  'effective length scale for agm', 'm',             &
                  missing_value=missing_value, range=(/-10.0,1.e10/))

  id_rhines_length = -1          
  id_rhines_length = register_diag_field ('ocean_model', 'rhines_length', &
                  Grd%tracer_axes(1:3), Time%model_time,                  &
                  'Rhines length approximated as eady_rate/beta', 'm',    &
                  missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm_qg = -1          
  id_agm_qg = register_diag_field ('ocean_model', 'agm_qg',  &
              Grd%tracer_axes(1:3), Time%model_time,         &
              'agm from QG theory', 'm^2/s',                 &
              missing_value=missing_value, range=(/-10.0,1.e10/))

  id_agm = -1
  id_agm = register_diag_field ('ocean_model', 'agm',         &
           Grd%tracer_axes(1:2), Time%model_time,             &
           'GM diffusivity at surface', 'm^2/sec',            &
           missing_value=missing_value, range=(/-10.0,1.e10/),&
           standard_name='ocean_tracer_bolus_laplacian_diffusivity')

  id_agm_grid_scaling = register_diag_field ('ocean_model','agm_grid_scaling',          &
                        Grd%tracer_axes(1:2), Time%model_time,                          &
                       'Scaling of AGM according to Delta(s)^2/(Delta(s)^2 + Rossby^2)',&
                        'dimensionless', missing_value=missing_value, range=(/-1e1,1e1/))

  id_aredi_grid_scaling = register_diag_field ('ocean_model','aredi_grid_scaling',        &
                        Grd%tracer_axes(1:2), Time%model_time,                            &
                       'Scaling of Aredi according to Delta(s)^2/(Delta(s)^2 + Rossby^2)',&
                        'dimensionless', missing_value=missing_value, range=(/-1e1,1e1/))

  id_agm_3d = -1
  id_agm_3d = register_diag_field ('ocean_model', 'agm_3d',      &
              Grd%tracer_axes(1:3), Time%model_time,             &
              '3d GM diffusivity', 'm^2/sec',                    &
              missing_value=missing_value, range=(/-10.0,1.e10/),&
              standard_name='ocean_tracer_bolus_laplacian_diffusivity')

  id_aredi = -1            
  id_aredi = register_diag_field ('ocean_model', 'aredi',       &
             Grd%tracer_axes(1:2), Time%model_time,             &
             'neutral diffusivity at k=1', 'm^2/sec',           &
             missing_value=missing_value, range=(/-10.0,1.e20/),&
             standard_name='ocean_tracer_epineutral_laplacian_diffusivity')

  id_aredi_3d = -1            
  id_aredi_3d = register_diag_field ('ocean_model', 'aredi_3d', &
                Grd%tracer_axes(1:3), Time%model_time,          &
                '3d neutral diffusivity', 'm^2/sec',            &
                missing_value=missing_value, range=(/-10.0,1.e20/))

  id_tx_trans_nrho_gm = register_diag_field ('ocean_model','tx_trans_nrho_gm', Dens%neutralrho_axes(1:3),  &
                        Time%model_time, 'T-cell i-mass transport from GM on neutral rho','Sv (10^9 kg/s)',&
                        missing_value=missing_value, range=(/-1e9,1e9/))
  id_ty_trans_nrho_gm = register_diag_field ('ocean_model','ty_trans_nrho_gm', Dens%neutralrho_axes(1:3),  &
                        Time%model_time, 'T-cell j-mass transport from GM on neutral rho','Sv (10^9 kg/s)',&
                        missing_value=missing_value, range=(/-1e9,1e9/))

  id_tx_trans_rho_gm = register_diag_field ('ocean_model','tx_trans_rho_gm', Dens%potrho_axes(1:3),   &
                       Time%model_time, 'T-cell i-mass transport from GM on pot_rho','Sv (10^9 kg/s)',&
                       missing_value=missing_value, range=(/-1e9,1e9/))
  id_ty_trans_rho_gm = register_diag_field ('ocean_model','ty_trans_rho_gm', Dens%potrho_axes(1:3),   &
                       Time%model_time, 'T-cell j-mass transport from GM on pot_rho','Sv (10^9 kg/s)',&
                       missing_value=missing_value, range=(/-1e9,1e9/))

  id_tx_trans_theta_gm = register_diag_field ('ocean_model','tx_trans_theta_gm', Dens%theta_axes(1:3),&
                       Time%model_time, 'T-cell i-mass transport from GM on theta','Sv (10^9 kg/s)',  &
                       missing_value=missing_value, range=(/-1e9,1e9/))
  id_ty_trans_theta_gm = register_diag_field ('ocean_model','ty_trans_theta_gm', Dens%theta_axes(1:3),&
                       Time%model_time, 'T-cell j-mass transport from GM on theta','Sv (10^9 kg/s)',  &
                       missing_value=missing_value, range=(/-1e9,1e9/))

  ! static fields 
  id_coriolis_param = register_static_field ('ocean_model', 'coriolis_param', Grd%tracer_axes(1:2), &
                                        'Coriolis frequency on T-cell for nphysics', '1/s',         &
                                         missing_value=missing_value, range=(/-10.0,10.0/))
  if (id_coriolis_param > 0) used = send_data (id_coriolis_param, coriolis_param(:,:), &
                                    Time%model_time, rmask=Grd%tmask(:,:,1),           &
                                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  id_beta_param = register_static_field ('ocean_model', 'beta_param', Grd%tracer_axes(1:2),&
                                        'Beta=df/dy on T-cell for nphysics', '1/s',        &
                                         missing_value=missing_value, range=(/-10.0,10.0/))
  if (id_beta_param > 0) used = send_data (id_beta_param, beta_param(:,:),   &
                                    Time%model_time, rmask=Grd%tmask(:,:,1), &
                                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  id_grid_length = register_static_field ('ocean_model', 'grid_length', Grd%tracer_axes(1:2),   &
                                        'Grid length scale used for nphysics calculations', 'm',&
                                         missing_value=missing_value, range=(/-10.0,1e10/))
  if (id_grid_length > 0) used = send_data (id_grid_length, grid_length(:,:), &
                                    Time%model_time, rmask=Grd%tmask(:,:,1),  &
                                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)



end subroutine ocean_nphysics_util_init
! </SUBROUTINE>  NAME="ocean_nphysics_util_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_coeff_init">
!
! <DESCRIPTION>
! Initialize the diffusivities used in neutral physics.
! Need to initialize them after the ocean_nphysics_util_init routine,
! since need to have the domain parameters known for passing the 
! array size information into the ocean_nphysics_coeff_init routine. 
! </DESCRIPTION>
!

subroutine ocean_nphysics_coeff_init(Time, Thickness, rossby_radius, rossby_radius_raw, &
                      bczone_radius, agm_array, aredi_array, ah_array)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  real, dimension(isd:,jsd:),   intent(inout) :: rossby_radius
  real, dimension(isd:,jsd:),   intent(inout) :: rossby_radius_raw
  real, dimension(isd:,jsd:),   intent(inout) :: bczone_radius
  real, dimension(isd:,jsd:,:), intent(inout) :: agm_array
  real, dimension(isd:,jsd:,:), intent(inout) :: aredi_array
  real, dimension(isd:,jsd:),   intent(inout) :: ah_array

  integer :: i,j,k
  integer :: i_delta, j_delta, k_delta
  integer :: id_restart
  real    :: ft, deltaX, deltaY, delta, delta_iso, delta_iso0, delta_min, A_max, A_max0
  character(len=64) :: file_name

  integer :: stdoutunit 
  stdoutunit=stdout() 

  ! initialization based on constant coefficients 
  agm_array_local(:,:,:) = agm*Grd%tmask(:,:,:)
  agm_array(:,:,:)       = agm*Grd%tmask(:,:,:)
  aredi_array(:,:,:)     = aredi*Grd%tmask(:,:,:)

  ! horizontal diffusivity in neutral boundary layer region 
  ah_array(:,:) = 0.0 
  if(neutral_horz_mix_bdy) then
      write (stdoutunit,'(a)') &
      ' Adding horizontal diffusivity in neutral boundary.'
      write (stdoutunit,'(a)') &
      ' This method is implemented only for the case with neutral_physics_simple=.false.'
      if(vel_micom_bdy > 0.0) then
         ah_array(:,:) = vel_micom_bdy*grid_length(:,:)
      else 
         ah_array(:,:) = ah_bdy
      endif 
      do j=jsc,jec
         write (stdoutunit,'(a,i4,a,e14.7,a)') &
         ' ah_array in neutral bdy layer at (isc,',j,',1)= ',ah_array(isc,j),' m^2/s'
      enddo
  endif


  ! Bryan-Lewis profile for aredi_array
  ! this profile is implemented for legacy purposes  
  if(bryan_lewis_aredi) then
    write(stdoutunit,'(/a)')' Using Bryan-Lewis depth profile for the Redi diffusivity.'
    write(stdoutunit,'(a)') ' The Bryan-Lewis profile has been tested ONLY with agm==0.'
    write(stdoutunit,'(a)') ' Vertical dependence to the neutral diffusivity has not been implemented'
    write(stdoutunit,'(a)') ' according to a physical theory giving a vertical structure depending'
    write(stdoutunit,'(a)') ' on the flow field. More theoretical work is required.'
    do k=1,nk
      aredi_array(:,:,k) = (ahb + (ahs - ahb)*exp(-Grd%zt(k)/500.0))
    enddo
    do k=1,nk
      write (stdoutunit,'(a,i3,e16.8)') '  k, diffusivity = ', k, aredi_array(isc,jsc,k)
    enddo
  endif

  ! Set diffusivity according to agm_lat_bands
  if(agm_lat_bands) then
      write(stdoutunit,'(1x,a)')    &
      '==>Note: Setting agm_array according to agm_lat_bands.'
      write(stdoutunit,'(a,e10.5)') &
      '      The ratio agm(south)/agm(north) is given by ',agm_lat_bands_boundary
      write(stdoutunit,'(a,e10.5)') &
      '      with the latitude separating the bands given by ',agm_lat_bands_ratio
      if(agm_lat_bands_boundary <= -90.0) then 
          write(stdoutunit,'(1x,a)') &
           '      Since agm_lat_bands_boundary <= -90, will default to globally constant agm value' 
      endif
      agm_array_local(:,:,:) = agm*Grd%tmask(:,:,:)
      if(agm_lat_bands_boundary > -90.) then 
          do j=jsc-Dom%yhalo,jec+Dom%yhalo
             do i=isc-Dom%xhalo,iec+Dom%xhalo
                if(Grd%yt(i,j) <= agm_lat_bands_boundary) then
                  agm_array_local(i,j,:) = agm*agm_lat_bands_ratio*Grd%tmask(i,j,:)
                endif 
             enddo
          enddo
      endif
  endif

  ! grid-scale dependent diffusivity suggested by that commonly used in MICOM
  ! vel_micom (m/s) sets velocity scale.
  ! space scale is set by grid size. 
  if(tracer_mix_micom) then
      do k=1,nk
        agm_array_local(:,:,k) = agm_micom(:,:)
      enddo
      do j=jsc,jec
         write (stdoutunit,'(a,i4,a,e14.7,a)') &
         ' Micom agm_array at (isc,',j,',1) = ',agm_array_local(isc,j,1),' m^2/s'
      enddo
      if(aredi_equal_agm) then
          aredi_array(:,:,:) = agm_array_local(:,:,:) 
          write (stdoutunit,'(a)') &
          ' aredi_array = agm_array as given by Micom grid scale dependent diffusivity'
      else
          write (stdoutunit,'(a)') ' Since (agm/=aredi), aredi_array=aredi'
      endif
  endif

  ! to have the halos filled 
  agm_array(:,:,:) = agm_array_local(:,:,:)

  ! make mandatory=.false. to facilitate backward compatibility with older restart files. 
  file_name = 'ocean_neutral.res.nc'
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'agm_array', agm_array, &
       domain=Dom%domain2d, mandatory=.false.) 
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'aredi_array', aredi_array, &
       domain=Dom%domain2d, mandatory=.false.)
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'rossby_radius', rossby_radius, &
       domain=Dom%domain2d, mandatory=.false.)
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'rossby_radius_raw', rossby_radius_raw, &
       domain=Dom%domain2d, mandatory=.false.)
  id_restart = register_restart_field(Nphysics_util_restart, file_name, 'bczone_radius', bczone_radius, &
       domain=Dom%domain2d, mandatory=.false.)

  if(Time%init .and. nphysics_util_zero_init) then 

     write(stdoutunit,'(1x,a)') ' Starting ocean_nphysics_util fields from raw initialization.'

  elseif(.NOT. file_exist('INPUT/ocean_neutral.res.nc')) then

     if (.NOT. Time%init) call mpp_error(FATAL,'Expecting file INPUT/ocean_neutral.res.nc to exist.&
             &This file was not found and Time%init=.false.')

  else

      call restore_state(Nphysics_util_restart)

      call mpp_update_domains(agm_array,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'agm_array being read from restart.'
      write (stdoutunit, *) 'checksum start agm_array', &
         mpp_chksum(agm_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))

      call mpp_update_domains(aredi_array,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'aredi_array being read from restart.'
      write (stdoutunit, *) 'checksum start aredi_array', &
         mpp_chksum(aredi_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))

      call mpp_update_domains(rossby_radius,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'rossby_radius being read from restart.'
      write (stdoutunit, *) 'checksum start rossby_radius', &
      mpp_chksum(rossby_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))

      call mpp_update_domains(rossby_radius_raw,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'rossby_radius_raw being read from restart.'
      write (stdoutunit, *) 'checksum start rossby_radius_raw', &
       mpp_chksum(rossby_radius_raw(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))

      call mpp_update_domains(bczone_radius,Dom%domain2d)
      write (stdoutunit,'(1x,a)') 'bczone_radius being read from restart.'
      write (stdoutunit, *) 'checksum start bczone_radius', &
        mpp_chksum(bczone_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))

  endif 

  ! for those cases with agm_array determined by 
  ! time invariant initialization methods.  
  if(.not. agm_closure) then
      if(.not. agm_read_restart) then
          write(stdoutunit,'(1x,a)') &
         'agm_closure=.false. and agm_read_restart=.false. => agm_array set to static profiles.'
          agm_array(:,:,:) = agm_array_local(:,:,:)
      else
          write(stdoutunit,'(1x,a)') &
          'agm_closure=.false. and agm_read_restart=.true. => agm_array set to restart values. '
      endif
  endif

  if(aredi_equal_agm) then 
      write(stdoutunit,'(1x,a)') &
      'aredi_equal_agm=.true. will force aredi_array to equal agm_array'
      aredi_array(:,:,:) = agm_array(:,:,:)
  else 
      write(stdoutunit,'(1x,a)') &
      'aredi_equal_agm=.false. allows aredi_array to differ from agm_array'
  endif 

  ! compute maximum stable neutral slope available for neutral diffusion
  i_delta = isc; j_delta = jsc; k_delta = 1
  delta_iso = 1e30
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if(Grd%tmask(i,j,k) /= 0.0) then 
          ft     = 0.5/(aredi_array(i,j,k)*dtime + epsln)
          deltaX = Grd%dxt(i,j)*Thickness%dzt(i,j,k)*ft
          deltaY = Grd%dyt(i,j)*Thickness%dzt(i,j,k)*ft
          if (delta_iso >= deltaX .or. delta_iso >= deltaY) then
            i_delta = i; j_delta = j; k_delta = k
            delta_iso = min(deltaX,deltaY)
          endif
        endif  
      enddo
    enddo
  enddo
  delta_iso  = delta_iso + 1.e-6*mpp_pe() ! to separate redundancies
  delta_iso0 = delta_iso
  call mpp_min (delta_iso)

  ! show most unstable location
  if (delta_iso == delta_iso0) then
     write(unit,'(a)')
     write(unit,'(a)')'---Neutral direction slope check I for linear stability of neutral diffusion---'
     write(unit,'(a,e14.7)')'With a neutral physics time step (secs) of ', dtime
     write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
     write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                   = ',Grd%xt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                   = ',Grd%yt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,') = ', &
                                           Thickness%dzt(i_delta,j_delta,k_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'aredi(',i_delta,',',j_delta,',',k_delta,') = ', &
                                           aredi_array(i_delta,j_delta,k_delta)
     write(unit,'(a,e14.7,a)')'delta_iso           = ',delta_iso,' is the maximum neutral direction slope' 
     write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
     write(unit,'(a)')'The namelist parameter smax should conservatively be <= delta_iso.'
     if(smax >= delta_iso) then 
        write(unit,'(a,f10.5,a)')'==> Warning: The namelist parameter smax= ',smax, ' is >= to delta_iso.'
        write(unit,'(a)')'Linear stability of the neutral diffusion scheme may be compromised.'
     endif
     write(unit,'(a)')
  endif


  ! Compute maximum diffusivity available given a maximum slope of smax
  i_delta = isc; j_delta = jsc; k_delta = 1
  ft = 0.5/(smax*dtime + epsln)
  delta_min = Thickness%dzt(i_delta,j_delta,k_delta)*Grd%dxt(i_delta,j_delta)
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if(Grd%tmask(i,j,k) /= 0.0) then        
          delta = min(Grd%dxt(i,j),Grd%dyt(i,j))*Thickness%dzt(i,j,k)
          if (delta_min > delta) then
            i_delta = i; j_delta = j; k_delta = k
            delta_min = delta
          endif
        endif   
      enddo
    enddo
  enddo
  A_max  = ft*delta_min + 1.e-6*mpp_pe() ! to separate redundancies
  A_max0 = A_max
  call mpp_min (A_max)

  ! show most unstable location
  if (A_max == A_max0) then
     write(unit,'(a)')
     write(unit,'(a)')'---Neutral direction slope check II for linear stability of neutral diffusion---'
     write(unit,'(a,e14.7)')'Assuming maximum Redi neutral diffusion slope of ', smax
     write(unit,'(a,e14.7)')'and neutral physics time step (secs) of ', dtime
     write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
     write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                  = ',&
     Grd%xt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                  = ',&
     Grd%yt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,')= ',&
     Thickness%dzt(i_delta,j_delta,k_delta)
     write(unit,'(a,e14.7,a)')'A_max      = ',A_max,' (m^2/sec) is the maximum neutral diffusivity'
     write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
     write(unit,'(a)')'Conservatively, neutral diffusivities used in the model should be less than A_max.'
     write(unit,'(a)')'--------------------------------------------------------------------------------'
     write(unit,'(a)')
  endif


end subroutine ocean_nphysics_coeff_init
! </SUBROUTINE>  NAME="ocean_nphysics_coeff_init"



!#######################################################################
! <SUBROUTINE NAME="tracer_derivs">
!
! <DESCRIPTION>
! Compute the tracer derivatives.
!
! Comments about smoothing drhodz:
!
! 1/ Tests in coupled 1-degree model showed extreme sensitivity 
! of MOC to smoothing.  GFDL users generally do NOT smooth, hence
! the default drhodz_smooth_vert=drhodz_smooth_horz=.false. 
!
! 2/ Smoothing the vertical derivative of drhodzb and drhodzh helps  
! is greatly needed for producing a regularized (i.e., well behaved)
! neutral slope vector.  
!
! 3/ An attempt was made to smooth dTdz and dSdz rather 
! than drhodz.  The resulting slope was smooth, but not as 
! smooth as when acting on drhodz itself.
!
! </DESCRIPTION>
!
subroutine tracer_derivs(Time, taum1, dtime, drhodT, drhodS, T_prog, dzwtr, &
                         dTdx, dTdy, dTdz, drhodzb, drhodzh)

  type(ocean_time_type),           intent(in)    :: Time
  integer,                         intent(in)    :: taum1
  real,                            intent(in)    :: dtime
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),    intent(in)    :: drhodS
  real, dimension(isd:,jsd:,0:),   intent(in)    :: dzwtr
  type(ocean_prog_tracer_type),    intent(in)    :: T_prog(:)

  type(tracer_3d_1_nk_type),       intent(inout) :: dTdx(:)
  type(tracer_3d_1_nk_type),       intent(inout) :: dTdy(:)
  type(tracer_3d_0_nk_type),       intent(inout) :: dTdz(:)
  real, dimension(isd:,jsd:,:,0:), intent(inout) :: drhodzb
  real, dimension(isd:,jsd:,:,0:), intent(inout) :: drhodzh


  integer :: i, j, k, m, n
  integer :: kp1, kbot
  integer :: kr, kpkr, km1pkr
  real    :: dTdz_ijk1, dTdz_ijk2
  real    :: drhodzb0_prev, drhodzb1_prev
  real    :: drhodzh0_prev, drhodzh1_prev
  real    :: tmpb0, tmpb1, tmph0, tmph1

  call mpp_clock_begin(id_clock_tracer_derivs)

  ! linear map to constant depth level horizontal derivatives (constant z-level);  
  ! this approach can introduce spurious extrema and so is not recommended, except
  ! perhaps for sigma-models (WARNING: neutral physics has not been tested 
  ! w/ sigma vertical coordinates).  
  if(horz_z_derivative) then 
      do n=1,num_prog_tracers
         wrk1(:,:,1) = T_prog(n)%field(:,:,1,taum1)
         do k=1,nk
            kp1 = min(k+1,nk)
            wrk1(:,:,2)          = T_prog(n)%field(:,:,k,taum1)
            dTdx(n)%field(:,:,k) = FDX_ZT(wrk1(:,:,1:2),k)*FMX(Grd%tmask(:,:,k))
            dTdy(n)%field(:,:,k) = FDY_ZT(wrk1(:,:,1:2),k)*FMY(Grd%tmask(:,:,k))
            wrk1(:,:,1)          = wrk1(:,:,2)
         enddo
      enddo

  ! horizontal derivatives taken along surfaces of 
  ! constant vertical coordinate (constant k-level)
  elseif(horz_s_derivative) then 
      do n=1,num_prog_tracers
         do k=1,nk
            kp1 = min(k+1,nk)
            wrk1(:,:,1)          = T_prog(n)%field(:,:,k,taum1)
            dTdx(n)%field(:,:,k) = FDX_T(wrk1(:,:,1))*FMX(Grd%tmask(:,:,k))
            dTdy(n)%field(:,:,k) = FDY_T(wrk1(:,:,1))*FMY(Grd%tmask(:,:,k))
         enddo
      enddo
  endif

  ! vertical derivative 
  do n=1,num_prog_tracers
     do k=1,nk
        kp1 = min(k+1,nk)
        do j=jsd,jed
           do i=isd,ied
              dTdz(n)%field(i,j,k) = Grd%tmask(i,j,kp1)*dzwtr(i,j,k) &
                *(T_prog(n)%field(i,j,k,taum1)-T_prog(n)%field(i,j,kp1,taum1))
           enddo
        enddo
     enddo
  enddo

  ! vertical neutral density derivative for use in fz_terms
  ! and fz_flux, and for use in fx_flux and fy_flux. 
  ! note that the derivative at k=nk vanishes by definition
  ! since these derivatives are at the bottom of tracer cell. 
  ! also note the use of -epsln_drhodz ensures the vertical 
  ! derivative is always < 0.  We also support the same 
  ! approach used in the mom4p0d code for legacy purposes. 
  if(drhodz_mom4p1) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied

               dTdz_ijk1 = dTdz(index_temp)%field(i,j,k) ; dTdz_ijk2 = dTdz(index_salt)%field(i,j,k)

               kr=0
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 + drhodS(i,j,kpkr)*dTdz_ijk2 
               drhodzb(i,j,k,kr) = min(drhodzb(i,j,k,kr), -epsln_drhodz)

               kr=1
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 + drhodS(i,j,kpkr)*dTdz_ijk2
               drhodzb(i,j,k,kr) = min(drhodzb(i,j,k,kr), -epsln_drhodz)


               kr=0
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dTdz(index_salt)%field(i,j,km1pkr)
               drhodzh(i,j,k,kr) = min(drhodzh(i,j,k,kr), -epsln_drhodz)

               kr=1
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dTdz(index_salt)%field(i,j,km1pkr)
               drhodzh(i,j,k,kr) = min(drhodzh(i,j,k,kr), -epsln_drhodz)

            enddo
         enddo
      enddo

  else 

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied

               dTdz_ijk1 = dTdz(index_temp)%field(i,j,k) ; dTdz_ijk2 = dTdz(index_salt)%field(i,j,k)

               kr=0
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 &
                                 + drhodS(i,j,kpkr)*dTdz_ijk2 -epsln

               kr=1
               kpkr = min(k+kr,nk)
               drhodzb(i,j,k,kr) = drhodT(i,j,kpkr)*dTdz_ijk1 &
                                 + drhodS(i,j,kpkr)*dTdz_ijk2 -epsln

               kr=0
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dTdz(index_salt)%field(i,j,km1pkr) -epsln

               kr=1
               km1pkr = max(k-1+kr,1)
               drhodzh(i,j,k,kr) = drhodT(i,j,k)*dTdz(index_temp)%field(i,j,km1pkr) &
                                 + drhodS(i,j,k)*dTdz(index_salt)%field(i,j,km1pkr) -epsln

            enddo
         enddo
      enddo
  endif

  ! vertically smooth the vertical derivative of density to 
  ! produce a smooth neutral slope vector for all flux components.
  if(drhodz_smooth_vert) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied

               drhodzb0_prev = onefourth*drhodzb(i,j,1,0)
               drhodzb1_prev = onefourth*drhodzb(i,j,1,1)
               drhodzh0_prev = onefourth*drhodzh(i,j,1,0)
               drhodzh1_prev = onefourth*drhodzh(i,j,1,1)

               kbot=Grd%kmt(i,j)
               if(kbot>1) then
                   do k=2,kbot-2

                      tmpb0            = drhodzb(i,j,k,0)
                      drhodzb(i,j,k,0) = drhodzb0_prev + onehalf*drhodzb(i,j,k,0) &
                                                       + onefourth*drhodzb(i,j,k+1,0)
                      drhodzb0_prev    = onefourth*tmpb0

                      tmpb1            = drhodzb(i,j,k,1)
                      drhodzb(i,j,k,1) = drhodzb1_prev + onehalf*drhodzb(i,j,k,1) &
                                                       + onefourth*drhodzb(i,j,k+1,1)
                      drhodzb1_prev    = onefourth*tmpb1

                      tmph0            = drhodzh(i,j,k,0)
                      drhodzh(i,j,k,0) = drhodzh0_prev + onehalf*drhodzh(i,j,k,0) &
                                                       + onefourth*drhodzh(i,j,k+1,0)
                      drhodzh0_prev    = onefourth*tmph0

                      tmph1            = drhodzh(i,j,k,1)
                      drhodzh(i,j,k,1) = drhodzh1_prev + onehalf*drhodzh(i,j,k,1) &
                                                       + onefourth*drhodzh(i,j,k+1,1)
                      drhodzh1_prev     = onefourth*tmph1

                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! horizontally smooth the vertical derivative of density to 
  ! produce a smooth neutral slope vector for all flux components.
  ! since drhodzb is needed only on the computational domain, there 
  ! is no need to update its values in the halos.  
  if(drhodz_smooth_horz) then 
      do k=1,nk-1
         drhodzb(:,:,k,0) = drhodzb(:,:,k,0) + dtime*LAP_T(drhodzb(:,:,k,0),smooth_lap(:,:))
         drhodzb(:,:,k,1) = drhodzb(:,:,k,1) + dtime*LAP_T(drhodzb(:,:,k,1),smooth_lap(:,:))
         drhodzh(:,:,k,0) = drhodzh(:,:,k,0) + dtime*LAP_T(drhodzh(:,:,k,0),smooth_lap(:,:))
         drhodzh(:,:,k,1) = drhodzh(:,:,k,1) + dtime*LAP_T(drhodzh(:,:,k,1),smooth_lap(:,:))
      enddo
      call mpp_update_domains(drhodzh(:,:,:,0), Dom%domain2d, complete=.false.) 
      call mpp_update_domains(drhodzh(:,:,:,1), Dom%domain2d, complete=.true.) 
  endif

  ! compute squared buoyancy frequency 
  wrk2(:,:,:)=0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk2(i,j,k) = -onehalf*grav*rho0r*(drhodzb(i,j,k,0)+drhodzb(i,j,k,1))
        enddo
     enddo
  enddo
  if(id_N2slope > 0) then 
      used = send_data (id_N2slope, wrk2(:,:,:),    &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk) 
  endif



  call mpp_clock_end(id_clock_tracer_derivs)

end subroutine tracer_derivs
! </SUBROUTINE> NAME="tracer_derivs"



!#######################################################################
! <SUBROUTINE NAME="neutral_slopes">
!
! <DESCRIPTION>
! Subroutine computes the neutral slopes for the triads associated 
! with the vertical flux component.  
!
! Array tensor_31 initially holds the x-slope used for flux component fz.
! Array tensor_32 initially holds the y-slope used for flux component fz.
!
! In subsequent calculations, these arrays will be multipied by the
! diffusivities.  
!
! No slope tapering is applied in this routine. 
!
! slopes are computed over k=1,nk-1, since the slope at k=nk 
! should be zero. 
!
! </DESCRIPTION>
!
subroutine neutral_slopes(Time, dTdx, dTdy, drhodT, drhodS, drhodzb, tensor_31, tensor_32)

  type(ocean_time_type),              intent(in)    :: Time
  type(tracer_3d_1_nk_type),          intent(in)    :: dTdx(:)
  type(tracer_3d_1_nk_type),          intent(in)    :: dTdy(:)
  real, dimension(isd:,jsd:,:),       intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),       intent(in)    :: drhodS
  real, dimension(isd:,jsd:,:,0:),    intent(in)    :: drhodzb
  real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: tensor_31
  real, dimension(isd:,jsd:,:,0:,0:), intent(inout) :: tensor_32

  integer :: i, j, k, kp1
  real :: drhodT_ijk, drhodS_ijk, drhodT_ijkp1, drhodS_ijkp1
  real :: tmask_ijkp1

  real :: dTdx_ijk1, dTdx_ijk2, dTdx_im1jk1, dTdx_im1jk2
  real :: dTdx_ijkp11, dTdx_ijkp12, dTdx_im1jkp11, dTdx_im1jkp12
  real :: dTdy_ijk1, dTdy_ijk2, dTdy_ijm1k1, dTdy_ijm1k2
  real :: dTdy_ijkp11, dTdy_ijkp12, dTdy_ijm1kp11, dTdy_ijm1kp12
  real :: drhodzbr_ijk0, drhodzbr_ijk1

  call mpp_clock_begin(id_clock_neutral_slopes)

  tensor_31(:,:,:,:,:) = 0.0
  tensor_32(:,:,:,:,:) = 0.0

  do k=1,nk-1
     kp1 = k+1
     do j=jsc,jec
        do i=isc,iec

           tmask_ijkp1    = Grd%tmask(i,j,kp1)

           drhodT_ijk     = drhodT(i,j,k)
           drhodS_ijk     = drhodS(i,j,k)
           drhodT_ijkp1   = drhodT(i,j,kp1)
           drhodS_ijkp1   = drhodS(i,j,kp1)

           dTdx_ijk1     = dTdx(index_temp)%field(i,j,k)
           dTdx_ijk2     = dTdx(index_salt)%field(i,j,k)
           dTdx_im1jk1   = dTdx(index_temp)%field(i-1,j,k)
           dTdx_im1jk2   = dTdx(index_salt)%field(i-1,j,k)
           dTdx_ijkp11   = dTdx(index_temp)%field(i,j,kp1)
           dTdx_ijkp12   = dTdx(index_salt)%field(i,j,kp1)
           dTdx_im1jkp11 = dTdx(index_temp)%field(i-1,j,kp1)
           dTdx_im1jkp12 = dTdx(index_salt)%field(i-1,j,kp1)

           dTdy_ijk1     = dTdy(index_temp)%field(i,j,k)
           dTdy_ijk2     = dTdy(index_salt)%field(i,j,k)
           dTdy_ijm1k1   = dTdy(index_temp)%field(i,j-1,k)
           dTdy_ijm1k2   = dTdy(index_salt)%field(i,j-1,k)
           dTdy_ijkp11   = dTdy(index_temp)%field(i,j,kp1)
           dTdy_ijkp12   = dTdy(index_salt)%field(i,j,kp1)
           dTdy_ijm1kp11 = dTdy(index_temp)%field(i,j-1,kp1)
           dTdy_ijm1kp12 = dTdy(index_salt)%field(i,j-1,kp1)

           drhodzbr_ijk0 = 1.0/drhodzb(i,j,k,0)
           drhodzbr_ijk1 = 1.0/drhodzb(i,j,k,1)


           ! ip=jq=0

           !   kr=0
           tensor_31(i,j,k,0,0) = -tmask_ijkp1                     &
                *(drhodT_ijk*dTdx_im1jk1 + drhodS_ijk*dTdx_im1jk2) &
                *drhodzbr_ijk0
           tensor_32(i,j,k,0,0) = -tmask_ijkp1                     &
                *(drhodT_ijk*dTdy_ijm1k1 + drhodS_ijk*dTdy_ijm1k2) &
                *drhodzbr_ijk0

           !   kr=1 
           tensor_31(i,j,k,0,1) = -tmask_ijkp1                             &
                *(drhodT_ijkp1*dTdx_im1jkp11 + drhodS_ijkp1*dTdx_im1jkp12) &
                *drhodzbr_ijk1

           tensor_32(i,j,k,0,1) = -tmask_ijkp1                             &
                *(drhodT_ijkp1*dTdy_ijm1kp11 + drhodS_ijkp1*dTdy_ijm1kp12) &
                *drhodzbr_ijk1


           ! ip=jq=1

           !   kr=0 
           tensor_31(i,j,k,1,0) = -tmask_ijkp1                 &
                *(drhodT_ijk*dTdx_ijk1 + drhodS_ijk*dTdx_ijk2) &
                *drhodzbr_ijk0

           tensor_32(i,j,k,1,0) = -tmask_ijkp1                 &
                *(drhodT_ijk*dTdy_ijk1 + drhodS_ijk*dTdy_ijk2) &
                *drhodzbr_ijk0
           !   kr=1 
           tensor_31(i,j,k,1,1) = -tmask_ijkp1                         &
                *(drhodT_ijkp1*dTdx_ijkp11 + drhodS_ijkp1*dTdx_ijkp12) &
                *drhodzbr_ijk1

           tensor_32(i,j,k,1,1) = -tmask_ijkp1                         &
                *(drhodT_ijkp1*dTdy_ijkp11 + drhodS_ijkp1*dTdy_ijkp12) &
                *drhodzbr_ijk1

        enddo
     enddo
  enddo


  ! send to diagnostic manager 
  if (id_slope31 > 0) then 
       wrk1(isc:iec,jsc:jec,:) = onefourth*                                     &
         (tensor_31(isc:iec,jsc:jec,:,0,0) + tensor_31(isc:iec,jsc:jec,:,0,1) + &
          tensor_31(isc:iec,jsc:jec,:,1,0) + tensor_31(isc:iec,jsc:jec,:,1,1))
       used = send_data (id_slope31, wrk1(:,:,:),      &
              Time%model_time, rmask=Grd%tmask(:,:,:), &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_slope32 > 0) then 
       wrk1(isc:iec,jsc:jec,:) = onefourth*                                     &
         (tensor_32(isc:iec,jsc:jec,:,0,0) + tensor_32(isc:iec,jsc:jec,:,0,1) + &
          tensor_32(isc:iec,jsc:jec,:,1,0) + tensor_32(isc:iec,jsc:jec,:,1,1))
       used = send_data (id_slope32, wrk1(:,:,:),      &
              Time%model_time, rmask=Grd%tmask(:,:,:), &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif  

  call mpp_clock_end(id_clock_neutral_slopes)

end subroutine neutral_slopes
! </SUBROUTINE> NAME="neutral_slopes"



!#######################################################################
! <SUBROUTINE NAME="compute_eady_rate">
!
! <DESCRIPTION>
! Finish computing eady growth rate.
! </DESCRIPTION>
!
subroutine compute_eady_rate(Time, Thickness, T_prog, Dens, eady_termx, eady_termy)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  real,  dimension(isd:,jsd:,:),intent(in) :: eady_termx
  real,  dimension(isd:,jsd:,:),intent(in) :: eady_termy

  integer :: i, j, k, m, kbot, tau
  real    :: eady_rate_prev, tmp
  real    :: eady_mld(isd:ied,jsd:jed)

  if(.not. agm_closure) return 

  call mpp_clock_begin(id_clock_compute_eady_rate)

  tau           = Time%tau 
  wrk1(:,:,:)   = Grd%tmask(:,:,:)   
  eady_mld(:,:) = 0.0


  ! raw Eady growth rate, with regularization not yet applied 
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           eady_rate(i,j,k) = Grd%tmask(i,j,k) &
            *sqrt((eady_termx(i,j,k)*count_x(i,j)+eady_termy(i,j,k)*count_y(i,j))*gravrho0r)
        enddo
     enddo
  enddo

  ! apply cap to the Eady growth rate
  if(agm_closure_eady_cap) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               eady_rate(i,j,k) = 2.0*eady_rate(i,j,k)*agm_growth_rate_max(i,j) &
                                  /(epsln+eady_rate(i,j,k)+agm_growth_rate_max(i,j))
            enddo
         enddo
      enddo
  endif


  ! vertically average Eady within surface mixed layer 
  if(agm_closure_eady_ave_mixed) then

      call calc_mixed_layer_depth(Thickness,                &
           T_prog(index_salt)%field(isd:ied,jsd:jed,:,tau), &
           T_prog(index_temp)%field(isd:ied,jsd:jed,:,tau), &
           Dens%rho(isd:ied,jsd:jed,:,tau),                 &
           Dens%pressure_at_depth(isd:ied,jsd:jed,:),       &
           Time%model_time, eady_mld)

      do j=jsc,jec
         do i=isc,iec
            eady_mld(i,j) = Grd%tmask(i,j,1)*min(eady_mld(i,j), Grd%ht(i,j))
         enddo
      enddo

      ! do not believe eady_rate at k=1, so always skip its contribution 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      do k=2,nk
         do j=jsc,jec
            do i=isc,iec
               if(Thickness%depth_zt(i,j,k) < eady_mld(i,j)) then 
                   wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%dzt(i,j,k)
                   wrk2_2d(i,j) = wrk2_2d(i,j) + eady_rate(i,j,k)*Thickness%dzt(i,j,k)
               endif
            enddo
         enddo
      enddo
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               if(Thickness%depth_zt(i,j,k) <= eady_mld(i,j)) then 
                   eady_rate(i,j,k) = Grd%tmask(i,j,k)*wrk2_2d(i,j)/(wrk1_2d(i,j)+epsln)
               endif
            enddo
         enddo
      enddo

  endif ! endif for agm_closure_eady_ave_mixed


  ! apply vertical 1-2-1 smoothing 
  if(agm_closure_eady_smooth_vert) then 
      do m=1,num_121_passes
         do j=jsc,jec
            do i=isc,iec
               eady_rate_prev = onefourth*eady_rate(i,j,1)
               kbot=Grd%kmt(i,j)
               if(kbot>1) then
                   do k=2,kbot-2
                      tmp              = eady_rate(i,j,k)
                      eady_rate(i,j,k) = eady_rate_prev + onehalf*eady_rate(i,j,k) &
                                                        + onefourth*eady_rate(i,j,k+1)
                      eady_rate_prev   = onefourth*tmp
                   enddo
               endif
            enddo
         enddo
      enddo
  endif 

  ! apply horizontal 1-2-1 smoothing 
  if(agm_closure_eady_smooth_horz) then 
      call mpp_update_domains (eady_rate(:,:,:), Dom%domain2d) 
      do k=1,nk
         eady_rate(:,:,k) = S2D(eady_rate(:,:,k))
      enddo
  endif 

  ! compute vertical average over specified depth 
  eady_rate_zave(:,:) = 0.0
  wrk1_2d(:,:)        = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           if(Thickness%depth_zwt(i,j,k) >= agm_closure_upper_depth .and. &
              Thickness%depth_zwt(i,j,k) <= agm_closure_lower_depth) then 
              eady_rate_zave(i,j) = eady_rate_zave(i,j) + eady_rate(i,j,k)*Thickness%dzwt(i,j,k)
              wrk1_2d(i,j)        = wrk1_2d(i,j) + Thickness%dzwt(i,j,k)
           endif
        enddo
     enddo
  enddo
  do j=jsc,jec
     do i=isc,iec
        eady_rate_zave(i,j) = Grd%tmask(i,j,1)*eady_rate_zave(i,j)/(wrk1_2d(i,j)+epsln)
     enddo
  enddo

  if (id_eady_rate > 0) used = send_data (id_eady_rate, eady_rate(:,:,:), &
                               Time%model_time, rmask=Grd%tmask(:,:,:),   &
                               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_eady_rate_zave > 0) used = send_data (id_eady_rate_zave, eady_rate_zave(:,:), &
                             Time%model_time, rmask=Grd%tmask(:,:,1),                  &
                             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_eady_mld > 0) used = send_data (id_eady_mld, eady_mld(:,:),   &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


  call mpp_clock_end(id_clock_compute_eady_rate)

end subroutine compute_eady_rate
! </SUBROUTINE> NAME="compute_eady_rate"



!#######################################################################
! <SUBROUTINE NAME="compute_baroclinicity">
!
! <DESCRIPTION>
! Finish computing baroclinicity, which is defined to be the vertically
! averaged magnitude of the horizontal density gradient.
! </DESCRIPTION>
!
subroutine compute_baroclinicity(model_time, baroclinic_termx, baroclinic_termy)

  type(time_type),              intent(in)    :: model_time
  real,  dimension(isd:,jsd:),  intent(in)    :: baroclinic_termx
  real,  dimension(isd:,jsd:),  intent(in)    :: baroclinic_termy

  integer :: i, j
  real    :: vertical_range

  if(.not. agm_closure) return 

  call mpp_clock_begin(id_clock_compute_baroclinicity)
  
  vertical_range = epsln + agm_closure_lower_depth - agm_closure_upper_depth

  do j=jsc,jec
     do i=isc,iec
        baroclinicity(i,j) = Grd%tmask(i,j,1)*                                       &
             (baroclinic_termx(i,j)*count_x(i,j)+baroclinic_termy(i,j)*count_y(i,j)) &
             /vertical_range
     enddo
  enddo

  if (id_baroclinicity > 0) used = send_data (id_baroclinicity, baroclinicity(:,:), &
                                   model_time, rmask=Grd%tmask(:,:,1),              &
                                   is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  call mpp_clock_end(id_clock_compute_baroclinicity)

end subroutine compute_baroclinicity
! </SUBROUTINE> NAME="compute_baroclinicity"


!#######################################################################
! <SUBROUTINE NAME="compute_rossby_radius">
!
! <DESCRIPTION>
! Subroutine computes the first baroclinic Rossby radius of deformation. 
! Employ WKB approach described by Chelton et al.  In particular, 
! use formulae (2.2), (2.3a) and (2.3b) from their paper. 
!
! Place a max and min value on the Rossby radius.
!
! Compute buoyancy frequency in terms of vertical gradient of 
! locally referenced potential density.  Place the reference point
! at the interface between the tracer cells, which is also where 
! the vertical derivative of neutral density is located.  This amounts 
! to a centered difference computation similar to that used by 
! Chelton et al. equation (B.4). 
! </DESCRIPTION>
!
subroutine compute_rossby_radius(Thickness, dTdz, model_time, drhodT, drhodS, &
                                 rossby_radius, rossby_radius_raw)

  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(tracer_3d_0_nk_type),     intent(in)    :: dTdz(:)
  type(time_type),               intent(in)    :: model_time
  real, dimension(isd:,jsd:,:),  intent(in)    :: drhodT
  real, dimension(isd:,jsd:,:),  intent(in)    :: drhodS
  real, dimension(isd:,jsd:),    intent(inout) :: rossby_radius
  real, dimension(isd:,jsd:),    intent(inout) :: rossby_radius_raw

  real     :: drhodzb_speed
  integer  :: i,j,k

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_clock_compute_rossby_radius)
 
  gravity_wave_speed(:,:) = 0.0
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           drhodzb_speed  = 0.5*( (drhodT(i,j,k) + drhodT(i,j,k+1))*dTdz(index_temp)%field(i,j,k) &
                                 +(drhodS(i,j,k) + drhodS(i,j,k+1))*dTdz(index_salt)%field(i,j,k) )
           gravity_wave_speed(i,j) = gravity_wave_speed(i,j) &
                                   + Thickness%dzwt(i,j,k)*sqrt(abs(gravrho0r*drhodzb_speed))
        enddo
     enddo
  enddo
  gravity_wave_speed(:,:) = gravity_wave_speed(:,:)/pi

  ! for Eden near the equator 
  do j=jsd,jed
     do i=isd,ied
        sqrt2betaCr(i,j) = sqrt(2.0*gravity_wave_speed(i,j)*beta_param(i,j))
     enddo
  enddo

  wrk1_2d(:,:)           = 0.0
  wrk2_2d(:,:)           = 0.0
  rossby_radius_raw(:,:) = 0.0
  do j=jsc,jec
     do i=isc,iec
        wrk1_2d(i,j)           = gravity_wave_speed(i,j)/(coriolis_param(i,j) + epsln)
        wrk2_2d(i,j)           = sqrt(gravity_wave_speed(i,j)/(2.0*beta_param(i,j)+epsln))
        rossby_radius_raw(i,j) = min(wrk1_2d(i,j), wrk2_2d(i,j))
     enddo
  enddo

  agm_grid_scaling(:,:) = 1.0
  if(agm_closure_grid_scaling) then 

      ! for backward bitwise compatibility 
      if(agm_closure_grid_scaling_power==2.0) then 
          do j=jsc,jec
             do i=isc,iec
                agm_grid_scaling(i,j) =  Grd%tmask(i,j,1) &
                *grid_length(i,j)**2/(grid_length(i,j)**2 + rossby_radius_raw(i,j)**2)
             enddo
          enddo
      else  
          do j=jsc,jec
             do i=isc,iec
                agm_grid_scaling(i,j) =  Grd%tmask(i,j,1)                 &
                     *grid_length(i,j)**agm_closure_grid_scaling_power    &
                     /( grid_length(i,j)**agm_closure_grid_scaling_power  &
                     + rossby_radius_raw(i,j)**agm_closure_grid_scaling_power)
             enddo
          enddo
      endif

  endif

  aredi_grid_scaling(:,:) = 1.0
  if(aredi_diffusivity_grid_scaling) then 
      do j=jsc,jec
         do i=isc,iec
            aredi_grid_scaling(i,j) = Grd%tmask(i,j,1)              &
                  *grid_length(i,j)**agm_closure_grid_scaling_power &
                /( grid_length(i,j)**agm_closure_grid_scaling_power &
                 + rossby_radius_raw(i,j)**agm_closure_grid_scaling_power)
         enddo
      enddo
      call mpp_update_domains (aredi_grid_scaling(:,:), Dom%domain2d)
  endif

  do j=jsc,jec
     do i=isc,iec
        rossby_radius(i,j) = min(rossby_radius_max,max(rossby_radius_min,rossby_radius_raw(i,j)))
     enddo
  enddo
   

  if (id_rossby_radius > 0) used = send_data (id_rossby_radius,rossby_radius_raw(:,:),&
                            model_time, rmask=Grd%tmask(:,:,1),                       &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_rossby > 0) used = send_data (id_rossby, rossby_radius(:,:), &
                            model_time, rmask=Grd%tmask(:,:,1),       &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_rossby_nonequator > 0) used = send_data (id_rossby_nonequator, wrk1_2d(:,:), &
                            model_time, rmask=Grd%tmask(:,:,1),                       &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_rossby_equator > 0) used = send_data (id_rossby_equator, wrk2_2d(:,:), &
                            model_time, rmask=Grd%tmask(:,:,1),                 &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_gw_speed > 0) used = send_data (id_gw_speed, gravity_wave_speed(:,:), &
                              model_time, rmask=Grd%tmask(:,:,1),              &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_sqrt2betaCr > 0) used = send_data (id_sqrt2betaCr, sqrt2betaCr(:,:), &
                              model_time, rmask=Grd%tmask(:,:,1),             &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_agm_grid_scaling > 0) used = send_data (id_agm_grid_scaling, agm_grid_scaling(:,:), &
                               model_time, rmask=Grd%tmask(:,:,1),                           &
                               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_aredi_grid_scaling > 0) used = send_data (id_aredi_grid_scaling, aredi_grid_scaling(:,:), &
                               model_time, rmask=Grd%tmask(:,:,1),                                 &
                               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


  if(debug_this_module) then 
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_nphysics_util'
      call write_timestamp(model_time)
      write (stdoutunit, *) 'checksum rossby_radius', &
      mpp_chksum(rossby_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
      write (stdoutunit, *) 'checksum rossby_radius_raw', &
      mpp_chksum(rossby_radius_raw(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
  endif

  call mpp_clock_end(id_clock_compute_rossby_radius)

end subroutine compute_rossby_radius
! </SUBROUTINE> NAME="compute_rossby_radius"



!#######################################################################
! <SUBROUTINE NAME="compute_bczone_radius">
!
! <DESCRIPTION>
! Subroutine computes the radius of the baroclinic zone in a manner 
! suggested by the Hadley Centre approach (Malcolm Roberts, personal 
! communication).  
!
! Algorithm is used in MOM3 and documented in the MOM3 Manual.
!
! </DESCRIPTION>
!
subroutine compute_bczone_radius(model_time, bczone_radius)

  type(time_type),            intent(in)    :: model_time
  real, dimension(isd:,jsd:), intent(inout) :: bczone_radius

  integer :: i,j
  real    :: n_zone, e_zone, s_zone, w_zone, fract, nstot, ewtot
  integer :: ip, jq

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if(.not. agm_closure) return
  if(.not. agm_closure_length_bczone) return 

  call mpp_clock_begin(id_clock_compute_bczone_radius)

  do j=jsc,jec
     do i=isc,iec
        bczone_rate(i,j) = eady_rate_zave(i,j)
     enddo
  enddo
  call mpp_update_domains (bczone_rate(:,:), BCzone_domain%domain2d)

  do j=jsc,jec
    do i=isc,iec
      if (bczone_rate(i,j) > agm_closure_bczone_crit_rate) then

        ! search northward 
        n_zone = bczone_dyt(i,j)
        do jq=j+1,j+bczone_max_pts
          if (bczone_rate(i,jq) > agm_closure_bczone_crit_rate) then
            n_zone = n_zone + bczone_dyt(i,jq)
          else
            exit
          endif
        enddo

        ! search southward 
        s_zone = bczone_dyt(i,j)
        do jq=j-1,j-bczone_max_pts,-1
          if (bczone_rate(i,jq) > agm_closure_bczone_crit_rate) then
            s_zone = s_zone + bczone_dyt(i,jq)
          else
            exit
          endif
        enddo

        ! search eastward 
        e_zone = bczone_dxt(i,j)
        do ip=i+1,i+bczone_max_pts
          if (bczone_rate(ip,j) > agm_closure_bczone_crit_rate) then
            e_zone = e_zone + bczone_dxt(ip,j)
          else
            exit
          endif
        enddo

        ! search westward 
        w_zone = bczone_dxt(i,j)
        do ip=i-1,i-bczone_max_pts,-1
          if (bczone_rate(ip,j) > agm_closure_bczone_crit_rate) then
            w_zone = w_zone + bczone_dxt(ip,j)
          else
            exit
          endif
        enddo

        ! total radius (subtraction accounts for double-counting central point)
        nstot=n_zone+s_zone-bczone_dyt(i,j)
        ewtot=e_zone+w_zone-bczone_dxt(i,j)

        if (nstot < ewtot) then 
          fract=min(n_zone,s_zone)/max(n_zone,s_zone) 
          bczone_radius(i,j)=fract*nstot 
        else   
          fract=min(e_zone,w_zone)/max(e_zone,w_zone)
          bczone_radius(i,j)=fract*ewtot
        endif
      endif
          
    enddo
  enddo

  if (id_bczone > 0) used = send_data (id_bczone, bczone_radius(:,:), &
                            model_time, rmask=Grd%tmask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


  if(debug_this_module) then 
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_nphysics_util'
      call write_timestamp(model_time)
      write (stdoutunit, *) &
           'checksum bczone_radius', mpp_chksum(bczone_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
  endif


  call mpp_clock_end(id_clock_compute_bczone_radius)

end subroutine compute_bczone_radius
! </SUBROUTINE> NAME="compute_bczone_radius"



!#######################################################################
! <SUBROUTINE NAME="compute_diffusivity">
!
! <DESCRIPTION>
! Subroutine computes flow dependent diffusivity.
! Allow for an added dimensionless tuning factor as well as a 
! minimum and maximum diffusivity. 
! </DESCRIPTION>
!
subroutine compute_diffusivity(model_time, ksurf_blayer, drhodz_zt, rossby_radius, &
                               bczone_radius, agm_array, aredi_array)

  type(time_type), intent(in)                 :: model_time
  integer, dimension(isd:,jsd:),intent(in)    :: ksurf_blayer
  real, dimension(isd:,jsd:,:), intent(in)    :: drhodz_zt
  real, dimension(isd:,jsd:),   intent(in)    :: rossby_radius
  real, dimension(isd:,jsd:),   intent(in)    :: bczone_radius
  real, dimension(isd:,jsd:,:), intent(inout) :: agm_array
  real, dimension(isd:,jsd:,:), intent(inout) :: aredi_array

  integer :: i, j, k, kp1
  real    :: denom, param, active_cells, tmp

  integer :: stdoutunit 
  stdoutunit=stdout() 

  call mpp_clock_begin(id_clock_compute_diffusivity)


  if(agm_closure) then 


     agm_growth_rate = 0.0
     agm_length      = 0.0
     wrk1            = 0.0 !growth_rate
     wrk2            = 0.0 !rhines length 
     wrk3            = 0.0 !agm_qg
     wrk4            = 0.0 !agm_fast
     wrk2_2d(:,:)    = 0.0 !drhodz_ref

     ! set the length scale 
     if(agm_closure_length_fixed .or. agm_closure_baroclinic) then 
         agm_length(:,:,:) = Grd%tmask(:,:,:)*agm_closure_length

     elseif(agm_closure_length_rossby) then 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  agm_length(i,j,k) = Grd%tmask(i,j,k)* &
                  2.0*grid_length(i,j)*rossby_radius(i,j)/(grid_length(i,j)+rossby_radius(i,j)+epsln)
               enddo
            enddo
         enddo

     elseif(agm_closure_length_bczone) then 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  agm_length(i,j,k) = Grd%tmask(i,j,k)* &
                  2.0*grid_length(i,j)*bczone_radius(i,j)/(grid_length(i,j)+bczone_radius(i,j)+epsln)
               enddo
            enddo
         enddo

     elseif(agm_closure_eden_greatbatch) then
         if(agm_closure_eden_length_const) then 
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      agm_length(i,j,k) = agm_closure_eden_length*Grd%tmask(i,j,k)
                   enddo
                enddo
             enddo
         else 
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      wrk2(i,j,k)       = eady_rate(i,j,k)/(beta_param(i,j)+epsln)
                      agm_length(i,j,k) = min(wrk2(i,j,k),rossby_radius(i,j))
                      agm_length(i,j,k) = Grd%tmask(i,j,k)* &
                      2.0*grid_length(i,j)*agm_length(i,j,k)/(grid_length(i,j)+agm_length(i,j,k)+epsln)
                   enddo
                enddo
             enddo
         endif
     endif

     if(agm_closure_length_cap) then 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  agm_length(i,j,k) = Grd%tmask(i,j,k)*min(agm_closure_length_max, agm_length(i,j,k))
               enddo
            enddo
         enddo
     endif


     ! set the growth rates 
     if(agm_closure_baroclinic) then 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp = gravrho0r_buoyr*baroclinicity(i,j)
                  wrk1(i,j,k) = 2.0*tmp*agm_growth_rate_max(i,j) &
                                /(tmp+agm_growth_rate_max(i,j)+epsln)  
               enddo
            enddo
         enddo

     elseif(agm_closure_eden_greatbatch) then 
         if(agm_closure_eden_gamma /= 0.0) then 
            do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      param = max(coriolis_param(i,j), sqrt2betaCr(i,j))
                      denom = sqrt(param**2 + agm_closure_eden_gamma*eady_rate(i,j,k)**2)+epsln
                      wrk1(i,j,k) = param*eady_rate(i,j,k)/denom
                   enddo
               enddo
            enddo
         else 
            do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      wrk1(i,j,k) = eady_rate(i,j,k)
                   enddo
               enddo
            enddo
         endif 
     else
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1(i,j,k) = eady_rate_zave(i,j)
               enddo
            enddo
         enddo
     endif


     ! diffusivity computed as coeff*(N/Nref)**2
     if(agm_closure_n2_scale) then 

         if(agm_closure_n2_scale_nref_cst) then 
             wrk2_2d(:,:) = rho0*grav_r*agm_closure_buoy_freq**2
         else 
             ! reference drhodz taken one level beneath blayer base,
             ! and no deeper than one cell from bottom.
             do j=jsc,jec
                do i=isc,iec
                   if(ksurf_blayer(i,j) > 0 .and. Grd%kmt(i,j) > 1) then 
                       k = min(Grd%kmt(i,j)-1, ksurf_blayer(i,j)+1)
                       wrk2_2d(i,j) = abs(drhodz_zt(i,j,k))
                   endif
                enddo
             enddo
         endif

         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  agm_growth_rate(i,j,k) = Grd%tmask(i,j,k)*wrk1(i,j,k)   ! computed just for diagnostics 
                  wrk3(i,j,k) = agm_closure_n2_scale_coeff*Grd%tmask(i,j,k) &
                                *abs(drhodz_zt(i,j,k))/(epsln+wrk2_2d(i,j))
                  wrk4(i,j,k) = agm_grid_scaling(i,j)*(wrk3(i,j,k) + agm_micom(i,j)) 
                  wrk4(i,j,k) = Grd%tmask(i,j,k)*max(agm_closure_min, min(agm_closure_max, wrk4(i,j,k)))      
               enddo
            enddo
         enddo


     ! diffusivity computed as growth_rate*length_scale**2
     else 

         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  agm_growth_rate(i,j,k) = Grd%tmask(i,j,k)*wrk1(i,j,k)
                  wrk3(i,j,k) = agm_closure_scaling*agm_growth_rate(i,j,k)*agm_length(i,j,k)**2 
                  wrk4(i,j,k) = agm_grid_scaling(i,j)*(wrk3(i,j,k) + agm_micom(i,j)) 
                  wrk4(i,j,k) = Grd%tmask(i,j,k)*max(agm_closure_min, min(agm_closure_max, wrk4(i,j,k)))      
               enddo
            enddo
         enddo

     endif


     ! time damping to get slowly evolving diffusivity 
     if(agm_smooth_time) then 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  agm_array(i,j,k) = &
                  Grd%tmask(i,j,k)*(agm_array(i,j,k) - gamma_damp*(agm_array(i,j,k)-wrk4(i,j,k)))
               enddo
            enddo
         enddo
     else
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  agm_array(i,j,k) = wrk4(i,j,k)            
               enddo
            enddo
         enddo
     endif
  

     ! spatial smoothing 
     if(agm_smooth_space) then 
         wrk1_2d(:,:) = 0.0
         call mpp_update_domains (agm_array(:,:,:), Dom%domain2d) 
         do k=1,nk

            do j=jsc,jec
               do i=isc,iec
                  if(Grd%tmask(i,j,k)==1.0) then 
                      active_cells = 4.0    +&
                        Grd%tmask(i-1,j,k)  +&
                        Grd%tmask(i+1,j,k)  +&
                        Grd%tmask(i,j-1,k)  +&
                        Grd%tmask(i,j+1,k)
                      if (active_cells > 4.0) then
                          wrk1_2d(i,j) = & 
                            (4.0*agm_array(i,j,k) +&
                            agm_array(i-1,j,k)    +&
                            agm_array(i+1,j,k)    +&
                            agm_array(i,j-1,k)    +&
                            agm_array(i,j+1,k)) / active_cells
                      else
                          wrk1_2d(i,j) = agm_array(i,j,k)
                      endif
                  endif
               enddo
            enddo
            do j=jsc,jec
               do i=isc,iec
                  agm_array(i,j,k) = wrk1_2d(i,j)*Grd%tmask(i,j,k)
               enddo
            enddo

         enddo
     endif

     ! need agm_array on full data domain 
     call mpp_update_domains (agm_array(:,:,:), Dom%domain2d) 


     if (id_agm_growth_rate > 0) used = send_data (id_agm_growth_rate, agm_growth_rate(:,:,:),&
                                 model_time, rmask=Grd%tmask(:,:,:),                          & 
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_agm_length > 0) used = send_data (id_agm_length, agm_length(:,:,:), &
                            model_time, rmask=Grd%tmask(:,:,:),                 &
                            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_rhines_length > 0) used = send_data (id_rhines_length, wrk2(:,:,:), &
                               model_time, rmask=Grd%tmask(:,:,:),              & 
                               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_agm_qg > 0)  used = send_data (id_agm_qg, wrk3(:,:,:),&
                         model_time, rmask=Grd%tmask(:,:,:),      &
                         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_growth_rate_baroclinic > 0) used = send_data (id_growth_rate_baroclinic, wrk1(:,:,:), &
                                        model_time, rmask=Grd%tmask(:,:,:),                       &
                                        is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if(id_N2_nblayer_base > 0) then 
         do j=jsd,jed
            do i=isd,ied
               wrk3_2d(i,j) = grav*rho0r*wrk2_2d(i,j)
            enddo
         enddo
         used = send_data (id_N2_nblayer_base, wrk2_2d(:,:), &
         model_time, rmask=Grd%tmask(:,:,1),                 &
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec) 
     endif

     if(id_N2_for_agm > 0) then 
         wrk4(:,:,:) = 0.0
         do k=1,nk
            do j=jsd,jed
               do i=isd,ied
                  wrk4(i,j,k) = grav*rho0r*abs(drhodz_zt(i,j,k))
               enddo
            enddo
         enddo
         used = send_data (id_N2_for_agm, wrk4(:,:,:), &
         model_time, rmask=Grd%tmask(:,:,:),           &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk) 
     endif


  endif ! endif of agm_closure 


  ! set redi diffusivity 
  if(aredi_diffusivity_grid_scaling) then 
      do k=1,nk
         aredi_array(:,:,k) = aredi*aredi_grid_scaling(:,:)*Grd%tmask(:,:,k) 
      enddo
  endif
  if(aredi_equal_agm) then 
      aredi_array(:,:,:) = agm_array(:,:,:) 
  endif


  ! diagnostics 
  if (id_agm > 0) used = send_data (id_agm, agm_array(:,:,1),&
                         model_time, rmask=Grd%tmask(:,:,1), &
                         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec) 

  if (id_agm_3d > 0) used = send_data (id_agm_3d, agm_array(:,:,:),&
                            model_time, rmask=Grd%tmask(:,:,:),    &
                            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk) 

  if (id_aredi > 0) used = send_data (id_aredi, aredi_array(:,:,1),&
                           model_time, rmask=Grd%tmask(:,:,1),     &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_aredi_3d > 0) used = send_data (id_aredi_3d, aredi_array(:,:,:),&
                              model_time, rmask=Grd%tmask(:,:,:),        &
                              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  if (id_ksurf_blayer >  0) then 
     wrk1_2d(:,:) = 0.0
     do j=jsd,jed
        do i=isd,ied
           wrk1_2d(i,j) = ksurf_blayer(i,j)
        enddo
     enddo
     used = send_data (id_ksurf_blayer, wrk1_2d(:,:), &
              model_time, rmask=Grd%tmask(:,:,1),     &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  if(debug_this_module) then 
      write(stdoutunit,*)' ' 
      write(stdoutunit,*) 'From ocean_nphysics_util'
      call write_timestamp(model_time)
      write (stdoutunit, *) &
           'checksum agm_array    ', &
          mpp_chksum(agm_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))
  endif


  call mpp_clock_end(id_clock_compute_diffusivity)

end subroutine compute_diffusivity
! </SUBROUTINE> NAME="compute_diffusivity"


!#######################################################################
! <SUBROUTINE NAME="transport_on_nrho_gm">
!
! <DESCRIPTION>
! Classify horizontal GM mass transport according to neutral density classes. 
!
! NOTE: This diagnostic works with transport integrated from bottom to 
! a particular cell depth. To get transport_on_nrho_gm, a remapping is 
! performed, rather than the binning done for trans_rho.  
!
! Code history 
! 2008: algorithm based (incorrectly) on transport_on_rho 
! 2009: algorithm corrected to be consistent with remapping 
!       used in tracer_on_rho algorithm
! </DESCRIPTION>
!
subroutine transport_on_nrho_gm (Time, Dens, tx_trans_gm, ty_trans_gm)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: tx_trans_gm
  real, dimension(isd:,jsd:,:), intent(in) :: ty_trans_gm
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_rho, neutralrho_nk
  real    :: work(isd:ied,jsd:jed,size(Dens%neutralrho_ref),2)
  real    :: W1, W2

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (transport_on_nrho_gm): module needs initialization ')
  endif 

  next_time = increment_time(Time%model_time, int(dtime), 0)

  if (need_data(id_tx_trans_nrho_gm,next_time) .or. need_data(id_ty_trans_nrho_gm,next_time)) then

      neutralrho_nk = size(Dens%neutralrho_ref(:))
      work(:,:,:,:) = 0.0

      ! for (i,j) points with neutralrho_ref < neutralrho(k=1),   work=0
      ! for (i,j) points with neutralrho_ref > neutralrho(k=kmt), work=0
      ! these assumptions mean there is no need to specially handle the endpoints,
      ! since the initial value for work is 0.

      ! interpolate trans_gm from k-levels to neutralrho_nk-levels
      do k_rho=1,neutralrho_nk
         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%neutralrho_ref(k_rho) >  Dens%neutralrho(i,j,k)  ) then
                      if( Dens%neutralrho_ref(k_rho) <= Dens%neutralrho(i,j,k+1)) then 
                          W1= Dens%neutralrho_ref(k_rho)- Dens%neutralrho(i,j,k)
                          W2= Dens%neutralrho(i,j,k+1)  - Dens%neutralrho_ref(k_rho)
                          work(i,j,k_rho,1) = (tx_trans_gm(i,j,k+1)*W1 +tx_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                          work(i,j,k_rho,2) = (ty_trans_gm(i,j,k+1)*W1 +ty_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

      do k_rho=1,neutralrho_nk
         do j=jsc,jec
            do i=isc,iec
               work(i,j,k_rho,1) = work(i,j,k_rho,1)*Grd%tmask(i,j,1)
               work(i,j,k_rho,2) = work(i,j,k_rho,2)*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_nrho_gm > 0) then 
          used = send_data (id_tx_trans_nrho_gm, work(:,:,:,1), Time%model_time, &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif
      if (id_ty_trans_nrho_gm > 0) then 
          used = send_data (id_ty_trans_nrho_gm, work(:,:,:,2), Time%model_time, &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=neutralrho_nk)
      endif

  endif

end subroutine transport_on_nrho_gm
! </SUBROUTINE> NAME="transport_on_nrho_gm"


!#######################################################################
! <SUBROUTINE NAME="transport_on_rho_gm">
!
! <DESCRIPTION>
! Classify horizontal GM mass transport according to potential density classes. 
!
! Algorithm based on linear interpolation of function on s-surfaces to 
! function on rho-surfaces.  
!
! Diagnostic makes sense when potrho is monotonically increasing with 
! depth, although the algorithm does not explicitly make this assumption.  
!
! NOTE: This diagnostic works with transport integrated from bottom to 
! a particular cell depth. To get transport_on_rho_gm, a remapping is 
! performed, rather than the binning done for trans_rho.  
!
! Code history 
!
! 2008: algorithm based (incorrectly) on transport_on_rho 
! 2009: algorithm corrected to be consistent with remapping 
!       used in tracer_on_rho algorithm
!
! </DESCRIPTION>
!
subroutine transport_on_rho_gm (Time, Dens, tx_trans_gm, ty_trans_gm)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  real, dimension(isd:,jsd:,:), intent(in) :: tx_trans_gm
  real, dimension(isd:,jsd:,:), intent(in) :: ty_trans_gm
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_rho, potrho_nk
  real    :: work(isd:ied,jsd:jed,size(Dens%potrho_ref),2)
  real    :: W1, W2

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (transport_on_rho_gm): module needs initialization ')
  endif 

  next_time = increment_time(Time%model_time, int(dtime), 0)
  if (need_data(id_tx_trans_rho_gm,next_time) .or. need_data(id_ty_trans_rho_gm,next_time)) then

      potrho_nk = size(Dens%potrho_ref(:))
      work(:,:,:,:) = 0.0

      ! for (i,j) points with potrho_ref < potrho(k=1),   work=0
      ! for (i,j) points with potrho_ref > potrho(k=kmt), work=0
      ! these assumptions mean there is no need to specially handle the endpoints,
      ! since the initial value for work is 0.

      ! interpolate trans_gm from k-levels to krho-levels
      do k_rho=1,potrho_nk
         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%potrho_ref(k_rho) >  Dens%potrho(i,j,k)  ) then
                      if( Dens%potrho_ref(k_rho) <= Dens%potrho(i,j,k+1)) then 
                          W1= Dens%potrho_ref(k_rho)- Dens%potrho(i,j,k)
                          W2= Dens%potrho(i,j,k+1)  - Dens%potrho_ref(k_rho)
                          work(i,j,k_rho,1) = (tx_trans_gm(i,j,k+1)*W1 +tx_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                          work(i,j,k_rho,2) = (ty_trans_gm(i,j,k+1)*W1 +ty_trans_gm(i,j,k)*W2) &
                                              /(W1 + W2 + epsln)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

      do k_rho=1,potrho_nk
         do j=jsc,jec
            do i=isc,iec
               work(i,j,k_rho,1) = work(i,j,k_rho,1)*Grd%tmask(i,j,1)
               work(i,j,k_rho,2) = work(i,j,k_rho,2)*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_rho_gm > 0) then 
          used = send_data (id_tx_trans_rho_gm, work(:,:,:,1), Time%model_time, &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)
      endif
      if (id_ty_trans_rho_gm > 0) then 
          used = send_data (id_ty_trans_rho_gm, work(:,:,:,2), Time%model_time, &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)
      endif

  endif


end subroutine transport_on_rho_gm
! </SUBROUTINE> NAME="transport_on_rho_gm"


!#######################################################################
! <SUBROUTINE NAME="transport_on_theta_gm">
!
! <DESCRIPTION>
! Classify horizontal GM mass transport according to potential temp classes. 
!
! Algorithm based on linear interpolation of function on s-surfaces to 
! function on rho-surfaces.  
!
! Diagnostic makes sense when potential temp is monotonically increasing 
! with depth, although the algorithm does not explicitly make this assumption.  
!
! NOTE: This diagnostic works with transport integrated from bottom to 
! a particular cell depth. To get transport_on_theta_gm, a remapping is 
! performed, rather than the binning done for trans_rho.  
!
! Code history 
!
! 2009: algorithm based on transport_on_rho_gm
!
! </DESCRIPTION>
!
subroutine transport_on_theta_gm (Time, Dens, Theta, tx_trans_gm, ty_trans_gm)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_prog_tracer_type), intent(in) :: Theta
  real, dimension(isd:,jsd:,:), intent(in) :: tx_trans_gm
  real, dimension(isd:,jsd:,:), intent(in) :: ty_trans_gm
  type(time_type)                          :: next_time 

  integer :: i, j, k, k_theta, theta_nk, tau
  real    :: work(isd:ied,jsd:jed,size(Dens%potrho_ref),2)
  real    :: W1, W2

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (transport_on_theta_gm): module needs initialization ')
  endif 

  tau = Time%tau
  next_time = increment_time(Time%model_time, int(dtime), 0)
  if (need_data(id_tx_trans_theta_gm,next_time) .or. need_data(id_ty_trans_theta_gm,next_time)) then

      theta_nk = size(Dens%theta_ref(:))
      work(:,:,:,:) = 0.0

      ! for (i,j) points with theta_ref > theta(k=1),   work=0
      ! for (i,j) points with theta_ref < theta(k=kmt), work=0
      ! these assumptions mean there is no need to specially handle the endpoints,
      ! since the initial value for work is 0.

      ! interpolate trans_gm from k-levels to theta-levels
      ! note sign change in the if-tests relative to transport on rho and nrho
      do k_theta=1,theta_nk
         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%theta_ref(k_theta) <  Theta%field(i,j,k,tau) ) then
                      if( Dens%theta_ref(k_theta) >= Theta%field(i,j,k+1,tau)) then 
                          W1= -Dens%theta_ref(k_theta)  + Theta%field(i,j,k,tau)
                          W2= -Theta%field(i,j,k+1,tau) + Dens%theta_ref(k_theta)
                          work(i,j,k_theta,1) = (tx_trans_gm(i,j,k+1)*W1 +tx_trans_gm(i,j,k)*W2) &
                                                /(W1 + W2 + epsln)
                          work(i,j,k_theta,2) = (ty_trans_gm(i,j,k+1)*W1 +ty_trans_gm(i,j,k)*W2) &
                                                /(W1 + W2 + epsln)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

      do k_theta=1,theta_nk
         do j=jsc,jec
            do i=isc,iec
               work(i,j,k_theta,1) = work(i,j,k_theta,1)*Grd%tmask(i,j,1)
               work(i,j,k_theta,2) = work(i,j,k_theta,2)*Grd%tmask(i,j,1)
            enddo
         enddo
      enddo

      if (id_tx_trans_theta_gm > 0) then 
          used = send_data (id_tx_trans_theta_gm, work(:,:,:,1), Time%model_time, &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=theta_nk)
      endif
      if (id_ty_trans_theta_gm > 0) then 
          used = send_data (id_ty_trans_theta_gm, work(:,:,:,2), Time%model_time, &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=theta_nk)
      endif

  endif


end subroutine transport_on_theta_gm
! </SUBROUTINE> NAME="transport_on_theta_gm"



!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_util_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_nphysics_util_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  call save_restart(Nphysics_util_restart, time_stamp)

end subroutine ocean_nphysics_util_restart
! </SUBROUTINE> NAME="ocean_nphysics_util_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysics_coeff_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysics_coeff_end(Time, agm_array, aredi_array, &
                    rossby_radius, rossby_radius_raw, bczone_radius)

  type(ocean_time_type),        intent(in) :: Time
  real, dimension(isd:,jsd:,:), intent(in) :: agm_array
  real, dimension(isd:,jsd:,:), intent(in) :: aredi_array
  real, dimension(isd:,jsd:),   intent(in) :: rossby_radius
  real, dimension(isd:,jsd:),   intent(in) :: rossby_radius_raw
  real, dimension(isd:,jsd:),   intent(in) :: bczone_radius 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  write(stdoutunit,*)' ' 
  write(stdoutunit,*) 'From ocean_nphysics_coeff_end: ending chksum'
  call write_timestamp(Time%model_time)

  write (stdoutunit, *) &
  'checksum ending agm_array',     mpp_chksum(agm_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))
  write (stdoutunit, *) &
  'checksum ending aredi_array',   mpp_chksum(aredi_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))
  write (stdoutunit, *) &
  'checksum ending rossby_radius', mpp_chksum(rossby_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
  write (stdoutunit, *) &
  'checksum ending rossby_radius_raw', mpp_chksum(rossby_radius_raw(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
  write (stdoutunit, *) &
  'checksum ending bczone_radius', mpp_chksum(bczone_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))

end subroutine ocean_nphysics_coeff_end
! </SUBROUTINE> NAME="ocean_nphysics_coeff_end"



end module ocean_nphysics_util_mod
