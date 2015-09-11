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
module ocean_nphysicsC_mod
! 
!<CONTACT EMAIL="stephen.griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted and density weighted time tendency for tracer 
! from Laplacian neutral diffusion + Laplacian skew-diffusion.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the cell thickness weighted and density 
! weighted tracer tendency from small angle Laplacian neutral diffusion
! plus Laplacian skew-diffusion.  The algorithms for neutral diffusion
! are based on mom4p0d methods.  The algorithm for neutral skewsion 
! are based on a projection onto a few of the lowest baroclinic 
! modes. This module is experimental, and should be used with caution. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, A. Gnanadesikan, R.C. Pacanowski, V. Larichev, 
! J.K. Dukowicz,  and R.D. Smith
! Isoneutral diffusion in a z-coordinate ocean model
! Journal of Physical Oceanography (1998) vol 28 pages 805-830
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! The Gent-McWilliams Skew-flux 
! Journal of Physical Oceanography (1998) vol 28 pages 831-841
! </REFERENCE>
!
! <REFERENCE>
! R. Ferrari, S.M. Griffies, A.J.G. Nurser, and G.K. Vallis  
! A boundary value problem for the parameterized mesoscale eddy transport
! Ocean Modelling, 2009.  
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! Fundamentals of Ocean Climate Models (2004)
! Princeton University Press 
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! Elements of MOM4p1 (2008)
! </REFERENCE>
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
! Gerdes, Koberle, and Willebrand
! The influence of numerical advection schemes on the results of ocean
! general circulation models, Climate Dynamics (1991), vol. 5, 
! pages 211--226. 
! </REFERENCE>
!
! <NOTE>
! Numerical implementation of the flux components follows the triad 
! approach documented in the references and implemented in MOM2 and MOM3.  
! The MOM4 algorithm accounts for partial bottom cells and generalized
! orthogonal horizontal coordinates.
! </NOTE> 
!
! <NOTE> 
! In steep neutral slope regions, neutral diffusive fluxes are tapered
! to zero with the tanh taper of Danabasoglu and McWilliams (1995) or the 
! quadratic scheme of Gerdes, Koberle, and Willebrand.  Tapering is 
! not required for the modal decomposed skew fluxes.  
! </NOTE> 
!
! </INFO>
!
!<NAMELIST NAME="ocean_nphysicsC_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  </DATA> 
!
!  <DATA NAME="epsln_bv_freq" UNITS="kg/m4" TYPE="real">
!  Minimum buoyancy frequency accepted for the computation of 
!  baroclinic modes. Default epsln_bv_freq=1e-10.  Note there 
!  is also a minimum drhodz set in ocean_density.F90 via the
!  nml epsln_drhodz in that module.  We provide yet another minimum
!  here in case we need an extra regularization for the amplitude
!  of the baroclinic modes.  
!  </DATA>
!
!  <DATA NAME="do_neutral_diffusion" TYPE="logical">
!  To compute tendency from neutral diffusion.
!  Default do_neutral_diffusion=.true. 
!  </DATA> 
!  <DATA NAME="do_gm_skewsion" TYPE="logical">
!  To compute tendency from GM skewsion. Default do_gm_skewsion=.true. 
!  </DATA> 
!  <DATA NAME="gm_skewsion_modes" TYPE="logical">
!  To compute tendency from GM skewsion using streamfunction established
!  by baroclinic modes. Default gm_skewsion_modes=.false.
!  </DATA> 
!  <DATA NAME="gm_skewsion_bvproblem" TYPE="logical">
!  To compute tendency from GM skewsion using streamfunction established
!  by a boundary value problem. Default gm_skewsion_bvproblem=.false.
!  </DATA> 
!
!  <DATA NAME="number_bc_modes" TYPE="integer">
!  The number of baroclinic modes used to construct the eddy induced 
!  streamfunction. Default number_bc_modes=1.
!  </DATA> 
!  <DATA NAME="bvp_bc_mode" TYPE="integer">
!  The particular baroclinic mode used to construct the BVP streamfunction.
!  If bvp_bc_mode=0, then will set bc_speed=0 when computing the BVP streamfunction.
!  Default bvp_bc_mode=1.  
!  </DATA> 
!  <DATA NAME="bvp_constant_speed" TYPE="logical">
!  For taking a constant speed to be used for the calculation 
!  of the BVP streamfunction. Default bvp_constant_speed=.false.  
!  </DATA> 
!  <DATA NAME="bvp_speed" UNITS="m/s" TYPE="real">
!  For setting the speed weighting the second order derivative operator 
!  in the BVP streamfunction method: 
!  c^2 = max[bvp_min_speed, (bvp_speed-c_mode)^2].
!  If bvp_constant_speed, then  c^2 = bvp_speed^2. 
!  Default bvp_speed=0.0, in which case c^2 = c_mode^2.  
!  </DATA> 
!  <DATA NAME="bvp_min_speed" UNITS="m/s" TYPE="real">
!  For setting a minimum speed for use with the calculation 
!  of the BVP streamfunction. We need  bvp_min_speed>0 to ensure
!  that the second order derivative operator contributes to the 
!  calculation of the streamfunction.  
!  Default bvp_min_speed=0.1.  
!  </DATA> 
! 
!  <DATA NAME="bv_freq_smooth_vert" TYPE="logical">
!  To smooth the buoyancy frequency for use in 
!  computing the baroclinic modes. Generally this field has already 
!  been smooted in ocean_density_mod, but we maintain the possibility of 
!  further smoothing here.  Default bv_freq_smooth_vert=.false.
!  </DATA>
!  <DATA NAME="num_121_passes" TYPE="integer">
!  The number of 121 passes used to smooth buoyancy frequency when 
!  bv_freq_smooth_vert=.true.  Default num_121_passes=1.
!  </DATA>
!  <DATA NAME="min_bc_speed"  UNITS="m/s"  TYPE="real">
!  The minimum speed used for computing the baroclinic modes. 
!  Default min_bc_speed=1e-6
!  </DATA> 
!
!  <DATA NAME="smooth_bc_modes" TYPE="logical">
!  For doing a vertical 1-2-1 smoothing on the baroclinic modes
!  prior to normalization.  This is useful to reduce noise.
!  Default smooth_bc_modes=.false.
!  </DATA> 
!
!  <DATA NAME="smooth_psi" TYPE="logical">
!  For doing a horizontal 1-2-1 smoothing on the psix and psiy fields. 
!  This is useful to reduce noise. Default smooth_psi=.true.
!  </DATA> 
!
!  <DATA NAME="regularize_psi" TYPE="logical">
!  To reduce the magnitude of psi in regions of weak stratification, 
!  using the slope = smax_psi to set the overall scale of the max allowed
!  for psi. Default regularize_psi=.true.
!  </DATA> 
!  <DATA NAME="smax_modes" TYPE="real">
!  Maximum slope used for setting the overall scale of a modal 
!  contribution to the parameterized transport.   
!  Default smax_psi=0.1.  
!  </DATA> 
! 
!  <DATA NAME="diffusion_all_explicit" TYPE="logical">
!  To compute all contributions from neutral diffusion explicitly in time, including
!  the K33 diagonal piece.  This approach is available only when have small time 
!  steps and/or running with just a single tracer.  It is for testing purposes. 
!  </DATA> 
!
!  <DATA NAME="neutral_physics_limit" TYPE="logical">
!  When tracer falls outside a specified range, revert to horizontal 
!  diffusive fluxes at this cell. This is an ad hoc and incomplete attempt
!  to maintain monotonicity with the neutral physics scheme.  
!  Default neutral_physics_limit=.true.
!  </DATA> 
!  <DATA NAME="tmask_neutral_on" TYPE="logical">
!  If .true. then this logical reduces the neutral diffusive fluxes to 
!  horizontal/vertical diffusion next to boundaries.  
!  This approach has been found to reduce spurious 
!  extrema resulting from truncation of triads used to compute 
!  a neutral flux component.   
!  Default tmask_neutral_on=.false.
!  </DATA> 
!
!  <DATA NAME="dm_taper" TYPE="logical">
!  Set to true to use the tanh tapering scheme of Danabasoglu and McWilliams.
!  Default is true. 
!  </DATA> 
!  <DATA NAME="gkw_taper" TYPE="logical">
!  Set to true to use the quadradic tapering scheme of Gerdes, Koberle, and Willebrand.
!  Default is false. 
!  </DATA> 
!
!  <DATA NAME="neutral_eddy_depth" TYPE="logical">
!  Compute eddy_depth according to depth over which eddies feel the ocean surface.
!  Default neutral_eddy_depth=.true. 
!  </DATA> 
!
!  <DATA NAME="turb_blayer_min" TYPE="real">
!  Minimum depth of a surface turbulent boundary layer
!  used in the transition of the neutral diffusion fluxes
!  to the surface.  Note that in mom4p0, 
!  turb_blayer_min was always set to zero. 
!  </DATA> 
!
!  <DATA NAME="transport_units" TYPE="character">
!  The units for writing out the transport.  Either in 
!  Sv (10^9 kg/s) or mks (kg/s). Note the mks unit is requested 
!  for CMIP5 purposes.
!  Default transport_units = 'Sv'. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,           only: epsln, pi, grav
use diag_manager_mod,        only: register_diag_field, register_static_field, send_data, need_data
use fms_mod,                 only: FATAL, WARNING, NOTE
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_io_mod,              only: string 
use mpp_domains_mod,         only: mpp_update_domains
use mpp_domains_mod,         only: CGRID_NE, WUPDATE, SUPDATE, EUPDATE, NUPDATE
use mpp_domains_mod,         only: cyclic_global_domain, global_data_domain 
use mpp_mod,                 only: mpp_error, mpp_chksum, stdout, stdlog
use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,        only: set_time, time_type, increment_time, operator ( + )

use ocean_domains_mod,         only: get_local_indices, set_ocean_domain
use ocean_nphysics_util_mod,   only: ocean_nphysics_coeff_init, ocean_nphysics_coeff_end
use ocean_nphysics_util_mod,   only: neutral_slopes, tracer_derivs 
use ocean_nphysics_util_mod,   only: compute_eady_rate, compute_baroclinicity
use ocean_nphysics_util_mod,   only: compute_rossby_radius, compute_bczone_radius
use ocean_nphysics_util_mod,   only: compute_diffusivity, ocean_nphysics_util_restart
use ocean_nphysics_util_mod,   only: transport_on_nrho_gm, transport_on_rho_gm, transport_on_theta_gm
use ocean_operators_mod,       only: FAX, FAY, FMX, FMY, BDX_ET, BDY_NT
use ocean_parameters_mod,      only: missing_value, onehalf, onefourth, oneeigth, DEPTH_BASED
use ocean_parameters_mod,      only: rho0r, rho0
use ocean_types_mod,           only: ocean_grid_type, ocean_domain_type, ocean_density_type
use ocean_types_mod,           only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,           only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,           only: tracer_2d_type, tracer_3d_0_nk_type, tracer_3d_1_nk_type 
use ocean_util_mod,            only: write_timestamp
use ocean_workspace_mod,       only: wrk1_2d, wrk1_v2d, wrk2_v2d
use ocean_workspace_mod,       only: wrk1, wrk2, wrk3, wrk4
use ocean_workspace_mod,       only: wrk1_v, wrk2_v, wrk3_v, wrk4_v

implicit none

public ocean_nphysicsC_init
public ocean_nphysicsC_end
public nphysicsC
public ocean_nphysicsC_restart

private fx_flux_ndiffuse
private fy_flux_ndiffuse
private fz_flux_ndiffuse
private fx_flux_gm
private fy_flux_gm
private fz_flux_gm
private fz_terms
private neutral_blayer
private compute_ndiffusion
private compute_gmskewsion
private baroclinic_modes 
private compute_psi_modes 
private compute_psi_bvp
private invtri_bvp

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

integer :: num_prog_tracers = 0

! clock ids
integer :: id_clock_ndiffuse
integer :: id_clock_gm_modes
integer :: id_clock_neutral_blayer
integer :: id_clock_fz_terms 
integer :: id_clock_fx_flux_ndiffuse
integer :: id_clock_fy_flux_ndiffuse
integer :: id_clock_fz_flux_ndiffuse
integer :: id_clock_fx_flux_gm
integer :: id_clock_fy_flux_gm
integer :: id_clock_fz_flux_gm
integer :: id_clock_bc_modes 
integer :: id_clock_psi_modes
integer :: id_clock_psi_bvp

! diagnostic manager ids
logical :: used
integer :: id_k33_explicit     =-1
integer :: id_N_squared        =-1
integer :: id_bv_freq          =-1
integer :: id_tx_trans_gm      =-1
integer :: id_ty_trans_gm      =-1
integer :: id_eddy_depth       =-1
integer :: id_psix_gm_modes    =-1
integer :: id_psiy_gm_modes    =-1
integer :: id_psix_gm_bvp      =-1
integer :: id_psiy_gm_bvp      =-1
integer :: id_gradx_rho        =-1
integer :: id_grady_rho        =-1
integer :: id_rescale_psi_x    =-1
integer :: id_rescale_psi_y    =-1
integer :: id_gm_eddy_ke_source=-1

integer                             :: id_neutral_rho_ndiffuse  ! tendency of neutral density from neutral diffuse
integer                             :: id_wdian_rho_ndiffuse    ! contribution to wdian from neutral diffuse
integer                             :: id_neutral_rho_gm_modes  ! tendency of neutral density from GM-modes 
integer                             :: id_wdian_rho_gm_modes    ! contribution to wdian from GM-modes
integer, dimension(:), allocatable  :: id_neutral_diffuse       ! tendency from neutral diffusion 
integer, dimension(:), allocatable  :: id_neutral_gm_modes      ! tendency from GM as projected into modes 
integer, dimension(:), allocatable  :: id_k33_implicit          ! K33 handled implicitly in time 
integer, dimension(:), allocatable  :: id_flux_x_ndiffuse       ! i-directed heat flux from neutral diffuse
integer, dimension(:), allocatable  :: id_flux_y_ndiffuse       ! j-directed heat flux from neutral diffuse 
integer, dimension(:), allocatable  :: id_flux_x_gm_modes       ! i-directed heat flux from neutral gm modes
integer, dimension(:), allocatable  :: id_flux_y_gm_modes       ! j-directed heat flux from neutral gm modes 
integer, dimension(:), allocatable  :: id_flux_x_ndiffuse_int_z ! vertically integrated i-directed diffuse flux 
integer, dimension(:), allocatable  :: id_flux_y_ndiffuse_int_z ! vertically integrated j-directed diffuse flux 
integer, dimension(:), allocatable  :: id_flux_x_gm_modes_int_z ! vertically integrated i-directed gm modes flux 
integer, dimension(:), allocatable  :: id_flux_y_gm_modes_int_z ! vertically integrated j-directed gm modes flux 

integer, dimension(:), allocatable  :: id_bc_mode  ! baroclinic modes
integer, dimension(:), allocatable  :: id_bc_speed ! baroclinic wave speeds 

#include <ocean_memory.h>

#ifdef MOM4_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed,nk,0:1) :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(isd:ied,jsd:jed,0:nk)   :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(isd:ied,jsd:jed,0:1)    :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(isd:ied,jsd:jed,nk) :: aredi_array !3D array of redi diffusivities (m^2/sec)     
real, dimension(isd:ied,jsd:jed,nk) :: agm_array   !3D array of gm diffusivities (m^2/sec)        
real, dimension(isd:ied,jsd:jed)    :: ah_array    !2D array of micom horizontal diffusivities (not used)

real, dimension(isd:ied,jsd:jed,nk,0:1) :: psix ! streamfunction x-component (m^2/sec) 
real, dimension(isd:ied,jsd:jed,nk,0:1) :: psiy ! streamfunction y-component (m^2/sec) 

real, dimension(isd:ied,jsd:jed,nk) :: bv_freq ! buoyancy frequency (1/sec) 

real, dimension(isd:ied,jsd:jed)    :: bczone_radius     !for bczone calculation (m) 
real, dimension(isd:ied,jsd:jed)    :: rossby_radius     !first baroclinic Rossby radius (m) 
real, dimension(isd:ied,jsd:jed)    :: rossby_radius_raw !first baroclinic Rossby radius (m) without max/min 
real, dimension(isd:ied,jsd:jed,nk) :: eady_termx        !rho_z*(S_x)^2 for computing Eady growth rate 
real, dimension(isd:ied,jsd:jed,nk) :: eady_termy        !rho_z*(S_y)^2 for computing Eady growth rate 
real, dimension(isd:ied,jsd:jed)    :: baroclinic_termx  !intermediate term for computing vert ave baroclinicity
real, dimension(isd:ied,jsd:jed)    :: baroclinic_termy  !intermediate term for computing vert ave baroclinicity 
real, dimension(isd:ied,jsd:jed)    :: grid_length       !grid length scale (m)

real, dimension(isd:ied,jsd:jed,nk)     :: drhodT       !drho/dtheta     (kg/m^3/C)
real, dimension(isd:ied,jsd:jed,nk)     :: drhodS       !drho/dsalinity  (kg/m^3/psu)
real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodzb      !vertical neutral density gradient (kg/m^3/m)
real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodzh      !vertical neutral density gradient (kg/m^3/m)

real, dimension(isd:ied,jsd:jed,nk)     :: K33_implicit !density weighted (kg/m^3) implicit in time 
                                                        !diagonal term in redi tensor (m^2/sec) 
real, dimension(isd:ied,jsd:jed,nk)     :: K33_explicit !density weighted (kg/m^3) explicit in time 
                                                        !diagonal term in redi tensor (m^2/sec) 

real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_31 !tracer independent portion of mixing tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: tensor_32 !tracer independent portion of mixing tensor

integer, dimension(isd:ied,jsd:jed) :: ksurf_blayer  ! k-value at base of surface nblayer 
real, dimension(isd:ied,jsd:jed)    :: eddy_depth    !max of depth(m) mesoscale eddies penetrate & kpp bldepth
real, dimension(isd:ied,jsd:jed,nk) :: tx_trans_gm   !diagnosing i-transport due to GM 
real, dimension(isd:ied,jsd:jed,nk) :: ty_trans_gm   !diagnosing j-transport due to GM 

#else

real, dimension(:,:,:,:),   allocatable :: delqc   !density weighted (kg/m^3) quarter cell thickness(m)
real, dimension(:,:,:),     allocatable :: dzwtr   !(1/dzwt)(m^-1)
real, dimension(:,:,:),     allocatable :: dtew    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),     allocatable :: dtns    !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),     allocatable :: dtwedyt !horizontal areas (m^2) of quarter cell
real, dimension(:,:,:),     allocatable :: dxtdtsn !horizontal areas (m^2) of quarter cell

real, dimension(:,:,:),     allocatable :: aredi_array !3D array of redi diffusivities (m^2/sec) 
real, dimension(:,:,:),     allocatable :: agm_array   !3D array of gm diffusivities (m^2/sec)     
real, dimension(:,:),       allocatable :: ah_array    !2D array of micom horizontal diffusivities (not used)

real, dimension(:,:,:,:), allocatable :: psix ! streamfunction x-component (m^2/sec) 
real, dimension(:,:,:,:), allocatable :: psiy ! streamfunction y-component (m^2/sec) 

real, dimension(:,:,:), allocatable :: bv_freq ! buoyancy frequency (1/sec) 

real, dimension(:,:),       allocatable :: bczone_radius     !for bzcone calculation (m) 
real, dimension(:,:),       allocatable :: rossby_radius     !first baroclinic Rossby radius (m) 
real, dimension(:,:),       allocatable :: rossby_radius_raw !first baroclinic Rossby radius (m) without max/min
real, dimension(:,:,:),     allocatable :: eady_termx        !rho_z*(S_x)^2 for computing Eady growth rate 
real, dimension(:,:,:),     allocatable :: eady_termy        !rho_z*(S_y)^2 for computing Eady growth rate 
real, dimension(:,:),       allocatable :: baroclinic_termx  !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: baroclinic_termy  !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: grid_length       !grid length scale (m)

real, dimension(:,:,:),     allocatable :: drhodT   !drho/dtheta     (kg/m^3/C)
real, dimension(:,:,:),     allocatable :: drhodS   !drho/dsalinity  (kg/m^3/psu)
real, dimension(:,:,:,:),   allocatable :: drhodzb  !vertical neutral density gradient (kg/m^3/m)
real, dimension(:,:,:,:),   allocatable :: drhodzh  !vertical neutral density gradient (kg/m^3/m)

real, dimension(:,:,:,:,:), allocatable :: tensor_31 !tracer independent portion of mixing tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_32 !tracer independent portion of mixing tensor

real, dimension(:,:,:),     allocatable :: K33_implicit !density weighted (kg/m^3) implicit in time 
                                                        !diagonal term in redi tensor (m^2/sec) 
real, dimension(:,:,:),     allocatable :: K33_explicit !density weighted (kg/m^3) explicit in time 

integer, dimension(:,:), allocatable :: ksurf_blayer    ! k-value at base of surface nblayer 
real, dimension(:,:),    allocatable :: eddy_depth      !max of depth (m) mesoscale eddies penetrate & kpp bldepth

real, dimension(:,:,:),  allocatable :: tx_trans_gm     !for diagnosing i-transport due to GM 
real, dimension(:,:,:),  allocatable :: ty_trans_gm     !for diagnosing j-transport due to GM 

#endif

real, dimension(:,:,:,:), allocatable :: bc_mode   !dimensionless normalized baroclinic mode
real, dimension(:,:,:),   allocatable :: bc_speed  !wave speed (m/s) for baroclinic modes 
real, dimension(:,:),     allocatable :: bc_speed2 !squared wave speed (m/s) for chosen baroclinic mode 

! introduce following derived types so that do not need to know num_prog_tracers at compile time 
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdx    ! tracer partial derivative (tracer/m)
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdy    ! tracer partial derivative (tracer/m)
type(tracer_3d_0_nk_type), dimension(:), allocatable  :: dTdz    ! tracer partial derivative (tracer/m)
type(tracer_2d_type),      dimension(:), allocatable  :: fz1     ! z-flux component for tracers at particular k 
type(tracer_2d_type),      dimension(:), allocatable  :: fz2     ! z-flux component for tracers at particular k 
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_x  ! i-flux component for tracers for all k
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_y  ! j-flux component for tracers for all k


integer :: index_temp
integer :: index_salt

character(len=128) :: version=&
     '$Id: ocean_nphysicsC.F90,v 1.1.2.28.4.1.2.7.30.1.2.1.16.1.32.1.4.1 2009/10/21 20:07:28 smg Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: module_is_initialized = .FALSE.

! time step settings 
real    :: dtime
real    :: two_dtime_inv

! constants
real :: pi_r
real :: grav_rho0_r
real :: grav_rho0r
real :: sqrt_grav

! vertical coordinate 
integer :: vert_coordinate_class

! lower and upper depth for vertically averaging ocean properties.
! read into ocean_nphysics_util_nml
real :: agm_closure_upper_depth
real :: agm_closure_lower_depth


!**************nml settings**************

! for turning on/off the schemes 
logical :: do_neutral_diffusion  = .true.
logical :: do_gm_skewsion        = .true.
logical :: gm_skewsion_modes     = .false.
logical :: gm_skewsion_bvproblem = .false.

! number of modes used to construct GM streamfunction 
! when gm_skewsion_modes=.true.
integer :: number_bc_modes = 1

! the particular baroclinic mode used to construct the BVP 
! streamfunction.  
integer :: bvp_bc_mode=1

! for setting the speed to constant with the BVP calculation.
logical :: bvp_constant_speed=.false. 
real    :: bvp_speed         = 0.0
real    :: bvp_min_speed     = 0.1

! for regularizing the streamfunction psi when computed from modes
logical :: regularize_psi=.true.
real    :: smax_psi=0.1 

! for setting minimum buoyancy frequency for computing modes
real    :: epsln_bv_freq=1.e-10

! min baroclinic wave speed (m/s)
real    :: min_bc_speed = 1.e-6

! for smoothing psix and psiy
logical :: smooth_psi = .true.    

! for smoothing baroclinic modes prior to normalization 
logical :: smooth_bc_modes = .false.    

! for smoothing the buoyancy frequency 
logical :: bv_freq_smooth_vert=.false.
integer :: num_121_passes     =1

! for setting the slope tapering methods 
real    :: smax                        ! set in ocean_neutral_util_nml
real    :: swidth                      ! set in ocean_neutral_util_nml
logical :: dm_taper         = .true.   ! tanh tapering scheme of Danabasoglu and McWilliams
real    :: swidthr                     ! inverse swidth  
real    :: smax_swidthr                ! useful combination of terms 
real    :: dm_taper_const   = 1.0      ! internally set to unity when dm_taper=.true.
logical :: gkw_taper        = .false.  ! quadratic tapering of Gerdes, Koberle, and Willebrand
real    :: gkw_taper_const  = 0.0      ! internally set to unity when gkw_taper=.true.

! for neutral blayer 
real :: turb_blayer_min = 0.0 ! metres (was always set to zero in mom4p0)

! for maintaining horizontal GM velocity constant with depth within neutral boundary layer 
logical :: neutral_linear_gm_taper=.true.   

! to compute an eddy depth, where shallower eddies feel the surface
logical :: neutral_eddy_depth=.true.        

! to reduce neutral fluxes to horz/vert diffusion next to model boundaries
logical :: tmask_neutral_on=.false.         

! to compute K33 explicitly in time. 
! diffusion_all_explicit=.false. for realistic simulations.
logical :: diffusion_all_explicit=.false.   

! revert to horizontal diffusion when tracer falls outside specified range 
logical :: neutral_physics_limit=.true.   

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_units='Sv' 
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 

! for the module as a whole 
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.

!**************end of nml settings**************

namelist /ocean_nphysicsC_nml/ use_this_module, debug_this_module, &
          do_neutral_diffusion, do_gm_skewsion,                    &
          neutral_physics_limit, neutral_eddy_depth,               &
          dm_taper, gkw_taper, tmask_neutral_on,                   &
          diffusion_all_explicit, turb_blayer_min,                 &
          gm_skewsion_modes, number_bc_modes,                      &
          gm_skewsion_bvproblem, bvp_bc_mode,                      &
          bvp_min_speed, bvp_speed, bvp_constant_speed,            &
          bv_freq_smooth_vert, num_121_passes,                     &
          min_bc_speed, smooth_psi, epsln_bv_freq,                 &
          regularize_psi, smax_psi, smooth_bc_modes, transport_units

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsC_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_nphysicsC_init(Grid, Domain, Time, Time_steps, Thickness, T_prog, &
           ver_coordinate_class, agm_closure_lower_dept, agm_closure_upper_dept,   &
           smx, swidt, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_thickness_type),   intent(in)           :: Thickness
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  integer,                      intent(in)           :: ver_coordinate_class
  real,                         intent(in)           :: agm_closure_lower_dept
  real,                         intent(in)           :: agm_closure_upper_dept
  real,                         intent(in)           :: smx 
  real,                         intent(in)           :: swidt
  logical,                      intent(in), optional :: debug

  integer :: ioun, io_status, ierr  
  integer :: i, j, n
  integer :: num_methods=0
  
  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysicsC_mod (ocean_nphysicsC_init):already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  num_prog_tracers = size(T_prog(:))
  dtime            = Time_steps%dtime_t
  Dom => Domain
  Grd => Grid

  vert_coordinate_class   = ver_coordinate_class
  agm_closure_lower_depth = agm_closure_lower_dept
  agm_closure_upper_depth = agm_closure_upper_dept
  smax                    = smx 
  swidth                  = swidt
  pi_r                    = 1.0/pi
  grav_rho0_r             = 1.0/(grav*rho0)
  grav_rho0r              = grav*rho0r
  sqrt_grav               = sqrt(grav)
 
  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_nphysicsC_nml,IOSTAT=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit,ocean_nphysicsC_nml)  
  write (stdlogunit,ocean_nphysicsC_nml)
  ierr = check_nml_error(io_status,'ocean_nphysicsC_nml')
  call close_file (ioun)

#ifndef MOM4_STATIC_ARRAYS
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  if(use_this_module) then 
    call mpp_error(NOTE, &
    '==> from ocean_nphysicsC_mod: USING ocean_nphysicsC.')
    write(stdoutunit,'(1x,a)')    &
    '==> Note from ocean_nphysicsC_mod: USING ocean_nphysicsC.'
  else 
    call mpp_error(NOTE, &
    '==> from ocean_nphysicsC_mod: NOT using ocean_nphysicsC_mod.')
    write(stdoutunit,'(1x,a)')    &
    '==> Note from ocean_nphysicsC_mod: NOT using ocean_nphysicsC.'
    return
  endif 

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_nphysicsC_mod with debug_this_module=.true.'  
  endif 

  write(stdoutunit,'(/1x,a,f10.2)') &
  '==> Note from ocean_nphysicsC_mod: using forward time step of (secs)', dtime 

  if(transport_units=='Sv') then
      transport_convert=1.0e-9 
      transport_dims   = 'Sv (10^9 kg/s)'
  else
      transport_convert=1.0
      transport_dims   = 'kg/s'
  endif


  ! Useful constants 
  two_dtime_inv = 0.5/dtime            !for explicit piece of K33 
  swidthr       = 1.0/(swidth + epsln) !for slope taper function when dm_taper used 
  smax_swidthr  = smax*swidthr         !for slope taper function when dm_taper used 

  if(do_neutral_diffusion) then 
      write(stdoutunit,'(1x,a)') &
      ' ==>Note from ocean_nphysicsC_mod: computing neutral diffusion acting on each tracer.' 
  else
      write(stdoutunit,'(1x,a)') &
      ' ==>Note from ocean_nphysicsC_mod: NOT computing neutral diffusion. Are you sure this is what you wish?' 
  endif 
  if(do_gm_skewsion) then 
      if(gm_skewsion_modes) then 
          write(stdoutunit,'(1x,a)') &
          ' ==>Note from ocean_nphysicsC_mod: computing GM skewsion as a projection onto baroclinic modes.' 
          write(stdoutunit,'(1x,a,i4,a)') &
          '    Using ',number_bc_modes,' baroclinic modes to compute the eddy induced streamfunction.'
      elseif(gm_skewsion_bvproblem) then 
          write(stdoutunit,'(1x,a)') &
          ' ==>Note from ocean_nphysicsC_mod: computing GM skewsion w/ streamfunction computed by boundary value problem.' 
          if(bvp_bc_mode > number_bc_modes) then 
              number_bc_modes=bvp_bc_mode
          endif 
          if(bvp_bc_mode == 0) then 
             write(stdoutunit,'(1x,a)') &
             ' ==>Note from ocean_nphysicsC_mod: setting baroclinic mode phase speed to zero for BVP streamfunction.'
          endif 
      endif
  else
      write(stdoutunit,'(1x,a)') &
           ' ==>Note from ocean_nphysicsC_mod: NOT computing GM skewsion. Are you sure this is what you wish?' 
  endif

  if(neutral_physics_limit) then
      write(stdoutunit,'(1x,a)') &
      ' ==>Note from ocean_nphysicsC_mod: neutral_physics_limit=.true.' 
      write(stdoutunit,'(7x,a)') &
      'Will revert to horizontal diffusion for points where tracer is outside specified range.' 
  endif

  if(dm_taper .and. .not. gkw_taper) then
      write(stdoutunit,'(1x,a)') &
      ' ==>Note from ocean_nphysicsC_mod: dm_taper=.true. Will use the tanh scheme ' 
      write(stdoutunit,'(7x,a)') &
      'of Danabasoglu and McWilliams to taper neutral diffusion in steep sloped regions'
      dm_taper_const =1.0
      gkw_taper_const=0.0
  endif
  if(gkw_taper .and. .not. dm_taper) then
      write(stdoutunit,'(1x,a)') &
      ' ==>Note from ocean_nphysicsC_mod: gkw_taper=.true. Will use the quadratic scheme' 
      write(stdoutunit,'(7x,a)') &
      'of Gerdes, Koberle, and Willebrand to taper neutral diffusion in steep sloped regions'
      dm_taper_const =0.0
      gkw_taper_const=1.0
  endif
  if(gkw_taper .and. dm_taper) then
      dm_taper_const =0.0
      gkw_taper_const=0.0
      call mpp_error(FATAL, &
      '==>Error from ocean_nphysicsC_mod: gkw_taper and dm_taper cannot both be set true--choose only one.')
  endif

  if(neutral_linear_gm_taper) then
      write(stdoutunit,'(1x,a)')' ==>Running with a nontrivial GM transport in steep neutral slope regions.'
  endif

  if(diffusion_all_explicit) then
      write(stdoutunit,'(1x,a)') &
      ' ==> Warning: Running w/ diffusion_all_explicit=.true., which means compute K33 contribution'
      write(stdoutunit,'(7x,a)') &
      'to neutral diffusion explicitly in time.  This method is stable ONLY if taking'
      write(stdoutunit,'(7x,a)') &
      'very small time steps and/or running with just a single tracer.'
  endif

  if(number_bc_modes >= nk) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_nphysicsC_mod: reduce number_bc_modes to less than number of vertical levels.')
  endif 


  do n=1,num_prog_tracers
    T_prog(n)%neutral_physics_limit = neutral_physics_limit
  enddo 

  allocate( dTdx(num_prog_tracers) )
  allocate( dTdy(num_prog_tracers) )
  allocate( dTdz(num_prog_tracers) )
  allocate( fz1(num_prog_tracers) )
  allocate( fz2(num_prog_tracers) )
  allocate( flux_x(num_prog_tracers) )
  allocate( flux_y(num_prog_tracers) )

  call set_ocean_domain(Dom_flux,Grid,xhalo=Dom%xhalo,yhalo=Dom%yhalo,name='flux dom neutral',maskmap=Dom%maskmap)

#ifndef MOM4_STATIC_ARRAYS
  allocate (dtew(isd:ied,jsd:jed,0:1))
  allocate (dtns(isd:ied,jsd:jed,0:1))
  allocate (dtwedyt(isd:ied,jsd:jed,0:1))
  allocate (dxtdtsn(isd:ied,jsd:jed,0:1))
  allocate (grid_length(isd:ied,jsd:jed))
  allocate (delqc(isd:ied,jsd:jed,nk,0:1))
  allocate (dzwtr(isd:ied,jsd:jed,0:nk))
  allocate (aredi_array(isd:ied,jsd:jed,nk)) 
  allocate (agm_array(isd:ied,jsd:jed,nk))   
  allocate (ah_array(isd:ied,jsd:jed))
  allocate (bczone_radius(isd:ied,jsd:jed))
  allocate (rossby_radius(isd:ied,jsd:jed))
  allocate (rossby_radius_raw(isd:ied,jsd:jed))
  allocate (eady_termx(isd:ied,jsd:jed,nk))
  allocate (eady_termy(isd:ied,jsd:jed,nk))
  allocate (baroclinic_termx(isd:ied,jsd:jed))
  allocate (baroclinic_termy(isd:ied,jsd:jed))
  allocate (tx_trans_gm(isd:ied,jsd:jed,nk))
  allocate (ty_trans_gm(isd:ied,jsd:jed,nk))

  allocate (psix(isd:ied,jsd:jed,nk,0:1))
  allocate (psiy(isd:ied,jsd:jed,nk,0:1))

  allocate (bv_freq(isd:ied,jsd:jed,nk))

  allocate (drhodT(isd:ied,jsd:jed,nk))
  allocate (drhodS(isd:ied,jsd:jed,nk))
  allocate (drhodzb(isd:ied,jsd:jed,nk,0:1))
  allocate (drhodzh(isd:ied,jsd:jed,nk,0:1))
  allocate (tensor_31(isd:ied,jsd:jed,nk,0:1,0:1))
  allocate (tensor_32(isd:ied,jsd:jed,nk,0:1,0:1))
  allocate (K33_implicit(isd:ied,jsd:jed,nk)) 
  allocate (K33_explicit(isd:ied,jsd:jed,nk)) 

  do n=1,num_prog_tracers
    allocate ( dTdx(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( dTdy(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( dTdz(n)%field(isd:ied,jsd:jed,0:nk) )
    allocate ( fz1(n)%field(isd:ied,jsd:jed) )
    allocate ( fz2(n)%field(isd:ied,jsd:jed) )
    allocate ( flux_x(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( flux_y(n)%field(isd:ied,jsd:jed,nk) )
  enddo 

  allocate(eddy_depth(isd:ied,jsd:jed))
  allocate(ksurf_blayer(isd:ied,jsd:jed))

#endif

  do n=1,num_prog_tracers 
    dTdx(n)%field(:,:,:)   = 0.0
    dTdy(n)%field(:,:,:)   = 0.0
    dTdz(n)%field(:,:,:)   = 0.0
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0
  enddo  

  psix               = 0.0
  psiy               = 0.0
  bv_freq            = 0.0
  tx_trans_gm        = 0.0
  ty_trans_gm        = 0.0
  eddy_depth         = 0.0
  ksurf_blayer       = 1

  dtew(:,:,0) = Grd%dtw(:,:)
  dtew(:,:,1) = Grd%dte(:,:)
  dtns(:,:,0) = Grd%dts(:,:)
  dtns(:,:,1) = Grd%dtn(:,:)

  dtwedyt(:,:,:) = 0.0
  dtwedyt(:,:,0) = Grd%dte(:,:)*Grd%dyt(:,:)
  do i=isc-1,iec
    dtwedyt(i,:,1) = Grd%dtw(i+1,:)*Grd%dyt(i+1,:)
  enddo

  dxtdtsn(:,:,:) = 0.0
  dxtdtsn(:,:,0) = Grd%dxt(:,:)*Grd%dtn(:,:)
  do j=jsc-1,jec
    dxtdtsn(:,j,1) = Grd%dxt(:,j+1)*Grd%dts(:,j+1)
  enddo

  grid_length(:,:) = 2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))

  eady_termx(:,:,:)            = 0.0
  eady_termy(:,:,:)            = 0.0
  baroclinic_termx(:,:)        = 0.0
  baroclinic_termy(:,:)        = 0.0

  rossby_radius(:,:)           = 0.0
  rossby_radius_raw(:,:)       = 0.0
  bczone_radius(:,:)           = 0.0
  agm_array(:,:,:)             = 0.0
  aredi_array(:,:,:)           = 0.0
  ah_array(:,:)                = 0.0

  call ocean_nphysics_coeff_init(Time, Thickness, rossby_radius, rossby_radius_raw, &
                         bczone_radius, agm_array, aredi_array, ah_array)

  drhodT(:,:,:)        = 0.0
  drhodS(:,:,:)        = 0.0
  drhodzb(:,:,:,:)     = 0.0
  drhodzh(:,:,:,:)     = 0.0
  tensor_31(:,:,:,:,:) = 0.0
  tensor_32(:,:,:,:,:) = 0.0
  K33_implicit(:,:,:)  = 0.0           
  K33_explicit(:,:,:)  = 0.0           

  allocate (bc_mode(isd:ied,jsd:jed,nk,number_bc_modes))
  allocate (bc_speed(isd:ied,jsd:jed,number_bc_modes))
  allocate (bc_speed2(isd:ied,jsd:jed))
  bc_mode   = 0.0
  bc_speed  = 0.0
  bc_speed2 = 0.0


  index_temp=-1;index_salt=-1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp == -1 .or. index_salt == -1) then 
     call mpp_error(FATAL, &
     '==>Error: temp and/or salt not identified in call to ocean_nphysicsC_init')
  endif 

  ! initialize clock ids 
  id_clock_ndiffuse         = mpp_clock_id('(Ocean neutral: diffuse)'   ,grain=CLOCK_ROUTINE)
  id_clock_gm_modes         = mpp_clock_id('(Ocean neutral: gm_modes)'  ,grain=CLOCK_ROUTINE)
  id_clock_neutral_blayer   = mpp_clock_id('(Ocean neutral: blayer)'    ,grain=CLOCK_ROUTINE)
  id_clock_fz_terms         = mpp_clock_id('(Ocean neutral: fz-terms)'  ,grain=CLOCK_ROUTINE)
  id_clock_fx_flux_ndiffuse = mpp_clock_id('(Ocean ndiffuse: fx-flux)'  ,grain=CLOCK_ROUTINE)
  id_clock_fy_flux_ndiffuse = mpp_clock_id('(Ocean ndiffuse: fy-flux)'  ,grain=CLOCK_ROUTINE)
  id_clock_fz_flux_ndiffuse = mpp_clock_id('(Ocean ndiffuse: fz-flux)'  ,grain=CLOCK_ROUTINE)
  id_clock_fx_flux_gm       = mpp_clock_id('(Ocean gm_modes: fx-flux)'  ,grain=CLOCK_ROUTINE)
  id_clock_fy_flux_gm       = mpp_clock_id('(Ocean gm_modes: fy-flux)'  ,grain=CLOCK_ROUTINE)
  id_clock_fz_flux_gm       = mpp_clock_id('(Ocean gm_modes: fz-flux)'  ,grain=CLOCK_ROUTINE)
  id_clock_bc_modes         = mpp_clock_id('(Ocean gm_modes: bc_modes)' ,grain=CLOCK_ROUTINE)
  id_clock_psi_modes        = mpp_clock_id('(Ocean gm_modes: psi_modes)',grain=CLOCK_ROUTINE)
  id_clock_psi_bvp          = mpp_clock_id('(Ocean gm_modes: psi_bvp)'  ,grain=CLOCK_ROUTINE)

  ! register fields for diagnostic output 

  id_rescale_psi_x = -1
  id_rescale_psi_x = register_diag_field ('ocean_model', 'rescale_psi_x',       &
                    Grd%tracer_axes(1:2), Time%model_time,                      &
                    'Dimensionless rescaling of mode-1 psix for GM', 'unitless',&
                    missing_value=missing_value, range=(/-1.0,1.e10/))

  id_rescale_psi_y = -1
  id_rescale_psi_y = register_diag_field ('ocean_model', 'rescale_psi_y',       &
                    Grd%tracer_axes(1:2), Time%model_time,                      &
                    'Dimensionless rescaling of mode-1 psiy for GM', 'unitless',&
                    missing_value=missing_value, range=(/-1.0,1.e10/))

  id_k33_explicit = -1
  id_k33_explicit = register_diag_field ('ocean_model', 'k33_explicit', &
                    Grd%tracer_axes_wt(1:3), Time%model_time,           &
                    'K33_explicit tensor element', 'm^2/sec',           &
                    missing_value=missing_value, range=(/-10.0,1.e20/))

  id_N_squared = -1
  id_N_squared = register_diag_field ('ocean_model', 'N_squared', &
            Grd%tracer_axes_wt(1:3), Time%model_time,             &
            'squared buoyancy frequency', '1/s^2',                &
            missing_value=missing_value, range=(/-10.0,1.e10/))

  id_bv_freq = -1
  id_bv_freq = register_diag_field ('ocean_model', 'bv_freq',          &
            Grd%tracer_axes_wt(1:3), Time%model_time,                  &
            'buoyancy frequency for computing baroclinic modes', '1/s',&
            missing_value=missing_value, range=(/-10.0,1.e10/))

  id_psix_gm_modes = register_diag_field ('ocean_model', 'psix_gm_modes', &
         Grd%tracer_axes_flux_y(1:3), Time%model_time,                    &
         'i-comp of modal decomposed GM streamfunction',                  &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))

  id_psiy_gm_modes = register_diag_field ('ocean_model', 'psiy_gm_modes', &
         Grd%tracer_axes_flux_x(1:3), Time%model_time,                    &
         'j-comp of modal decomposed GM streamfunction',                  &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))

  id_psix_gm_bvp = register_diag_field ('ocean_model', 'psix_gm_bvp', &
         Grd%tracer_axes_flux_y(1:3), Time%model_time,                &
         'i-comp of GM streamfunction computed via a BVP',            &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))

  id_psiy_gm_bvp = register_diag_field ('ocean_model', 'psiy_gm_bvp', &
         Grd%tracer_axes_flux_x(1:3), Time%model_time,                &
         'j-comp of GM streamfunction computed via a BVP',            &
         'm^2/sec', missing_value=missing_value, range=(/-1.e10,1e10/))

  id_tx_trans_gm = -1
  id_tx_trans_gm = register_diag_field ('ocean_model', 'tx_trans_gm',      &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time,           &
                   'T-cell mass i-transport from GM',trim(transport_dims), &
                   missing_value=missing_value, range=(/-1e8,1.e8/))

  id_ty_trans_gm = -1
  id_ty_trans_gm = register_diag_field ('ocean_model', 'ty_trans_gm',     &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time,          &
                   'T-cell mass j-transport from GM',trim(transport_dims),&
                   missing_value=missing_value, range=(/-1e8,1.e8/))

  id_eddy_depth = -1
  id_eddy_depth = register_diag_field ('ocean_model','eddy_depth', &
                  Grd%tracer_axes(1:2), Time%model_time,           &
                  'eddy penetration depth', 'm',                   &
                  missing_value=missing_value, range=(/-1.0,1.e8/))

  id_neutral_rho_ndiffuse = -1
  id_neutral_rho_ndiffuse = register_diag_field ('ocean_model', 'neutral_rho_ndiffuse',     &
                          Grd%tracer_axes(1:3), Time%model_time,                            &
                          'update of neutral density from explicit in time neutral diffuse',&
                          'rho*rho_dz/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_wdian_rho_ndiffuse = -1
  id_wdian_rho_ndiffuse = register_diag_field ('ocean_model', 'wdian_rho_ndiffuse',                &
                          Grd%tracer_axes(1:3), Time%model_time,                                   &
                          'dianeutral velocity component due to explicit in time  neutral diffuse',&
                          'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_neutral_rho_gm_modes = -1
  id_neutral_rho_gm_modes = register_diag_field ('ocean_model', 'neutral_rho_gm_modes',     &
                          Grd%tracer_axes(1:3), Time%model_time,                            &
                          'update of neutral density from GM skewsion projected into modes',&
                          'rho*rho_dz/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))
  id_wdian_rho_gm_modes = -1
  id_wdian_rho_gm_modes = register_diag_field ('ocean_model', 'wdian_rho_gm_modes',               &
                          Grd%tracer_axes(1:3), Time%model_time,                                  &
                          'dianeutral velocity component due to GM skewsion projected into modes',&
                          'm/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_gradx_rho = -1
  id_gradx_rho = register_diag_field ('ocean_model', 'gradx_rho', &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,       &
              'd(rho)/dx for neutral density', 'kg/m^4',          &
              missing_value=missing_value, range=(/-1e6,1e6/))

  id_grady_rho = -1
  id_grady_rho = register_diag_field ('ocean_model', 'grady_rho', &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,       &
              'd(rho)/dy for neutral density', 'kg/m^4',          &
              missing_value=missing_value, range=(/-1e6,1e6/))

  id_gm_eddy_ke_source=-1
  id_gm_eddy_ke_source = register_diag_field ('ocean_model', 'gm_eddy_ke_source',                        &
               Grd%tracer_axes(1:2), Time%model_time,                                                    &
               'vertical sum of rho*dz*[(N*psi)^2 + (c*psi_z)^2]/agm = eddy ke source from BVP', 'W/m^2',&
               missing_value=missing_value, range=(/-1.e1,1.e15/),                                       &
               standard_name='tendency_of_ocean_eddy_kinetic_energy_content_due_to_bolus_transport')


  allocate (id_k33_implicit(num_prog_tracers))
  allocate (id_neutral_diffuse(num_prog_tracers))
  allocate (id_neutral_gm_modes(num_prog_tracers))
  id_k33_implicit    = -1
  id_neutral_diffuse = -1
  id_neutral_gm_modes= -1
  do n=1,num_prog_tracers
     id_k33_implicit(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_k33_implicit', &
                          Grd%tracer_axes_wt(1:3), Time%model_time,                                  &
                          'K33_implicit tensor element for '//trim(T_prog(n)%name),                  &
                          'm^2/sec', missing_value=missing_value, range=(/-10.0,1.e10/))

     if (T_prog(n)%name == 'temp') then
       id_neutral_diffuse(n) = register_diag_field ('ocean_model', 'neutral_diffuse_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                        &
                               'rho*dzt*cp*explicit neutral diffusion tendency (heating)',                   &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
       id_neutral_gm_modes(n) = register_diag_field ('ocean_model', 'neutral_gm_modes_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                          &
                               'rho*dzt*cp*GM as projected into modes tendency (heating)',                     &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
     else 
       id_neutral_diffuse(n) = register_diag_field ('ocean_model', 'neutral_diffuse_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                        &
                               'rho*dzt*explicit neutral diffuseion tendency for '//trim(T_prog(n)%name),    &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
       id_neutral_gm_modes(n) = register_diag_field ('ocean_model', 'neutral_gm_modes_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time,                                          &
                               'rho*dzt*GM as projected into modes tendency for '//trim(T_prog(n)%name),       &
                               trim(T_prog(n)%flux_units), missing_value=missing_value, range=(/-1.e10,1.e10/))
     endif 
  enddo 

  allocate (id_flux_x_ndiffuse(num_prog_tracers))
  allocate (id_flux_y_ndiffuse(num_prog_tracers))
  allocate (id_flux_x_gm_modes(num_prog_tracers))
  allocate (id_flux_y_gm_modes(num_prog_tracers))
  allocate (id_flux_x_ndiffuse_int_z(num_prog_tracers))
  allocate (id_flux_y_ndiffuse_int_z(num_prog_tracers))
  allocate (id_flux_x_gm_modes_int_z(num_prog_tracers))
  allocate (id_flux_y_gm_modes_int_z(num_prog_tracers))
  id_flux_x_ndiffuse = -1
  id_flux_y_ndiffuse = -1
  id_flux_x_gm_modes = -1
  id_flux_y_gm_modes = -1
  id_flux_x_ndiffuse_int_z = -1
  id_flux_y_ndiffuse_int_z = -1
  id_flux_x_gm_modes_int_z = -1
  id_flux_y_gm_modes_int_z = -1

  do n=1,num_prog_tracers
     if(n == index_temp) then 
         id_flux_x_ndiffuse(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_ndiffuse',               &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,          &
              'cp*ndiffuse_xflux*dyt*rho_dzt*temp',                  &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y_ndiffuse(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_ndiffuse',               &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,          &
              'cp*ndiffuse_yflux*dxt*rho_dzt*temp',                  &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_x_gm_modes(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_gm_modes',               &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,          &
              'cp*gm_modes_xflux*dyt*rho_dzt*temp',                  &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y_gm_modes(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_gm_modes',               & 
              Grd%tracer_axes_flux_y(1:3), Time%model_time,          &
              'cp*gm_modes_yflux*dxt*rho_dzt*temp',                  &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))

         id_flux_x_ndiffuse_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_ndiffuse_int_z',               &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                &
              'z-integral cp*ndiffuse_xflux*dyt*rho_dzt*temp',             &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y_ndiffuse_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_ndiffuse_int_z',               &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                &
              'z-integral cp*ndiffuse_yflux*dxt*rho_dzt*temp',             &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_x_gm_modes_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_gm_modes_int_z',               &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                &
              'z-integral cp*gm_modes_xflux*dyt*rho_dzt*temp',             &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
         id_flux_y_gm_modes_int_z(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_gm_modes_int_z',               &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                &
              'z-integral cp*gm_modes_yflux*dyt*rho_dzt*temp',             &
              'Watt', missing_value=missing_value, range=(/-1.e18,1.e18/))
     else
         id_flux_x_ndiffuse(n) = register_diag_field ('ocean_model',         &
              trim(T_prog(n)%name)//'_xflux_ndiffuse',                       &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,                  &
              'ndiffuse_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value, &
              range=(/-1.e18,1.e18/))
         id_flux_y_ndiffuse(n) = register_diag_field ('ocean_model',         &
              trim(T_prog(n)%name)//'_yflux_ndiffuse',                       &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,                  &
              'ndiffuse_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value, &
              range=(/-1.e18,1.e18/))
         id_flux_x_gm_modes(n) = register_diag_field ('ocean_model',         &
              trim(T_prog(n)%name)//'_xflux_gm_modes',                       &
              Grd%tracer_axes_flux_x(1:3), Time%model_time,                  &
              'gm_modes_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value, &
              range=(/-1.e18,1.e18/))
         id_flux_y_gm_modes(n) = register_diag_field ('ocean_model',         &
              trim(T_prog(n)%name)//'_yflux_gm_modes',                       &
              Grd%tracer_axes_flux_y(1:3), Time%model_time,                  &
              'gm_modes_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value, &
              range=(/-1.e18,1.e18/))
 
         id_flux_x_ndiffuse_int_z(n) = register_diag_field ('ocean_model',              &
              trim(T_prog(n)%name)//'_xflux_ndiffuse_int_z',                            &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                             &
              'z-integral ndiffuse_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,            &
              range=(/-1.e18,1.e18/))
         id_flux_y_ndiffuse_int_z(n) = register_diag_field ('ocean_model',             &
              trim(T_prog(n)%name)//'_yflux_ndiffuse_int_z',                           &
              Grd%tracer_axes_flux_y(1:2), Time%model_time,                            &
              'z-integral ndiffuse_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name),&
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,           &
              range=(/-1.e18,1.e18/))
         id_flux_x_gm_modes_int_z(n) = register_diag_field ('ocean_model',              &
              trim(T_prog(n)%name)//'_xflux_gm_modes_int_z',                            &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                             &
              'z-integral gm_modes_xflux*dyt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,            & 
              range=(/-1.e18,1.e18/))
         id_flux_y_gm_modes_int_z(n) = register_diag_field ('ocean_model',              &
              trim(T_prog(n)%name)//'_yflux_gm_modes_int_z',                            &
              Grd%tracer_axes_flux_x(1:2), Time%model_time,                             &
              'z-integral gm_modes_yflux*dxt*rho_dzt*tracer for'//trim(T_prog(n)%name), &
              trim(T_prog(n)%units)//' kg/sec', missing_value=missing_value,            &
              range=(/-1.e18,1.e18/))

     endif

  enddo

  allocate (id_bc_mode(number_bc_modes))
  allocate (id_bc_speed(number_bc_modes))
  id_bc_mode  = -1
  id_bc_speed = -1
  do n=1,number_bc_modes
         id_bc_mode(n) = register_diag_field ('ocean_model',     &
              'baroclinic_mode_'//string(n),                     &
              Grd%tracer_axes(1:3), Time%model_time,             &
              'Dimensionless baroclinic mode number '//string(n),&
              'Dimensionless', missing_value=missing_value, range=(/-1.e2,1.e2/))
         id_bc_speed(n) = register_diag_field ('ocean_model',  &
              'baroclinic_speed_'//string(n),                  &
              Grd%tracer_axes(1:2), Time%model_time,           &
              'Speed for baroclinic mode number '//string(n),  &
              'm/s', missing_value=missing_value, range=(/-1.e2,1.e2/))
  enddo


end subroutine ocean_nphysicsC_init
! </SUBROUTINE>  NAME="ocean_nphysicsC_init"


!#######################################################################
! <FUNCTION NAME="nphysicsC">
!
! <DESCRIPTION>
! This function computes the thickness weighted and density weighted
! time tendency for tracer from neutral physics.  Full discussion
! and details are provided by Griffies (2008). 
!
! Here is a brief summary.  
!
!---How the neutral diffusive flux components are computed:
!
! The vertical flux component is split into diagonal (3,3) and 
! off-diagonal (3,1) and (3,2) terms. The off-diagonal (3,1) and (3,2) 
! terms are included explicitly in time. The main contribution from the 
! (3,3) term to the time tendency is included implicitly in time 
! along with the usual contribution from diapycnal processes 
! (vertical mixing schemes).  This is the K33_implicit term.
! This approach is necessary with high vertical resolution, as 
! noted by Cox (1987).  However, splitting the vertical flux into 
! an implicit and explicit piece compromises the 
! integrity of the vertical flux component (see Griffies et al. 1998).
! So to minimize the disparity engendered by this split, the portion of 
! K33 that can be stably included explicitly in time is computed along 
! with the (3,1) and (3,2) terms. 
! 
! All other terms in the mixing tensor are included explicitly in time
! using a forward time step as required for temporal stability of 
! numerical diffusive processes.  
!
! The off-diagonal terms in the horizontal flux components, and all terms
! in the vertical flux component, are tapered in regions of steep neutral
! slope according to the requirements of linear stability.  MOM4 allows for 
! choice of two tapering schemes:
! (a) the tanh taper of Danabasoglu and McWilliams (1995)
! (b) the quadratic scheme of Gerdes, Koberle, and Willebrand (1991)
! Linear stability is far less stringent on the diagonal (1,1) and (2,2)
! part of the horizontal flux.  Indeed, these terms in practice need
! not be tapered in steep sloped regions. 
!
!---How the skew diffusive flux components are computed:
!
! The skew flux components are purely off-diagonal.  
! They are computed based on a vector streamfunction which 
! is built from a sum of baroclinic modes. 
! It is this part of the calculation that differs from 
! ocean_nphysicsA and ocean_nphysicsB.
! </DESCRIPTION>
!
subroutine nphysicsC (Time, Thickness, Dens, rho, T_prog, &
                      surf_blthick, gm_diffusivity, baroclinic_rossby_radius)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  real, dimension(isd:,jsd:,:), intent(in)    :: rho
  real, dimension(isd:,jsd:),   intent(in)    :: surf_blthick
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  real, dimension(isd:,jsd:,:), intent(inout) :: gm_diffusivity
  real, dimension(isd:,jsd:),   intent(inout) :: baroclinic_rossby_radius

  integer :: i, j, k
  integer :: tau, taum1

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysicsC (neutral_physics): needs initialization')
  endif 

  if (size(T_prog(:)) /= num_prog_tracers) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysicsC (neutral_physics): inconsistent size of T_prog')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  ! time dependent delqc geometric factor 
  do k=1,nk
    delqc(:,:,k,0) = Grd%fracdz(k,0)*Thickness%rho_dzt(:,:,k,tau)
    delqc(:,:,k,1) = Grd%fracdz(k,1)*Thickness%rho_dzt(:,:,k,tau)
  enddo

  ! time dependent inverse dzwt 
  do k=0,nk
    dzwtr(:,:,k) = 1.0/Thickness%dzwt(:,:,k) 
  enddo

  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           drhodT(i,j,k) = Dens%drhodT(i,j,k)
           drhodS(i,j,k) = Dens%drhodS(i,j,k)
        enddo
     enddo
  enddo

  call tracer_derivs(Time, taum1, dtime, drhodT, drhodS, T_prog, dzwtr, &
                     dTdx, dTdy, dTdz, drhodzb, drhodzh) 

  call neutral_slopes(Time, dTdx, dTdy, drhodT, drhodS, drhodzb, tensor_31, tensor_32)

  call neutral_blayer(Time, Thickness, surf_blthick)
 
  call fz_terms(Time, Thickness, T_prog, rho)

  if(do_neutral_diffusion) then 
     call compute_ndiffusion(Time, Thickness, Dens, T_prog)
  endif 

  ! gm_diffusivity passed to neutral_physics for use in computing form drag
  gm_diffusivity(:,:,:) = 0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
          gm_diffusivity(i,j,k) = agm_array(i,j,k)
        enddo
     enddo
  enddo

  if(do_gm_skewsion) then 
      call baroclinic_modes(Time, Thickness, Dens) 
      if(gm_skewsion_modes) then 
          call compute_psi_modes(Time, Dens, Thickness, T_prog)
      else 
          call compute_psi_bvp(Time, Dens, Thickness, T_prog)
      endif
      call compute_gmskewsion(Time, Thickness, Dens, T_prog)
  endif
  
  ! compute Eady growth rate and baroclinicity for use in next time step 
  call compute_eady_rate(Time, Thickness, T_prog, Dens, eady_termx, eady_termy)
  call compute_baroclinicity(Time%model_time, baroclinic_termx, baroclinic_termy)

  ! compute rossby radius for use in next time step 
  call compute_rossby_radius(Thickness, dTdz, Time%model_time, drhodT, drhodS, &
                             rossby_radius, rossby_radius_raw)
  baroclinic_rossby_radius(:,:) = rossby_radius_raw(:,:)

  ! compute width of baroclinic zone for use in next time step 
  call compute_bczone_radius(Time%model_time, bczone_radius)

  ! update closure-based diffusivity for next time step 
  call compute_diffusivity(Time%model_time, ksurf_blayer, Dens%drhodz_zt, rossby_radius, &
                           bczone_radius, agm_array, aredi_array)


end subroutine nphysicsC
! </FUNCTION> NAME="nphysicsC"


!#######################################################################
! <SUBROUTINE NAME="neutral_blayer">
!
! <DESCRIPTION>
! Subroutine computes the boundary layer as determined by 
! 1. depth within which typical mesoscale eddies are partially outcropped
! 2. depth within which vertical mixing scheme (e.g., kpp) computes a boundary layer
!
! Determine depth over which mesoscale eddies feel the ocean 
! surface.  This depth is a function of the neutral slope 
! and the Rossby radius.  This depth is called "eddy_depth".
! The algorithm for computing this depth is taken from 
! the appendix to Large etal, 1997 JPO vol 27, 2418-2447. 
! 
! In addition to considering mesoscale eddy lengths,
! include the possibility that the diabatic vertical
! mixing (e.g., KPP) produces a mixed layer depth that is 
! deeper than the depth that mesoscale eddies feel the ocean 
! surface.  Include this surf_blthick in the considerations so 
! to determine the depth of this generalized "boundary layer" 
! and the neutral slope at the base of the boundary layer. 
!
! Note: Only consider surface boundary layers here.  
!
! This subroutine is a modification of that in ocean_nphysicsA. 
! Here, we only compute the eddy_depth based on the
! algorithm in Large etal.  We do not compute an eddy
! depth which is also a function of smax. that is, we 
! remove the ocean_nphysicsA portion of the calculation 
! that sits inside the neutral_linear_gm_taper if-test.   
! 
! </DESCRIPTION>
!
subroutine neutral_blayer(Time, Thickness, surf_blthick)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  real, dimension(isd:,jsd:),   intent(in) :: surf_blthick

  logical :: slopeok(isc:iec)
  logical :: within_interior 
  integer :: i, j, k, ip, jq, kb, kr, kk, kkpkr
  real :: depth, slope, absslope 
  real :: thick_31(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_32(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_13(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_23(isd:ied,jsd:jed,0:1,0:1)

  eddy_depth  = 0.0
  if(.not. neutral_eddy_depth) return

  call mpp_clock_begin(id_clock_neutral_blayer)

  thick_31    = Grd%zw(1)
  thick_32    = Grd%zw(1)
  thick_13    = Grd%zw(1)
  thick_23    = Grd%zw(1)

  ! 31-triads 
  do ip=0,1
     do kr=0,1
        do j=jsc,jec
           slopeok(:)=.false.
           kkloop31a:     do kk=1,nk-1
              do i=isc,iec 

                 absslope = abs(tensor_31(i,j,kk,ip,kr))
                 depth    = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                 if(depth > Thickness%depth_zwt(i,j,kk) .and. .not. slopeok(i)) then
                     thick_31(i,j,ip,kr) = Thickness%depth_zwt(i,j,kk)
                 else 
                     slopeok(i)=.true.
                 endif
              enddo
              if(all(slopeok(:))) exit kkloop31a

           enddo kkloop31a
        enddo
     enddo
  enddo

  ! 13-triads 
  do kr=0,1
     do ip=0,1 
        do j=jsc,jec
           slopeok(:)=.false.

           kkloop13a:     do kk=1,nk
              kkpkr = min(kk+kr,nk)
              do i=isc,iec

                 slope = -Grd%tmask(i+ip,j,kkpkr)                            &
                      *(drhodT(i+ip,j,kk)*dTdx(index_temp)%field(i,j,kk)   + &
                        drhodS(i+ip,j,kk)*dTdx(index_salt)%field(i,j,kk))    &
                      /(drhodT(i+ip,j,kk)*dTdz(index_temp)%field(i+ip,j,kk-1+kr) +&
                        drhodS(i+ip,j,kk)*dTdz(index_salt)%field(i+ip,j,kk-1+kr) - epsln)

                 absslope = abs(slope)  

                 depth = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                 if(depth > Thickness%depth_zwt(i,j,kk) .and. .not.slopeok(i)) then 
                     thick_13(i,j,ip,kr) = Thickness%depth_zwt(i,j,kk)
                 else 
                     slopeok(i)=.true.
                 endif

              enddo
              if(all(slopeok(:))) exit kkloop13a
           enddo kkloop13a

        enddo
     enddo
  enddo

  ! 32-triads 
  do jq=0,1
     do kr=0,1
        do j=jsc,jec
           slopeok(:)=.false.
           kkloop32a:      do kk=1,nk-1
              do i=isc,iec 

                 absslope = abs(tensor_32(i,j,kk,jq,kr))
                 depth    = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                 if(depth > Thickness%depth_zwt(i,j,kk) .and. .not. slopeok(i) ) then
                     thick_32(i,j,jq,kr) = Thickness%depth_zwt(i,j,kk)
                 else 
                     slopeok(i) = .true.
                 endif
              enddo
              if(all(slopeok(:))) exit kkloop32a
           enddo kkloop32a
        enddo
     enddo
  enddo

  ! 23-triads 
  do kr=0,1
     do jq=0,1  
        do j=jsc,jec
           slopeok(:)=.false.

           kkloop23a:     do kk=1,nk
              do i=isc,iec
                 kkpkr = min(kk+kr,nk)

                 slope = -Grd%tmask(i,j+jq,kkpkr)&
                      *(drhodT(i,j+jq,kk)*dTdy(index_temp)%field(i,j,kk)+ &
                        drhodS(i,j+jq,kk)*dTdy(index_salt)%field(i,j,kk)) & 
                      /(drhodT(i,j+jq,kk)*dTdz(index_temp)%field(i,j+jq,kk-1+kr) +&
                        drhodS(i,j+jq,kk)*dTdz(index_salt)%field(i,j+jq,kk-1+kr) - epsln)  
                 absslope = abs(slope)  

                 depth = max(rossby_radius(i,j)*absslope, max(turb_blayer_min,surf_blthick(i,j)))
                 if(depth > Thickness%depth_zwt(i,j,kk) .and. .not.slopeok(i)) then 
                     thick_23(i,j,jq,kr) = Thickness%depth_zwt(i,j,kk) 
                 else
                     slopeok(i)=.true.
                 endif

              enddo
              if(all(slopeok(:))) exit kkloop23a
           enddo kkloop23a
        enddo
     enddo
  enddo

  ! max of triad depth defines eddy_depth.
  do j=jsc,jec
     do i=isc,iec
        eddy_depth(i,j) = max(Thickness%depth_zwt(i,j,1), &
                    thick_31(i,j,0,0), thick_31(i,j,0,1), &
                    thick_31(i,j,1,0), thick_31(i,j,1,1), &
                    thick_32(i,j,0,0), thick_32(i,j,0,1), &
                    thick_32(i,j,1,0), thick_32(i,j,1,1), &
                    thick_13(i,j,0,0), thick_13(i,j,0,1), &
                    thick_13(i,j,1,0), thick_13(i,j,1,1), &
                    thick_23(i,j,0,0), thick_23(i,j,0,1), &
                    thick_23(i,j,1,0), thick_23(i,j,1,1), &
                    turb_blayer_min)
     enddo
  enddo
  call mpp_update_domains(eddy_depth, Dom%domain2d) 


  ! save the k-level surface transition layer. 
  ! initialize ksurf_blayer to 1 rather than 0, to avoid 
  ! accessing the zeroth element of an array.   
  ksurf_blayer(:,:) = 1
  wrk1_2d(:,:)      = 0.0
  do j=jsd,jed
     do i=isd,ied
        kb=Grd%kmt(i,j)
        within_interior=.false.
        ksurf_blayer_loop: do k=1,kb
           if(Thickness%depth_zwt(i,j,k) > eddy_depth(i,j)) then 
               ksurf_blayer(i,j) = k
               wrk1_2d(i,j)      = float(k)
               exit ksurf_blayer_loop
           endif
        enddo ksurf_blayer_loop
     enddo
  enddo

  if (id_eddy_depth >  0) then 
       used = send_data (id_eddy_depth, eddy_depth(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1),   &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


  call mpp_clock_end(id_clock_neutral_blayer)

end subroutine neutral_blayer
! </SUBROUTINE> NAME="neutral_blayer"


!#######################################################################
! <SUBROUTINE NAME="compute_ndiffusion">
!
! <DESCRIPTION>
! Subroutine to compute the tendency from neutral diffusion.
! </DESCRIPTION>
!
subroutine compute_ndiffusion(Time, Thickness, Dens, T_prog)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  real, dimension(isd:ied,jsd:jed) :: tmp_flux
  real    :: temporary 
  integer :: i,j,k,n,tau

  call mpp_clock_begin(id_clock_ndiffuse)
 

  tau=Time%tau 

  do n=1,num_prog_tracers  
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0
    T_prog(n)%wrk1(:,:,:)  = 0.0
  enddo 
  
  do k=1,nk
     call fx_flux_ndiffuse(Time, Thickness, T_prog, k)
     call fy_flux_ndiffuse(Time, Thickness, T_prog, k)
  enddo
  if (Grd%tripolar) then 
      do n=1,num_prog_tracers
         call mpp_update_domains(flux_x(n)%field(:,:,:), flux_y(n)%field(:,:,:), Dom_flux%domain2d, &
              gridtype=CGRID_NE, complete=T_prog(n)%complete) 
      enddo
  endif

  do k=1,nk
     call mpp_clock_begin(id_clock_fz_flux_ndiffuse)
     call fz_flux_ndiffuse(T_prog,k)
     call mpp_clock_end(id_clock_fz_flux_ndiffuse)
     do n=1,num_prog_tracers
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%wrk1(i,j,k) =                                                 & 
                   Grd%tmask(i,j,k)                                                   &
                   *(fz1(n)%field(i,j)-fz2(n)%field(i,j)                              &
                   +( flux_x(n)%field(i,j,k)-flux_x(n)%field(i-1,j,k)                 &
                     +flux_y(n)%field(i,j,k)-flux_y(n)%field(i,j-1,k) )*Grd%datr(i,j) &
                    )
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
              fz1(n)%field(i,j) = fz2(n)%field(i,j)
           enddo
        enddo

     enddo
  enddo


  ! diagnostics 
  do n=1,num_prog_tracers

     if(id_neutral_diffuse(n) > 0) then 
        used = send_data (id_neutral_diffuse(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion, &
               Time%model_time, rmask=Grd%tmask(:,:,:),                                      &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

     ! send fluxes to diag_manager 
     ! minus sign is due to a MOM-convention for physics fluxes  
     if(id_flux_x_ndiffuse(n) > 0) then 
        used = send_data (id_flux_x_ndiffuse(n), -1.0*T_prog(n)%conversion*flux_x(n)%field(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:),                                            &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_flux_y_ndiffuse(n) > 0) then 
        used = send_data (id_flux_y_ndiffuse(n), -1.0*T_prog(n)%conversion*flux_y(n)%field(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:),                                            &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_flux_x_ndiffuse_int_z(n) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) + flux_x(n)%field(isc:iec,jsc:jec,k)
        enddo        
        used = send_data (id_flux_x_ndiffuse_int_z(n), -1.0*T_prog(n)%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1),                                         &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif
     if(id_flux_y_ndiffuse_int_z(n) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) + flux_y(n)%field(isc:iec,jsc:jec,k)
        enddo        
        used = send_data (id_flux_y_ndiffuse_int_z(n), -1.0*T_prog(n)%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1),                                         &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif
     
  enddo   ! n=1,num_prog_tracers


  if(id_neutral_rho_ndiffuse > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*                     &
                    (drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k)  &
                    +drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k))
            enddo
         enddo
      enddo
      used = send_data (id_neutral_rho_ndiffuse, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),         &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_wdian_rho_ndiffuse > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               temporary   = Dens%drhodz_zt(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
               wrk1(i,j,k) = Grd%tmask(i,j,k)*                              &
                             (drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k)  &
                             +drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)) &
                             /(epsln+temporary)
            enddo
         enddo
      enddo
      used = send_data (id_wdian_rho_ndiffuse, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),       &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif



  if (id_k33_explicit > 0) used = send_data (id_k33_explicit, rho0r*K33_explicit(:,:,:), &
                                  Time%model_time, rmask=Grd%tmask(:,:,:),               &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  call mpp_clock_end(id_clock_ndiffuse)

end subroutine compute_ndiffusion
! </SUBROUTINE> NAME="compute_ndiffusion"


!#######################################################################
! <SUBROUTINE NAME="compute_gmskewsion">
!
! <DESCRIPTION>
! Subroutine to compute the tendency from GM skewsion, as determined
! by projecting GM streamfunction onto baroclinic modes.  
! </DESCRIPTION>
!
subroutine compute_gmskewsion(Time, Thickness, Dens, T_prog)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  real, dimension(isd:ied,jsd:jed) :: tmp_flux
  real    :: temporary 
  real    :: upsilonx, upsilony
  real    :: d_upsilonx_dz, d_upsilony_dz
  real    :: term1, term2, term3
  integer :: i,j,k,n,tau,kp1

  call mpp_clock_begin(id_clock_gm_modes)

  tau=Time%tau 

  do n=1,num_prog_tracers  
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0
    T_prog(n)%wrk1(:,:,:)  = 0.0
  enddo 
  
  call fx_flux_gm
  call fy_flux_gm
  if (Grd%tripolar) then 
      do n=1,num_prog_tracers
         call mpp_update_domains(flux_x(n)%field(:,:,:), flux_y(n)%field(:,:,:), Dom_flux%domain2d, &
              gridtype=CGRID_NE, complete=T_prog(n)%complete) 
      enddo
  endif

  do k=1,nk
     call mpp_clock_begin(id_clock_fz_flux_gm)
     call fz_flux_gm(T_prog,k)
     call mpp_clock_end(id_clock_fz_flux_gm)
     do n=1,num_prog_tracers
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%wrk1(i,j,k) =                                                 & 
                   Grd%tmask(i,j,k)                                                   &
                   *(fz1(n)%field(i,j)-fz2(n)%field(i,j)                              &
                   +( flux_x(n)%field(i,j,k)-flux_x(n)%field(i-1,j,k)                 &
                     +flux_y(n)%field(i,j,k)-flux_y(n)%field(i,j-1,k) )*Grd%datr(i,j) &
                    )
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)
              fz1(n)%field(i,j) = fz2(n)%field(i,j)
           enddo
        enddo

     enddo
  enddo


  ! diagnostics 
  do n=1,num_prog_tracers

     if(id_neutral_gm_modes(n) > 0) then 
        used = send_data (id_neutral_gm_modes(n), T_prog(n)%wrk1(:,:,:)*T_prog(n)%conversion, &
               Time%model_time, rmask=Grd%tmask(:,:,:),                                       &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

     ! send fluxes to diag_manager 
     ! minus sign is due to a MOM-convention for physics fluxes  
     if(id_flux_x_gm_modes(n) > 0) then 
        used = send_data (id_flux_x_gm_modes(n), -1.0*T_prog(n)%conversion*flux_x(n)%field(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:),                                            &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_flux_y_gm_modes(n) > 0) then 
        used = send_data (id_flux_y_gm_modes(n), -1.0*T_prog(n)%conversion*flux_y(n)%field(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:),                                            &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_flux_x_gm_modes_int_z(n) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) + flux_x(n)%field(isc:iec,jsc:jec,k)
        enddo        
        used = send_data (id_flux_x_gm_modes_int_z(n), -1.0*T_prog(n)%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1),                                         &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif
     if(id_flux_y_gm_modes_int_z(n) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) + flux_y(n)%field(isc:iec,jsc:jec,k)
        enddo        
        used = send_data (id_flux_y_gm_modes_int_z(n), -1.0*T_prog(n)%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1),                                         &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif
     
  enddo   ! n=1,num_prog_tracers


  if(id_neutral_rho_gm_modes > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = Grd%tmask(i,j,k)*                     &
                    (drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k)  &
                    +drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k))
            enddo
         enddo
      enddo
      used = send_data (id_neutral_rho_gm_modes, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),         &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_wdian_rho_gm_modes > 0) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               temporary   = Dens%drhodz_zt(i,j,k)*Thickness%rho_dzt(i,j,k,tau)
               wrk1(i,j,k) = Grd%tmask(i,j,k)*                              &
                             (drhodT(i,j,k)*T_prog(index_temp)%wrk1(i,j,k)  &
                             +drhodS(i,j,k)*T_prog(index_salt)%wrk1(i,j,k)) &
                             /(epsln+temporary)
            enddo
         enddo
      enddo
      used = send_data (id_wdian_rho_gm_modes, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:),       &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  if(id_gm_eddy_ke_source > 0) then 
      wrk1_2d(:,:) = 0.0
      do k=1,nk
         kp1 = min(nk,k+1) 
         do j=jsc,jec
            do i=isc,iec
               if(agm_array(i,j,k) > 0.0) then 
                   upsilonx = -onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))
                   upsilony =  onehalf*(psix(i,j,k,0)+psix(i,j,k,1))
                   d_upsilonx_dz = -onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1)-psiy(i,j,kp1,0)-psiy(i,j,kp1,1)) &
                        /(epsln+Thickness%dzt(i,j,k))
                   d_upsilony_dz =  onehalf*(psix(i,j,k,0)+psix(i,j,k,1)-psix(i,j,kp1,0)-psix(i,j,kp1,1)) &
                        /(epsln+Thickness%dzt(i,j,k))
                   term1 = Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)*Thickness%rho_dzt(i,j,k,tau)/agm_array(i,j,k)
                   term2 = bv_freq(i,j,k)*bv_freq(i,j,k)*(upsilonx**2 + upsilony**2)
                   term3 = bc_speed2(i,j)*(d_upsilonx_dz**2 + d_upsilony_dz**2)
                   wrk1_2d(i,j) = wrk1_2d(i,j) + term1*(term2 + term3) 
               endif
            enddo
         enddo
      enddo
      used = send_data (id_gm_eddy_ke_source, wrk1_2d(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1),       &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif



  call mpp_clock_end(id_clock_gm_modes)

end subroutine compute_gmskewsion
! </SUBROUTINE> NAME="compute_gmskewsion"



!#######################################################################
! <SUBROUTINE NAME="baroclinic_modes">
!
! <DESCRIPTION>
! Subroutine computes the baroclinic wave speeds and the dimensionless
! baroclinic mode eigenfunction for the vertical velocity baroclinic 
! modes.  These modes vanish at the surface and the bottom.  We use
! the Chelton etal WKB analytic formulae for the speeds and modes.
!
! The baroclinic modes are dimensionless, and normalized over the 
! depth of the ocean, from free surface to bottom.  
!
! The speeds are m/sec.
!
! </DESCRIPTION>
!
subroutine baroclinic_modes(Time, Thickness, Dens)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_density_type),   intent(in) :: Dens

  integer :: i,j,k,m,n,kbot
  real    :: drhodz_prev, tmpdrhodz
  real    :: norm, argument 
  real    :: bc_mode_prev, tmp_bcmode

  call mpp_clock_begin(id_clock_bc_modes)


  ! place buoyancy frequency at T-point into bv_freq.
  ! set minimum value for frequency.  
  bv_freq=0.0
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           bv_freq(i,j,k) = max(epsln_bv_freq,sqrt(abs(grav*rho0r*Dens%drhodz_zt(i,j,k))))
        enddo
     enddo
  enddo

  ! vertically smooth buoyancy frequency, even beyond that in ocean_density_mod
  if(bv_freq_smooth_vert) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied
               drhodz_prev = onefourth*bv_freq(i,j,1) 
               kbot=Grd%kmt(i,j)
               if(kbot>1) then
                   do k=2,kbot-2
                      tmpdrhodz      = bv_freq(i,j,k)
                      bv_freq(i,j,k) = drhodz_prev + onehalf*bv_freq(i,j,k) + onefourth*bv_freq(i,j,k+1)
                      drhodz_prev    = onefourth*tmpdrhodz
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! compute first baroclinic wave speed 
  bc_speed(:,:,:) = 0.0
  n=1
  do j=jsd,jed
     do i=isd,ied
        kbot = Grd%kmt(i,j)  
        if(kbot > 1) then 
            do k=1,kbot  
               bc_speed(i,j,n) = bc_speed(i,j,n) + Thickness%dzt(i,j,k)*bv_freq(i,j,k) 
            enddo
            bc_speed(i,j,n) = max(min_bc_speed,bc_speed(i,j,n)*pi_r)
        endif
     enddo
  enddo

  ! compute argument of sine used in defining the baroclinic modes 
  wrk2(:,:,:) = 0.0
  do j=jsd,jed
     do i=isd,ied
        kbot = Grd%kmt(i,j)  
        if(kbot > 1) then 
            wrk2(i,j,kbot) = bv_freq(i,j,kbot) &
             *(Thickness%depth_zwt(i,j,kbot)-Thickness%depth_zt(i,j,kbot)) 
            do k=kbot-1,1,-1
               wrk2(i,j,k) = wrk2(i,j,k+1)                                                  &
               + (bv_freq(i,j,k+1)*(Thickness%depth_zt(i,j,k+1)-Thickness%depth_zwt(i,j,k)) &
               +  bv_freq(i,j,k)  *(Thickness%depth_zwt(i,j,k) -Thickness%depth_zt(i,j,k))) 
            enddo
            do k=1,kbot
               wrk2(i,j,k) = wrk2(i,j,k)/(epsln + bc_speed(i,j,1))  
            enddo
        endif
     enddo
  enddo

  ! compute baroclinic modes 
  do n=1,number_bc_modes

     do j=jsd,jed
        do i=isd,ied

           kbot = Grd%kmt(i,j)  
           bc_mode(i,j,:,n) = 0.0

           if(kbot > 1) then 

               ! wave speeds 
               bc_speed(i,j,n) = bc_speed(i,j,1)/n 

               ! dimensionless baroclinic mode
               do k=1,kbot
                  bc_mode(i,j,k,n) = sin(wrk2(i,j,k)*n)  
               enddo

               ! vertically smooth bc_mode prior to normalization 
               if(smooth_bc_modes) then
                   do m=1,num_121_passes
                      bc_mode_prev        = onefourth*bc_mode(i,j,1,n) 
                      do k=2,kbot-2
                         tmp_bcmode       = bc_mode(i,j,k,n) 
                         bc_mode(i,j,k,n) = bc_mode_prev + onehalf*bc_mode(i,j,k,n) + onefourth*bc_mode(i,j,k+1,n)
                         bc_mode_prev     = onefourth*tmp_bcmode
                      enddo
                   enddo
               endif

               ! normalize the baroclinic mode 
               norm = 0.0
               do k=1,kbot
                  norm = norm + Thickness%dzt(i,j,k)*(bv_freq(i,j,k)*bc_mode(i,j,k,n))**2
               enddo
               norm = sqrt_grav/(epsln+sqrt(norm))
               do k=1,kbot
                  bc_mode(i,j,k,n) = bc_mode(i,j,k,n)*norm
               enddo

           endif

        enddo
     enddo
  enddo

  ! for use in BVP approach 
  if(bvp_bc_mode == 0 .or. bvp_constant_speed) then 
      bc_speed2(:,:) = bvp_speed*bvp_speed*Grd%tmask(:,:,1)
  else
      do j=jsd,jed
         do i=isd,ied
            bc_speed2(i,j) = Grd%tmask(i,j,1) &
                *max(bvp_min_speed**2,(bvp_speed-bc_speed(i,j,bvp_bc_mode))**2)
         enddo
      enddo
  endif
 
  ! to avoid problems with spurious init values 
  if(Time%itt0 <= 1) bc_mode=0.0


  ! diagnostics 
  if(id_bv_freq > 0) then 
      used = send_data (id_bv_freq, bv_freq(:,:,:),&
           Time%model_time, rmask=Grd%tmask(:,:,:),&
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk) 
  endif 

  do n=1,number_bc_modes
     if(id_bc_mode(n) > 0) then 
        used = send_data (id_bc_mode(n), bc_mode(:,:,:,n), &
               Time%model_time, rmask=Grd%tmask(:,:,:),    &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_bc_speed(n) > 0) then 
        used = send_data (id_bc_speed(n), bc_speed(:,:,n), &
               Time%model_time, rmask=Grd%tmask(:,:,1),    &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif
  enddo

  call mpp_clock_end(id_clock_bc_modes)


end subroutine baroclinic_modes
! </SUBROUTINE> NAME="baroclinic_modes"


!#######################################################################
! <SUBROUTINE NAME="compute_psi_modes">
!
! <DESCRIPTION>
! Compute vector streamfunction as projection onto baroclinic modes.
!
! Units of psi are m^2/sec
!
! </DESCRIPTION>
!
subroutine compute_psi_modes(Time, Dens, Thickness, T_prog)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  real    :: factor1, factor2
  real    :: gradx, grady 
  real    :: active_cells
  real    :: project_x(0:1), project_y(0:1)
  real    :: tmpx, tmpy
  integer :: i, j, k, n, kbot 
  integer :: ip, jq
  integer :: tau

  call mpp_clock_begin(id_clock_psi_modes)

  tau = Time%tau

  psix(:,:,:,:)   = 0.0
  psiy(:,:,:,:)   = 0.0

  ! compute horizontal density gradients 
  ! compute max(agm_array)
  wrk1_v(:,:,:,:) = 0.0
  wrk2_v(:,:,:,:) = 0.0
  wrk3_v(:,:,:,:) = 0.0
  wrk4_v(:,:,:,:) = 0.0
  wrk1_2d(:,:)    = 0.0
  do k=1,nk
     do j=jsc-1,jec
        do i=isc-1,iec   
           do ip=0,1 
              jq=ip
              gradx = drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k) &
                     +drhodS(i+ip,j,k)*dTdx(index_salt)%field(i,j,k) 
              grady = drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k) &
                     +drhodS(i,j+jq,k)*dTdy(index_salt)%field(i,j,k)
              wrk1_v(i,j,k,ip+1) = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)*gradx
              wrk2_v(i,j,k,jq+1) = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)*grady
              wrk3_v(i,j,k,ip+1) = Thickness%dzt(i,j,k)*wrk1_v(i,j,k,ip+1)*agm_array(i,j,k) 
              wrk4_v(i,j,k,jq+1) = Thickness%dzt(i,j,k)*wrk2_v(i,j,k,jq+1)*agm_array(i,j,k) 
           enddo
           if(agm_array(i,j,k) > wrk1_2d(i,j)) wrk1_2d(i,j) = agm_array(i,j,k)          
        enddo
     enddo
  enddo

  ! for psi regularization 
  wrk1_2d(:,:) = wrk1_2d(:,:)*smax_psi


  do n=1,number_bc_modes

     do j=jsc-1,jec
        do i=isc-1,iec

           kbot = Grd%kmt(i,j)
           if(kbot > 1) then  

               ! projection of agm*grad_rho onto baroclinic modes, 
               ! integrated over the depth of an ocean column. 
               ! project_x and project_y have dimensions kg/(m*sec) 
               project_x(:)=0.0 
               project_y(:)=0.0 
               do k=1,kbot 
                  do ip=0,1 
                     jq=ip
                     project_x(ip) = project_x(ip) + bc_mode(i,j,k,n)*wrk3_v(i,j,k,ip+1)
                     project_y(jq) = project_y(jq) + bc_mode(i,j,k,n)*wrk4_v(i,j,k,jq+1)
                  enddo
               enddo

               ! compute the vector streamfunction (m^2/sec)
               do ip=0,1 
                  jq=ip
                  do k=1,kbot
                     psix(i,j,k,jq) = psix(i,j,k,jq) - rho0r*bc_mode(i,j,k,n)*project_y(jq)
                     psiy(i,j,k,ip) = psiy(i,j,k,ip) + rho0r*bc_mode(i,j,k,n)*project_x(ip)
                  enddo
               enddo  ! do ip=0,1

           endif   !if(kbot > 1)

        enddo   ! do j
     enddo      ! do i

  enddo   !n=1,number_bc_modes


  ! smooth psi to reduce potentials for checkerboard noise 
  if(smooth_psi) then 

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

      do ip=0,1
         jq=ip
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec

                  if(Grd%tmask(i,j,k)==1.0) then 

                      active_cells = 4.0       +&
                           Grd%tmask(i-1,j,k)  +&
                           Grd%tmask(i+1,j,k)  +&
                           Grd%tmask(i,j-1,k)  +&
                           Grd%tmask(i,j+1,k)

                      if (active_cells > 4.0) then
                          wrk1(i,j,k) = & 
                               (4.0*psix(i,j,k,jq) +&
                               psix(i-1,j,k,jq)    +&
                               psix(i+1,j,k,jq)    +&
                               psix(i,j-1,k,jq)    +&
                               psix(i,j+1,k,jq)) / active_cells
                          wrk2(i,j,k) =  &
                               (4.0*psiy(i,j,k,ip) +&
                               psiy(i-1,j,k,ip)    +&
                               psiy(i+1,j,k,ip)    +&
                               psiy(i,j-1,k,ip)    +&
                               psiy(i,j+1,k,ip)) / active_cells

                      else
                          wrk1(i,j,k) = psix(i,j,k,jq)
                          wrk2(i,j,k) = psiy(i,j,k,ip)
                      endif

                  endif

               enddo
            enddo

            do j=jsc,jec
               do i=isc,iec
                  psix(i,j,k,jq) = wrk1(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                  psiy(i,j,k,ip) = wrk2(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
               enddo
            enddo

         enddo
      enddo

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

  endif

  ! regularize the vector streamfunction to keep magnitude 
  ! under control in regions of weak vertical stratification.
  wrk1_v2d(:,:,:) = 1.0
  wrk2_v2d(:,:,:) = 1.0 
  if(regularize_psi) then 

      do ip=0,1
         jq=ip
         do j=jsc-1,jec
            do i=isc-1,iec

               kbot=Grd%kmt(i,j)
               if(kbot > 1) then 

                   wrk1_v2d(i,j,ip+1) = wrk1_2d(i,j)
                   wrk2_v2d(i,j,jq+1) = wrk1_2d(i,j)
                   do k=1,kbot
                      tmpx = abs(psix(i,j,k,jq))
                      tmpy = abs(psiy(i,j,k,ip))
                      if(tmpx > wrk1_v2d(i,j,jq+1)) wrk1_v2d(i,j,jq+1) = tmpx
                      if(tmpy > wrk2_v2d(i,j,ip+1)) wrk2_v2d(i,j,ip+1) = tmpy
                   enddo
                   wrk1_v2d(i,j,jq+1) = wrk1_2d(i,j)/wrk1_v2d(i,j,jq+1)
                   wrk2_v2d(i,j,ip+1) = wrk1_2d(i,j)/wrk2_v2d(i,j,ip+1)

                   do k=1,kbot
                     psix(i,j,k,jq) = psix(i,j,k,jq)*wrk1_v2d(i,j,jq+1) 
                     psiy(i,j,k,ip) = psiy(i,j,k,ip)*wrk2_v2d(i,j,ip+1) 
                   enddo

               endif

            enddo
         enddo
      enddo

  endif

  if(id_rescale_psi_x > 0) then 
      used = send_data (id_rescale_psi_x, onehalf*(wrk1_v2d(:,:,1)+wrk1_v2d(:,:,2)), &
           Time%model_time, rmask=Grd%tmask(:,:,1),                                  &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  if(id_rescale_psi_y > 0) then 
      used = send_data (id_rescale_psi_y, onehalf*(wrk2_v2d(:,:,1)+wrk2_v2d(:,:,2)), &
           Time%model_time, rmask=Grd%tmask(:,:,1),                                  &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  if (id_psix_gm_modes > 0) then 
      used = send_data (id_psix_gm_modes, onehalf*(psix(:,:,:,0)+psix(:,:,:,1)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                           &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_psiy_gm_modes > 0) then 
      used = send_data (id_psiy_gm_modes, onehalf*(psiy(:,:,:,0)+psiy(:,:,:,1)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                           &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 

  if(id_gradx_rho > 0) then 
      used = send_data (id_gradx_rho, onehalf*(wrk1_v(:,:,:,1)+wrk1_v(:,:,:,2)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                            &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_grady_rho > 0) then 
      used = send_data (id_grady_rho, onehalf*(wrk2_v(:,:,:,1)+wrk2_v(:,:,:,2)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                            &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


  ! diagnose mass transports for sending to diagnostic
  wrk1_v = 0.0
  if(vert_coordinate_class==DEPTH_BASED) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = &
                -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*rho0
               wrk1_v(i,j,k,2) = &
                 transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*rho0
            enddo
         enddo
      enddo
  else
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = &
                -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*Dens%rho(i,j,k,tau)
               wrk1_v(i,j,k,2) = &
                 transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*Dens%rho(i,j,k,tau)
            enddo
         enddo
      enddo
  endif

  ! change signs to agree with convention in nphysicsA and nphysicsB
  call transport_on_rho_gm (Time, Dens, &
       -1.0*wrk1_v(:,:,:,1), -1.0*wrk1_v(:,:,:,2))
  call transport_on_nrho_gm(Time, Dens, &
       -1.0*wrk1_v(:,:,:,1), -1.0*wrk1_v(:,:,:,2))
  call transport_on_theta_gm(Time, Dens, T_prog(index_temp), &
       -1.0*wrk1_v(:,:,:,1), -1.0*wrk1_v(:,:,:,2))

  ! change signs to agree with convention in nphysicsA and nphysicsB
  if (id_tx_trans_gm > 0) then 
      used = send_data (id_tx_trans_gm, -1.0*wrk1_v(:,:,:,1), &
           Time%model_time, rmask=Grd%tmask(:,:,:),           &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_ty_trans_gm > 0) then 
      used = send_data (id_ty_trans_gm, -1.0*wrk1_v(:,:,:,2), &
           Time%model_time, rmask=Grd%tmask(:,:,:),           &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


  call mpp_clock_end(id_clock_psi_modes)


end subroutine compute_psi_modes
! </SUBROUTINE> NAME="compute_psi_modes"



!#######################################################################
! <SUBROUTINE NAME="compute_psi_bvp">
!
! <DESCRIPTION>
! Compute vector streamfunction by solving a boundary value problem. 
!
! psi is centered on bottom of tracer cell; for example, 
! psi(k=1)=psi at bottom of tracer cell k=1.
! psi vanishes at the ocean surface: psi(k=0)=0
! and ocean bottom: psi(k=kmt)=0.
!
! We solve for psi(k=1,kmt-1) using a tridiagonal solver from 
! Section 2.4 of Press etal 1986.
!
! Units of psi are m^2/sec
!
! </DESCRIPTION>
!
subroutine compute_psi_bvp(Time, Dens, Thickness, T_prog)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_density_type),     intent(in) :: Dens
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  real    :: factor1, factor2
  real    :: gradx, grady 
  real    :: active_cells
  real    :: project_x(0:1), project_y(0:1)
  real    :: tmpx, tmpy
  real    :: dzwtr
  integer :: i, j, k, n, kp1, kbot 
  integer :: ip, jq
  integer :: tau

  call mpp_clock_begin(id_clock_psi_bvp)

  tau = Time%tau

  ! initialize to zero 
  psix(:,:,:,:) = 0.0
  psiy(:,:,:,:) = 0.0


  ! compute horizontal density gradients and RHS source terms 
  wrk1_v(:,:,:,:) = 0.0
  wrk2_v(:,:,:,:) = 0.0
  wrk3_v(:,:,:,:) = 0.0
  wrk4_v(:,:,:,:) = 0.0
  wrk1_2d(:,:)    = 0.0
  do k=1,nk-1
     do j=jsc-1,jec
        do i=isc-1,iec   
           do ip=0,1 
              jq=ip
              gradx = drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k) &
                     +drhodS(i+ip,j,k)*dTdx(index_salt)%field(i,j,k) 
              grady = drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k) &
                     +drhodS(i,j+jq,k)*dTdy(index_salt)%field(i,j,k)
              wrk1_v(i,j,k,ip+1) = Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)*gradx
              wrk2_v(i,j,k,jq+1) = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)*grady
              wrk3_v(i,j,k,ip+1) = Grd%tmask(i,j,k+1)*grav_rho0r*wrk1_v(i,j,k,ip+1)*agm_array(i,j,k) 
              wrk4_v(i,j,k,jq+1) = Grd%tmask(i,j,k+1)*grav_rho0r*wrk2_v(i,j,k,jq+1)*agm_array(i,j,k) 
           enddo
        enddo
     enddo
  enddo

  ! intermediate diagnostics 
  if(id_gradx_rho > 0) then 
      used = send_data (id_gradx_rho, onehalf*(wrk1_v(:,:,:,1)+wrk1_v(:,:,:,2)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                            &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_grady_rho > 0) then 
      used = send_data (id_grady_rho, onehalf*(wrk2_v(:,:,:,1)+wrk2_v(:,:,:,2)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                            &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


  ! a,b,c fields (using Press etal nomenclature) to invert tridiagonal matrix 
  wrk1(:,:,:) = 0.0   ! a
  wrk2(:,:,:) = 0.0   ! b
  wrk3(:,:,:) = 0.0   ! c 
  do k=1,nk-1
     do j=jsc-1,jec
        do i=isc-1,iec    
           dzwtr       = Grd%tmask(i,j,k)*Grd%tmask(i,j,k+1)/(epsln+Thickness%dzwt(i,j,k))          
           wrk1(i,j,k) = bc_speed2(i,j)*dzwtr/(epsln+Thickness%dzt(i,j,k))
           wrk3(i,j,k) = bc_speed2(i,j)*dzwtr/(epsln+Thickness%dzt(i,j,k+1))
           wrk2(i,j,k) = -Grd%tmask(i,j,k)*Grd%tmask(i,j,k+1) &
                          *(bv_freq(i,j,k)*bv_freq(i,j,k) + wrk1(i,j,k) + wrk3(i,j,k))
        enddo
     enddo
  enddo


  ! invert to solve for upsilon_x, and set psiy=-upsilon_x
  wrk1_v = 0.0
  call invtri_bvp(Grd%tmask, wrk1, wrk2, wrk3, wrk3_v, wrk1_v)
  psiy(:,:,:,:) = -wrk1_v(:,:,:,:)    

  ! invert to solve for upsilon_y, and set psix=upsilon_y
  wrk2_v = 0.0
  call invtri_bvp(Grd%tmask, wrk1, wrk2, wrk3, wrk4_v, wrk2_v)    
  psix(:,:,:,:) = wrk2_v(:,:,:,:)    

  ! smooth psi to reduce potentials for checkerboard noise 
  if(smooth_psi) then 

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

      do ip=0,1
         jq=ip
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec

                  if(Grd%tmask(i,j,k)==1.0) then 

                      active_cells = 4.0       +&
                           Grd%tmask(i-1,j,k)  +&
                           Grd%tmask(i+1,j,k)  +&
                           Grd%tmask(i,j-1,k)  +&
                           Grd%tmask(i,j+1,k)

                      if (active_cells > 4.0) then
                          wrk1(i,j,k) = & 
                               (4.0*psix(i,j,k,jq) +&
                               psix(i-1,j,k,jq)    +&
                               psix(i+1,j,k,jq)    +&
                               psix(i,j-1,k,jq)    +&
                               psix(i,j+1,k,jq)) / active_cells
                          wrk2(i,j,k) =  &
                               (4.0*psiy(i,j,k,ip) +&
                               psiy(i-1,j,k,ip)    +&
                               psiy(i+1,j,k,ip)    +&
                               psiy(i,j-1,k,ip)    +&
                               psiy(i,j+1,k,ip)) / active_cells

                      else
                          wrk1(i,j,k) = psix(i,j,k,jq)
                          wrk2(i,j,k) = psiy(i,j,k,ip)
                      endif

                  endif

               enddo
            enddo

            do j=jsc,jec
               do i=isc,iec
                  psix(i,j,k,jq) = wrk1(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
                  psiy(i,j,k,ip) = wrk2(i,j,k)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)
               enddo
            enddo

         enddo
      enddo

      call mpp_update_domains(psix(:,:,:,:), Dom%domain2d) 
      call mpp_update_domains(psiy(:,:,:,:), Dom%domain2d) 

  endif


  ! diagnostics 
  if (id_psix_gm_bvp > 0) then 
      used = send_data (id_psix_gm_bvp, onehalf*(psix(:,:,:,0)+psix(:,:,:,1)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                           &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_psiy_gm_bvp > 0) then 
      used = send_data (id_psiy_gm_bvp, onehalf*(psiy(:,:,:,0)+psiy(:,:,:,1)), &
             Time%model_time, rmask=Grd%tmask(:,:,:),                           &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 



  ! diagnose mass transports for sending to diagnostic
  wrk1_v = 0.0
  if(vert_coordinate_class==DEPTH_BASED) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = &
                -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*rho0
               wrk1_v(i,j,k,2) = &
                 transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*rho0
            enddo
         enddo
      enddo
  else
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,1) = &
                -transport_convert*onehalf*(psiy(i,j,k,0)+psiy(i,j,k,1))*Grd%dyte(i,j)*Dens%rho(i,j,k,tau)
               wrk1_v(i,j,k,2) = &
                 transport_convert*onehalf*(psix(i,j,k,0)+psix(i,j,k,1))*Grd%dxtn(i,j)*Dens%rho(i,j,k,tau)
            enddo
         enddo
      enddo
  endif

  ! change signs to agree with convention in nphysicsA and nphysicsB
  call transport_on_rho_gm (Time, Dens, &
       -1.0*wrk1_v(:,:,:,1), -1.0*wrk1_v(:,:,:,2))
  call transport_on_nrho_gm(Time, Dens, &
       -1.0*wrk1_v(:,:,:,1), -1.0*wrk1_v(:,:,:,2))
  call transport_on_theta_gm(Time, Dens, T_prog(index_temp), &
       -1.0*wrk1_v(:,:,:,1), -1.0*wrk1_v(:,:,:,2))


  ! change signs to agree with convention in nphysicsA and nphysicsB
  if (id_tx_trans_gm > 0) then 
      used = send_data (id_tx_trans_gm, -1.0*wrk1_v(:,:,:,1), &
           Time%model_time, rmask=Grd%tmask(:,:,:),           &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_ty_trans_gm > 0) then 
      used = send_data (id_ty_trans_gm, -1.0*wrk1_v(:,:,:,2), &
           Time%model_time, rmask=Grd%tmask(:,:,:),           &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


  call mpp_clock_end(id_clock_psi_bvp)


end subroutine compute_psi_bvp
! </SUBROUTINE> NAME="compute_psi_bvp"


!#######################################################################
! <SUBROUTINE NAME="fz_terms">
!
! <DESCRIPTION>
! Subroutine computes the tracer independent pieces of the vertical 
! flux component. As a result of this routine, 
! Array tensor_31 = x-diffusivity*slope (m^2/sec) for fz
! Array tensor_32 = y-diffusivity*slope (m^2/sec) for fz 
!
! K33 is the (3,3) term in small angle Redi diffusion tensor.
! It is broken into an explicit in time piece and implicit 
! in time piece.  It is weighted by density for non-Boussinesq
! and rho0 for Boussinesq.  
!
! K33 has units (kg/m^3)*m^2/sec.
!
! Also will compute the squared Eady growth rate, with the maximum
! slope contributing to this growth rate set by smax.
!
! This routine is nearly the same as in ocean_nphysicsA, except 
! for the following changes: 
! 1/ the routine here removes all pieces related to GM-skewsion.  
! 2/ the routine here uses Thickness%depth_zwt rather than Grd%zt. 
!
! </DESCRIPTION>
!
subroutine fz_terms(Time, Thickness, T_prog, rho)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(num_prog_tracers)
  real, dimension(isd:,jsd:,:), intent(in)    :: rho

  integer :: i, j, k, n, ip, jq, kr, tau

  real :: baroclinic_triad, K33, K33_crit
  real :: sumx(0:1), sumy(0:1)

  real :: slope, absslope, depth_ratio
  real :: taperA, taperB
  real :: taper_slope, taper_slope2

  real :: aredi_scalar
  real :: absdrhodzb(0:1)
  real :: mindelqc(0:1,0:1), delqc_ijk_1, delqc_ijkp1_0

  call mpp_clock_begin(id_clock_fz_terms)

  tau = Time%tau

  eady_termx(:,:,:)     = 0.0
  eady_termy(:,:,:)     = 0.0
  baroclinic_termx(:,:) = 0.0
  baroclinic_termy(:,:) = 0.0

  ! Main loop. Short ip,jq,kr loops are explicitly unrolled
  ! in order to expose independence and allow vectorisation

  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec

           aredi_scalar  = aredi_array(i,j,k)

           absdrhodzb(0) = abs(drhodzb(i,j,k,0))
           absdrhodzb(1) = abs(drhodzb(i,j,k,1))

           delqc_ijk_1      = delqc(i,j,k,1)
           delqc_ijkp1_0    = delqc(i,j,k+1,0)

           !mindelqc(ip,kr) = min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr)) 
           ip=0 ; kr=0
           mindelqc(ip,kr)  = min(delqc(i-1,j,k,1),delqc_ijk_1) 
           ip=0 ; kr=1
           mindelqc(ip,kr)  = min(delqc(i-1,j,k+1,0),delqc_ijkp1_0) 
           ip=1 ; kr=0
           mindelqc(ip,kr)  = min(delqc_ijk_1,delqc(i+1,j,k,1)) 
           ip=1 ; kr=1
           mindelqc(ip,kr)  = min(delqc_ijkp1_0,delqc(i+1,j,k+1,0)) 

           baroclinic_triad = 0.0

           ip=0   
           sumx(ip) = 0.0

           kr=0
           slope    = tensor_31(i,j,k,ip,kr)
           absslope = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope
           sumx(ip)               = sumx(ip)         + mindelqc(ip,kr)*taper_slope2
           eady_termx(i,j,k)      = eady_termx(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)

           kr=1
           slope    = tensor_31(i,j,k,ip,kr)
           absslope = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope
           sumx(ip)               = sumx(ip)         + mindelqc(ip,kr)*taper_slope2
           eady_termx(i,j,k)      = eady_termx(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)
           sumx(ip)               = sumx(ip)*dtew(i,j,ip)

           ip=1   
           sumx(ip) = 0.0

           kr=0
           slope        = tensor_31(i,j,k,ip,kr)
           absslope     = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope
           sumx(ip)               = sumx(ip)         + mindelqc(ip,kr)*taper_slope2
           eady_termx(i,j,k)      = eady_termx(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)

           kr=1
           slope        = tensor_31(i,j,k,ip,kr)
           absslope     = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope
           sumx(ip)               = sumx(ip)         + mindelqc(ip,kr)*taper_slope2
           eady_termx(i,j,k)      = eady_termx(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)
           sumx(ip)               = sumx(ip)*dtew(i,j,ip)

           if(Thickness%depth_zwt(i,j,k) >= agm_closure_upper_depth .and. &
              Thickness%depth_zwt(i,j,k) <= agm_closure_lower_depth) then 
              baroclinic_termx(i,j) = baroclinic_termx(i,j) + baroclinic_triad*Thickness%dzwt(i,j,k)
           endif
   
           !mindelqc(jq,kr) = min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr)) 
           jq=0 ; kr=0
           mindelqc(jq,kr)  = min(delqc(i,j-1,k,1),delqc_ijk_1) 
           jq=0 ; kr=1
           mindelqc(jq,kr)  = min(delqc(i,j-1,k+1,0),delqc_ijkp1_0) 
           jq=1 ; kr=0
           mindelqc(jq,kr)  = min(delqc_ijk_1,delqc(i,j+1,k,1)) 
           jq=1 ; kr=1
           mindelqc(jq,kr)  = min(delqc_ijkp1_0,delqc(i,j+1,k+1,0)) 

           baroclinic_triad = 0.0 

           jq=0   
           sumy(jq) = 0.0

           kr=0
           slope      = tensor_32(i,j,k,jq,kr)
           absslope   = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope
           sumy(jq)               = sumy(jq)         + mindelqc(jq,kr)*taper_slope2
           eady_termy(i,j,k)      = eady_termy(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)

           kr=1
           slope        = tensor_32(i,j,k,jq,kr)
           absslope     = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope
           sumy(jq)               = sumy(jq)         + mindelqc(jq,kr)*taper_slope2
           eady_termy(i,j,k)      = eady_termy(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)
           sumy(jq)               = sumy(jq)*dtns(i,j,jq)

           jq=1   
           sumy(jq) = 0.0

           kr=0
           slope    = tensor_32(i,j,k,jq,kr)
           absslope = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope
           sumy(jq)               = sumy(jq)         + mindelqc(jq,kr)*taper_slope2
           eady_termy(i,j,k)      = eady_termy(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)

           kr=1
           slope    = tensor_32(i,j,k,jq,kr)
           absslope = abs(slope)

           ! taper for steep slope regions 
           if(absslope > smax) then 
               taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
               taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
           else 
               taperA = 1.0
           endif

           ! taper for grid depths shallower than eddy_depth
           if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
               depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
               taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
           else 
               taperB          = 1.0
           endif

           ! taper times slope for use with neutral diffusion
           taper_slope  = taperA*taperB*slope
           taper_slope2 = slope*taper_slope

           tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope
           sumy(jq)               = sumy(jq)         + mindelqc(jq,kr)*taper_slope2
           eady_termy(i,j,k)      = eady_termy(i,j,k)+ absdrhodzb(kr)*absslope*absslope
           baroclinic_triad       = baroclinic_triad + absdrhodzb(kr)*abs(taper_slope)

           sumy(jq) = sumy(jq)*dtns(i,j,jq)

           if(Thickness%depth_zwt(i,j,k) >= agm_closure_upper_depth .and. &
              Thickness%depth_zwt(i,j,k) <= agm_closure_lower_depth) then 
               baroclinic_termy(i,j) = baroclinic_termy(i,j) + baroclinic_triad*Thickness%dzwt(i,j,k)
           endif

           K33 = aredi_scalar*Grd%tmask(i,j,k+1)*(Grd%dxtr(i,j)*(sumx(0)+sumx(1)) &
                + Grd%dytr(i,j)*(sumy(0)+sumy(1)))*dzwtr(i,j,k)

           ! determine part of K33 included explicitly in time and that part implicitly. 
           ! explicit calculation is more accurate than implicit, so aim to compute as much 
           ! as stably possible via the explicit method. 
           K33_explicit(i,j,k) = K33
           K33_implicit(i,j,k) = 0.0 
           if(.not. diffusion_all_explicit) then 
               K33_crit = two_dtime_inv*Thickness%rho_dzt(i,j,k,tau)*Thickness%dzt(i,j,k)
               if(K33 >= K33_crit) then 
                   K33_explicit(i,j,k) = K33_crit
                   K33_implicit(i,j,k) = K33-K33_crit
               endif
           endif

        enddo
     enddo
  enddo

  if(tmask_neutral_on) then
      do j=jsc,jec
         do i=isc,iec
            K33_implicit(i,j,1) = 0.0  
            K33_explicit(i,j,1) = 0.0  
            if(Grd%kmt(i,j) > 1) then 
                k = Grd%kmt(i,j)-1
                K33_implicit(i,j,k) = 0.0
                K33_explicit(i,j,k) = 0.0
            endif
         enddo
      enddo
  endif

  do n=1,num_prog_tracers 
     T_prog(n)%K33_implicit(:,:,:) = K33_implicit(:,:,:)
  enddo
  if(neutral_physics_limit) then 
      do n=1,num_prog_tracers 
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if(T_prog(n)%tmask_limit(i,j,k)==1.0) T_prog(n)%K33_implicit(i,j,k) = 0.0
               enddo
            enddo
         enddo
      enddo
  endif


  do n=1,num_prog_tracers
     if (id_k33_implicit(n) > 0) then
         used = send_data (id_k33_implicit(n), rho0r*T_prog(n)%K33_implicit(:,:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,:), &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
  enddo

  if(id_N_squared > 0) then 
      wrk2(:,:,:)=0.0
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               wrk2(i,j,k) = 1.0/(epsln + rho(i,j,k)) 
               wrk2(i,j,k) = -onehalf*grav*(drhodzb(i,j,k,0)+drhodzb(i,j,k,1))*wrk2(i,j,k) 
            enddo
         enddo
      enddo
      used = send_data (id_N_squared, wrk2(:,:,:),  &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk) 
  endif


  call mpp_clock_end(id_clock_fz_terms)

end subroutine fz_terms
! </SUBROUTINE>  NAME="fz_terms"


!#######################################################################
! <SUBROUTINE NAME="fx_flux_ndiffuse">
!
! <DESCRIPTION>
! Subroutine computes the zonal neutral diffusion tracer flux component.
! Compute this component for all tracers at level k.
!
! fx has physical dimensions (area*diffusivity*density*tracer gradient)
!
! This routine is the same as that in ocean_nphysicsA, except 
! for the following changes: 
! 1/ the routine here removes all pieces related to GM-skewsion.  
! 2/ the routine here uses Thickness%depth_zwt rather than Grd%zt. 
! 3/ ah_array is removed.
!
! </DESCRIPTION>
!
subroutine fx_flux_ndiffuse(Time, Thickness, T_prog, k)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)
  integer,                      intent(in) :: k

  integer :: i, j, n, ip, tau
  integer :: kr, kpkr
  real    :: tensor_11(isd:ied,jsd:jed,0:1,0:1)
  real    :: tensor_13(isd:ied,jsd:jed,0:1,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)
  real    :: taperA, taperB
  real    :: slope, absslope, taper_slope
  real    :: depth_ratio
  real    :: drhodx, drhodz 

  call mpp_clock_begin(id_clock_fx_flux_ndiffuse)

  tau = Time%tau

  ! initialize arrays to zero 
  tensor_11(:,:,:,:) = 0.0
  tensor_13(:,:,:,:) = 0.0

  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do ip=0,1 
        do j=jsc,jec
           do i=isc-1,iec

              drhodz   = drhodzh(i+ip,j,k,kr)
              drhodx   = Grd%tmask(i+ip,j,kpkr)                            &
                        *( drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k)  &
                          +drhodS(i+ip,j,k)*dTdx(index_salt)%field(i,j,k))
              slope    = -drhodx/drhodz
              absslope = abs(slope)

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
                  depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! fill tensor components 
              tensor_11(i,j,ip,kr) = aredi_array(i,j,k) 
              tensor_13(i,j,ip,kr) = aredi_array(i,j,k)*taper_slope

           enddo
        enddo
     enddo
  enddo

  ! tracer-dependent part of the calculation
  do n=1,num_prog_tracers
     sumz(:,:,:) = 0.0
     do kr=0,1
        do ip=0,1
           do j=jsc,jec
              do i=isc-1,iec
                sumz(i,j,kr) = sumz(i,j,kr) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) &
                *(tensor_11(i,j,ip,kr)*dTdx(n)%field(i,j,k) + tensor_13(i,j,ip,kr)*dTdz(n)%field(i+ip,j,k-1+kr)) &
                * min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc,jec
        do i=isc-1,iec
           flux_x(n)%field(i,j,k) = Grd%dxter(i,j)*(sumz(i,j,0)+sumz(i,j,1))
        enddo
     enddo
  enddo

  ! apply some masks 
  if(tmask_neutral_on) then
     do n=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(k==Grd%kmt(i,j)) then
                flux_x(n)%field(i,j,k) = aredi_array(i,j,k)*dTdx(n)%field(i,j,k)*Grd%dyte(i,j)* &
                                          min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i+1,j,k,tau))
              endif
           enddo
        enddo
     enddo
     if(k==1) then
         do n=1,num_prog_tracers
            flux_x(n)%field(:,:,k) = aredi_array(:,:,k)*dTdx(n)%field(:,:,k)*Grd%dyte(:,:) &
                                      *FMX(Thickness%rho_dzt(:,:,k,tau))
         enddo
     endif
  endif

 if(neutral_physics_limit) then 
     do n=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(T_prog(n)%tmask_limit(i,j,k)==1.0) then 
                  flux_x(n)%field(i,j,k) = aredi_array(i,j,k)*dTdx(n)%field(i,j,k)*Grd%dyte(i,j) &
                                            *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i+1,j,k,tau))
              endif
           enddo
        enddo
     enddo
 endif

 call mpp_clock_end(id_clock_fx_flux_ndiffuse) 

end subroutine fx_flux_ndiffuse
! </SUBROUTINE> NAME="fx_flux_ndiffuse"



!#######################################################################
! <SUBROUTINE NAME="fy_flux_ndiffuse">
!
! <DESCRIPTION>
! Subroutine computes the meridional neutral diffusion tracer flux component.
! Compute this component for all tracers at level k.
!
! fy has physical dimensions (area*diffusivity*density*tracer gradient)
!
! This routine is the same as that in ocean_nphysicsA, except 
! for the following changes: 
! 1/ the routine here removes all pieces related to GM-skewsion.  
! 2/ the routine here uses Thickness%depth_zwt rather than Grd%zt. 
! 3/ ah_array is removed.
!
! </DESCRIPTION>
!
subroutine fy_flux_ndiffuse(Time, Thickness, T_prog, k)

  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)
  integer,                      intent(in) :: k

  integer :: i, j, n, jq, tau
  integer :: kr, kpkr
  real :: tensor_22(isd:ied,jsd:jed,0:1,0:1)
  real :: tensor_23(isd:ied,jsd:jed,0:1,0:1)
  real :: sumz(isd:ied,jsd:jed,0:1)
  real :: taperA, taperB
  real :: slope, absslope, taper_slope
  real :: depth_ratio
  real :: drhody, drhodz 

  call mpp_clock_begin(id_clock_fy_flux_ndiffuse)

  tau = Time%tau 

  ! initialize arrays to zero 
  tensor_22(:,:,:,:) = 0.0
  tensor_23(:,:,:,:) = 0.0

  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do jq=0,1  
       do j=jsc-1,jec
          do i=isc,iec
 
              drhodz   = drhodzh(i,j+jq,k,kr) 
              drhody   = Grd%tmask(i,j+jq,kpkr)                            &
                        *( drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k)  &
                          +drhodS(i,j+jq,k)*dTdy(index_salt)%field(i,j,k)) 
              slope    = -drhody/drhodz
              absslope = abs(slope)  

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Thickness%depth_zwt(i,j,k)) then                                           
                  depth_ratio     = Thickness%depth_zwt(i,j,k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! fill tensor components 
              tensor_22(i,j,jq,kr) = aredi_array(i,j,k)
              tensor_23(i,j,jq,kr) = aredi_array(i,j,k)*taper_slope

           enddo
        enddo
     enddo
  enddo

  ! tracer-dependent part of the calculation
  do n=1,num_prog_tracers
     sumz(:,:,:) = 0.0
     do kr=0,1
        do jq=0,1
           do j=jsc-1,jec
              do i=isc,iec
                 sumz(i,j,kr) = sumz(i,j,kr) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) &
                 *(tensor_22(i,j,jq,kr)*dTdy(n)%field(i,j,k) + tensor_23(i,j,jq,kr)*dTdz(n)%field(i,j+jq,k-1+kr)) &
                 * min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc-1,jec
        do i=isc,iec
           flux_y(n)%field(i,j,k) = Grd%dytnr(i,j)*(sumz(i,j,0)+sumz(i,j,1))
        enddo
     enddo
  enddo


  ! apply some masks 
  if(tmask_neutral_on) then
     do n=1,num_prog_tracers
        do j=jsc-1,jec
           do i=isc,iec
              if(k==Grd%kmt(i,j)) then
                 flux_y(n)%field(i,j,k) = aredi_array(i,j,k)*dTdy(n)%field(i,j,k)*Grd%dxtn(i,j) &
                               *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i,j+1,k,tau))
              endif
           enddo
        enddo
     enddo
     if(k==1) then
         do n=1,num_prog_tracers
            flux_y(n)%field(:,:,k) = &
            aredi_array(:,:,k)*dTdy(n)%field(:,:,k)*Grd%dxtn(:,:)*FMY(Thickness%rho_dzt(:,:,k,tau))
         enddo
     endif
  endif

  if(neutral_physics_limit) then 
      do n=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(T_prog(n)%tmask_limit(i,j,k)==1.0) then 
                   flux_y(n)%field(i,j,k) = aredi_array(i,j,k)*dTdy(n)%field(i,j,k)*Grd%dxtn(i,j) &
                                   *min(Thickness%rho_dzt(i,j,k,tau),Thickness%rho_dzt(i,j+1,k,tau))
               endif
            enddo
         enddo
      enddo
  endif

  call mpp_clock_end(id_clock_fy_flux_ndiffuse)

end subroutine fy_flux_ndiffuse
! </SUBROUTINE> NAME="fy_flux_ndiffuse"


!#######################################################################
! <SUBROUTINE NAME="fz_flux_ndiffuse">
!
! <DESCRIPTION>
! Subroutine computes the vertical neutral diffusion tracer flux component.
! Compute this component for all tracers at level k.
! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
!
! fz has physical dimensions (density*diffusivity*tracer gradient)
!
! This is nearly the same as the subroutine in ocean_nphysicsA. 
!
! </DESCRIPTION>
!
subroutine fz_flux_ndiffuse(T_prog, k)

  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)
  integer, intent(in) :: k

  integer :: i, j, n, ip, jq, kr
  real :: sumx_0, sumx_1, sumy_0, sumy_1
  real :: temparray31(isc:iec,jsc:jec,0:1,0:1)
  real :: temparray32(isc:iec,jsc:jec,0:1,0:1)

  if(tmask_neutral_on .and. k==1) then
     do n=1,num_prog_tracers
       fz2(n)%field(:,:) = 0.0 
      enddo
     return 
  endif

  if(k==nk) then 

     do n=1,num_prog_tracers
        do j=jsc,jec
           do i=isc,iec
              fz2(n)%field(i,j) = 0.0
           enddo
        enddo
     enddo
     return  

  elseif(k < nk) then  

      ! tracer-independent part of the calculation 
      do kr=0,1
         do ip=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray31(i,j,ip,kr) = tensor_31(i,j,k,ip,kr)*dtew(i,j,ip) &
                       *min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr))
               enddo
            enddo
         enddo
         do jq=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray32(i,j,jq,kr) = tensor_32(i,j,k,jq,kr)*dtns(i,j,jq) &
                       *min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr))  
               enddo
            enddo
         enddo
      enddo

      ! tracer-dependent part of the calculation  
      do n=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec
               sumx_0 =  temparray31(i,j,0,0)*dTdx(n)%field(i-1,j,k) &
                      +  temparray31(i,j,0,1)*dTdx(n)%field(i-1,j,k+1)
               sumx_1 =  temparray31(i,j,1,0)*dTdx(n)%field(i,j,k) &
                      +  temparray31(i,j,1,1)*dTdx(n)%field(i,j,k+1)
               sumy_0 =  temparray32(i,j,0,0)*dTdy(n)%field(i,j-1,k) &
                      +  temparray32(i,j,0,1)*dTdy(n)%field(i,j-1,k+1)
               sumy_1 =  temparray32(i,j,1,0)*dTdy(n)%field(i,j,k) &
                      +  temparray32(i,j,1,1)*dTdy(n)%field(i,j,k+1)

               ! compute time explicit portion of the vertical flux
               fz2(n)%field(i,j) =   Grd%tmask(i,j,k+1) &
                    *( Grd%dxtr(i,j)*(sumx_0+sumx_1) &
                      +Grd%dytr(i,j)*(sumy_0+sumy_1)) &
                    *dzwtr(i,j,k) &
                    + K33_explicit(i,j,k)*dTdz(n)%field(i,j,k)

            enddo
         enddo
      enddo
      
      if(tmask_neutral_on) then
          do n=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(k==(Grd%kmt(i,j)-1))  fz2(n)%field(i,j) = 0.0
                enddo
             enddo
          enddo
      endif

      if(neutral_physics_limit) then 
          do n=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(T_prog(n)%tmask_limit(i,j,k)==1.0) fz2(n)%field(i,j) = 0.0  
                enddo
             enddo
          enddo
      endif

  endif  !if-test for k-level 


end subroutine fz_flux_ndiffuse
! </SUBROUTINE> NAME="fz_flux_ndiffuse"


!#######################################################################
! <SUBROUTINE NAME="fx_flux_gm">
!
! <DESCRIPTION>
! Subroutine computes the zonal GM tracer flux component.
! Compute this component for all tracers at level k.
!
! fx has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fx_flux_gm

  integer :: i, j, k, n
  integer :: ip, kr
  real    :: tensor_13(isd:ied,jsd:jed,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)

  call mpp_clock_begin(id_clock_fx_flux_gm)

  do k=1,nk

     ! tracer-independent part of the calculation 
     tensor_13(:,:,:) = 0.0
     do ip=0,1 
        do j=jsc,jec
           do i=isc-1,iec
              tensor_13(i,j,ip) = -psiy(i,j,k,ip)
           enddo
        enddo
     enddo

     ! tracer-dependent part of the calculation
     do n=1,num_prog_tracers 
        sumz(:,:,:) = 0.0
        do kr=0,1
           do ip=0,1
              do j=jsc,jec
                 do i=isc-1,iec
                    sumz(i,j,kr) = sumz(i,j,kr) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) &
                         *tensor_13(i,j,ip)*dTdz(n)%field(i+ip,j,k-1+kr)                              &
                         *min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
                 enddo
              enddo
           enddo
        enddo
        do j=jsc,jec
           do i=isc-1,iec
              flux_x(n)%field(i,j,k) = Grd%dxter(i,j)*(sumz(i,j,0)+sumz(i,j,1))
           enddo
        enddo
     enddo

  enddo

  call mpp_clock_end(id_clock_fx_flux_gm)
  

end subroutine fx_flux_gm
! </SUBROUTINE> NAME="fx_flux_gm"



!#######################################################################
! <SUBROUTINE NAME="fy_flux_gm">
!
! <DESCRIPTION>
! Subroutine computes the meridional GM tracer flux component.
! Compute this component for all tracers at level k.
!
! fy has physical dimensions (area*diffusivity*density*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fy_flux_gm

  integer :: i, j, k, n
  integer :: jq, kr
  real    :: tensor_23(isd:ied,jsd:jed,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1)

  call mpp_clock_begin(id_clock_fy_flux_gm)

  do k=1,nk

     ! tracer-independent part of the calculation 
     tensor_23(:,:,:) = 0.0
     do jq=0,1  
        do j=jsc-1,jec
           do i=isc,iec
              tensor_23(i,j,jq) = psix(i,j,k,jq)
           enddo
        enddo
     enddo

     ! tracer-dependent part of the calculation
     do n=1,num_prog_tracers

        sumz(:,:,:) = 0.0
        do kr=0,1
           do jq=0,1
              do j=jsc-1,jec
                 do i=isc,iec
                    sumz(i,j,kr) = sumz(i,j,kr) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) &
                         *tensor_23(i,j,jq)*dTdz(n)%field(i,j+jq,k-1+kr)                              &
                         *min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
                 enddo
              enddo
           enddo
        enddo
        do j=jsc-1,jec
           do i=isc,iec
              flux_y(n)%field(i,j,k) = Grd%dytnr(i,j)*(sumz(i,j,0)+sumz(i,j,1))
           enddo
        enddo

     enddo

  enddo

  call mpp_clock_end(id_clock_fy_flux_gm)


end subroutine fy_flux_gm
! </SUBROUTINE> NAME="fy_flux_gm"


!#######################################################################
! <SUBROUTINE NAME="fz_flux_gm">
!
! <DESCRIPTION>
! Subroutine computes the vertical GM tracer flux component.
! Compute this component for all tracers at level k.
! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
!
! fz has physical dimensions (density*diffusivity*tracer gradient)
!
! </DESCRIPTION>
!
subroutine fz_flux_gm(T_prog, k)

  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)

  integer :: i, j, k, n, ip, jq, kr
  real :: sumx_0, sumx_1, sumy_0, sumy_1
  real :: temparray31(isc:iec,jsc:jec,0:1,0:1)
  real :: temparray32(isc:iec,jsc:jec,0:1,0:1)
  real :: tensor_31(isd:ied,jsd:jed,0:1)
  real :: tensor_32(isd:ied,jsd:jed,0:1)

  if(k==nk) then 

      do n=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec
               fz2(n)%field(i,j) = 0.0
            enddo
         enddo
      enddo
      return 

  elseif(k < nk) then  

      ! tracer independent portion of calculation 
      tensor_31   = 0.0
      tensor_32   = 0.0
      temparray31 = 0.0
      temparray32 = 0.0
      do ip=0,1
         jq=ip  
         do j=jsc-1,jec
            do i=isc-1,iec
               tensor_31(i,j,ip) =  psiy(i,j,k,ip)
               tensor_32(i,j,jq) = -psix(i,j,k,jq)
            enddo
         enddo
      enddo

      do kr=0,1

         do ip=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray31(i,j,ip,kr) = tensor_31(i,j,ip)*dtew(i,j,ip) &
                       *min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr))
               enddo
            enddo
         enddo
         do jq=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray32(i,j,jq,kr) = tensor_32(i,j,jq)*dtns(i,j,jq) &
                       *min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr))  
               enddo
            enddo
         enddo

      enddo

      ! tracer dependent portion 
      do n=1,num_prog_tracers 
         do j=jsc,jec
            do i=isc,iec
               sumx_0 =  temparray31(i,j,0,0)*dTdx(n)%field(i-1,j,k) &
                      +  temparray31(i,j,0,1)*dTdx(n)%field(i-1,j,k+1)
               sumx_1 =  temparray31(i,j,1,0)*dTdx(n)%field(i,j,k)   &
                      +  temparray31(i,j,1,1)*dTdx(n)%field(i,j,k+1)
               sumy_0 =  temparray32(i,j,0,0)*dTdy(n)%field(i,j-1,k) &
                      +  temparray32(i,j,0,1)*dTdy(n)%field(i,j-1,k+1)
               sumy_1 =  temparray32(i,j,1,0)*dTdy(n)%field(i,j,k)   &
                      +  temparray32(i,j,1,1)*dTdy(n)%field(i,j,k+1)

               fz2(n)%field(i,j) =  Grd%tmask(i,j,k+1)                                 &
                    *( Grd%dxtr(i,j)*(sumx_0+sumx_1) + Grd%dytr(i,j)*(sumy_0+sumy_1) ) &
                    *dzwtr(i,j,k)
            enddo
         enddo
      enddo

      if(neutral_physics_limit) then 
          do n=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(T_prog(n)%tmask_limit(i,j,k)==1.0) fz2(n)%field(i,j) = 0.0  
                enddo
             enddo
          enddo
      endif

  endif


end subroutine fz_flux_gm
! </SUBROUTINE> NAME="fz_flux_gm"


!#######################################################################
! <SUBROUTINE NAME="invtri_bvp">
!
! <DESCRIPTION>
! Solve the vertical diffusion equation implicitly using the
! method of inverting a tridiagonal matrix as described in
! Numerical Recipes in Fortran, The art of Scientific Computing,
! Second Edition, Press, Teukolsky, Vetterling, Flannery, 1992
! pages 42,43.
!
! enforce upsilon(k=kmt) = 0 via use of mask(k+1).  
!
! </DESCRIPTION>
!
subroutine invtri_bvp (mask, a, b, c, r, upsilon)

  real, intent(in), dimension(isd:,jsd:,:)      :: mask
  real, intent(in), dimension(isd:,jsd:,:)      :: a
  real, intent(in), dimension(isd:,jsd:,:)      :: b
  real, intent(in), dimension(isd:,jsd:,:)      :: c
  real, intent(in), dimension(isd:,jsd:,:,:)    :: r
  real, intent(inout), dimension(isd:,jsd:,:,:) :: upsilon

  real, dimension(isd:ied)    :: bet
  real, dimension(isd:ied,nk) :: gam

  integer :: i, j, k
  integer :: ipjq 

  bet = 0.0
  gam = 0.0
  
  do ipjq=1,2
     do j=jsc-1,jec


        ! decomposition and forward substitution
        k=1
        do i=isc-1,iec
           bet(i)              = mask(i,j,k+1)/(b(i,j,k) + epsln)
           upsilon(i,j,k,ipjq) = r(i,j,k,ipjq)*bet(i)
        enddo
        do k=2,nk-1
           do i=isc-1,iec
              gam(i,k)            = c(i,j,k-1)*bet(i)
              bet(i)              = mask(i,j,k+1)/(b(i,j,k) - a(i,j,k)*gam(i,k) + epsln)
              upsilon(i,j,k,ipjq) = (r(i,j,k,ipjq) - a(i,j,k)*upsilon(i,j,k-1,ipjq))*bet(i)
           enddo
        enddo

        ! back substitution
        do k=nk-1,1,-1
           do i=isc-1,iec
              upsilon(i,j,k,ipjq) = upsilon(i,j,k,ipjq) - gam(i,k+1)*upsilon(i,j,k+1,ipjq)
           enddo
        enddo


     enddo  ! enddo for j
  enddo     ! enddo for ipjq


end subroutine invtri_bvp
! </SUBROUTINE> NAME="invtri_bvp">




!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsC_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_nphysicsC_restart(time_stamp)
  character(len=*), intent(in), optional :: time_stamp

  if(.not. use_this_module) return

  call ocean_nphysics_util_restart(time_stamp)

end subroutine ocean_nphysicsC_restart
! </SUBROUTINE> NAME="ocean_nphysicsC_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsC_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_nphysicsC_end(Time)

  type(ocean_time_type), intent(in) :: Time
  
  if(.not. use_this_module) return

  call ocean_nphysics_coeff_end(Time, agm_array, aredi_array, rossby_radius, rossby_radius_raw, bczone_radius)

end subroutine ocean_nphysicsC_end
! </SUBROUTINE> NAME="ocean_nphysicsC_end"


end module ocean_nphysicsC_mod
