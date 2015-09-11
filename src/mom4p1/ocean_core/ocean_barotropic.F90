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
module ocean_barotropic_mod
!
! <CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies (mom4 algorithms)
! </CONTACT>
!
! <CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski (mom3 algorithm)
! </CONTACT>
!
! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang (OBC)
! </CONTACT>
!
! <CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> Martin Schmidt (OBC)
! </CONTACT>
!
! <CONTACT EMAIL="Harper.Simmons@noaa.gov"> Harper Simmons (tides)
! </CONTACT>
!
! <REVIEWER EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
! </REVIEWER>
!
!<OVERVIEW>
! Update the vertically integrated dynamics using a 
! split-explicit algorithm.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module time steps the vertically integrated dynamics. 
!
! Two explicit time stepping schemes are available:
!
! A. Leap-frog with optional Robert-Asselin time filter 
!
! B. Predictor-Corrector with adjustable dissipation
!
! Both use a split-explicit method.
!
! There is no rigid lid available in mom4. 
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, R.C. Pacanowski, R.M. Schmidt, and V. Balaji
! Tracer Conservation with an Explicit Free Surface Method for 
! Z-coordinate Ocean Models
! Monthly Weather Review (2001) vol 129 pages 1081--1098
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies
! Fundamentals of Ocean Climate Models
! Princeton University Press (2004)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4 (2004)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies: Elements of MOM4p1 (2006)
! </REFERENCE>
!
! </INFO>
!<NAMELIST NAME="ocean_barotropic_nml">
!
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="zero_tendency" TYPE="logical">
!  If true, will not integrate the barotropic fields.
!  </DATA> 
!  <DATA NAME="zero_forcing_bt" TYPE="logical">
!  Will set to zero all of the terms forcing the barotropic velocity.  
!  </DATA> 
!  <DATA NAME="zero_nonlinear_forcing_bt" TYPE="logical">
!  Will set to zero the nonlinear forcing terms, leaving only the smf and bmf 
!  terms to force the barotropic velocity. 
!  </DATA> 
!  <DATA NAME="zero_eta_ic" TYPE="logical">
!  To initialize eta_t to zero. 
!  </DATA> 
!  <DATA NAME="zero_eta_t" TYPE="logical">
!  To maintain eta_t at zero, but to allow other fields to evolve
!  For debugging.  Default zero_eta_t=.false.
!  </DATA> 
!  <DATA NAME="zero_eta_u" TYPE="logical">
!  To maintain eta_u at zero, but to allow other fields to evolve
!  For debugging.  Default zero_eta_u=.false.
!  </DATA> 
!  <DATA NAME="zero_eta_tendency" TYPE="logical">
!  To maintain deta_dt at zero.  For debugging. Default zero_eta_t=.false.
!  </DATA> 
!
!  <DATA NAME="ideal_initial_eta" TYPE="logical">
!  To initialize eta_t to an ideal profile.  This overrides 
!  all other initialization that may have occurred.  
!  Default=.false. 
!  </DATA> 
!  <DATA NAME="ideal_initial_eta_amplitude" TYPE="real" UNITS="metre">
!  Amplitude for initializing eta with an ideal profile.
!  Default ideal_initial_eta_amplitude = 5.0
!  </DATA> 
!  <DATA NAME="ideal_initial_eta_xwidth" TYPE="real" UNITS="metre">
!  Width in x-direction for sine-wave profile. 
!  Default xwidth=100e3
!  </DATA> 
!  <DATA NAME="ideal_initial_eta_ywidth" TYPE="real" UNITS="metre">
!  Width in y-direction for sine-wave profile. 
!  Default ywidth=100e3
!  </DATA> 
!
!  <DATA NAME="barotropic_time_stepping_mom4p0" TYPE="logical">
!  Use the general approach from mom4p0, in which the eta_t and 
!  pbot_t fields are updated with a big time step. This approach 
!  can be used with either barotropic predictor-corrector or 
!  barotropic leap-frog. 
!  Default barotropic_time_stepping_mom4p0=.false. 
!  </DATA> 
!
!  <DATA NAME="barotropic_time_stepping_mom4p1" TYPE="logical">
!  Use the alternative approach in which we assume the barotropic
!  scheme is a predictor-corrector.  In this way, the eta_t and 
!  pbot_t fields are updated with a time average. This approach 
!  is only available when using barotropic_pred_corr=.true. 
!  Default barotropic_time_stepping_mom4p1=.false. 
!  </DATA> 
!
!  <DATA NAME="barotropic_leap_frog" TYPE="logical">
!  Use leap-frog scheme for barotropic time stepping. 
!  Not the recommended method, since it requires smaller time 
!  steps.  It is maintained for legacy purposes.  
!  Default barotropic_leap_frog=.false.
!  </DATA> 
!  <DATA NAME="robert_asselin_bt" TYPE="real">
!  Robert time filter for use with leap-frog scheme for barotropic.
!  </DATA> 
!
!  <DATA NAME="barotropic_pred_corr" TYPE="logical">
!  Use preditor-corrector for barotropic time stepping. 
!  This is the recommended method.
!  Default barotropic_pred_corr=.true.
!  </DATA> 
!  <DATA NAME="pred_corr_gamma" UNITS="dimensionless" TYPE="real">
!  Dimensionless dissipation parameter for the preditor-corrector
!  scheme.  Setting pred_corr_gamma=0.0 reduces the scheme to a 
!  forward-backward, but it has been found to be unstable.  
!  So pred_corr_gamma > 0.0 is recommended.  Note that 
!  pred_corr_gamma > 0.25 may be over-dissipated and so may 
!  go unstable. 
!  </DATA> 
!
!  <DATA NAME="smooth_eta_t_bt_laplacian" TYPE="logical">
!  For spatially smoothing the eta_t field at each barotropic 
!  time step using a Laplacian operator.  May not be necessary when running 
!  with barotropic_pred_corr=.true. and pred_corr_gamma > 0.0, since 
!  predictor-corrector has dissipation from pred_corr_gamma > 0.0.
!  Applicable just for DEPTH_BASED vertical coordinates. 
!  </DATA> 
!  <DATA NAME="smooth_eta_t_bt_biharmonic" TYPE="logical">
!  For spatially smoothing the eta_t field at each barotropic 
!  time step using a biharmonic operator. May not be necessary when running
!  with barotropic_pred_corr=.true. and pred_corr_gamma > 0.0, since 
!  predictor-corrector has dissipation from pred_corr_gamma > 0.0. 
!  Applicable just for DEPTH_BASED vertical coordinates. 
!  WARNING: This operator is NOT positive definite, and so can 
!  produce spurious extrema.  It is not recommended just for this 
!  reason.
!  </DATA> 
!  <DATA NAME="smooth_eta_t_laplacian" TYPE="logical">
!  For spatially smoothing the eta_t field on the big time step
!  by using a laplacian operator. For compatibility 
!  and global conservation, must also introduce a mixing 
!  to the thickness weighted tracer concentration in the k=1 cell. 
!  Applicable just for DEPTH_BASED vertical coordinates. 
!  </DATA> 
!  <DATA NAME="smooth_eta_t_biharmonic" TYPE="logical">
!  For spatially smoothing the eta_t field on the big time step
!  by using a biharmonic operator. For compatibility 
!  and global conservation, must also introduce a mixing 
!  to the thickness weighted tracer concentration in the k=1 cell. 
!  Applicable just for DEPTH_BASED vertical coordinates. 
!  WARNING: This operator is NOT positive definite, and so can 
!  produce spurious extrema.  It is not recommended just for this 
!  reason.
!  </DATA> 
!  <DATA NAME="eta_offset" UNITS="metre" TYPE="real">
!  Uniform offset for use in determining the filter 
!  acting on tracer when smoothing the surface height. 
!  Default eta_offset=1e-12.
!  </DATA> 
!
!  <DATA NAME="smooth_eta_diag_laplacian" TYPE="logical">
!  For spatially smoothing the diagnosed eta_t field
!  using a laplacian operator. Default 
!  smooth_eta_diag_laplacian=.true.
!  </DATA> 
!  <DATA NAME="smooth_eta_diag_biharmonic" TYPE="logical">
!  For spatially smoothing the diagnosed eta_t field
!  using a biharmonic operator.  Default 
!  smooth_eta_diag_biharmonic=.false.
!  </DATA> 
!  <DATA NAME="vel_micom_lap_diag" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in the Laplacian smoothing of diagnosed surface height.
!  </DATA>
!  <DATA NAME="vel_micom_bih_diag" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM biharmonic mixing 
!  coefficient used in the biharmonic smoothing of diagnosed surface height.
!  </DATA>
!
!  <DATA NAME="smooth_anompb_bt_laplacian" TYPE="logical">
!  For spatially smoothing anomalous pbot_t at each barotropic 
!  time step using a Laplacian operator.  May not be necessary when running 
!  with barotropic_pred_corr=.true. and pred_corr_gamma > 0.0, since 
!  predictor-corrector has dissipation from pred_corr_gamma > 0.0.
!  Applicable just for PRESSURE_BASED vertical coordinates. 
!  </DATA> 
!  <DATA NAME="smooth_anompb_bt_biharmonic" TYPE="logical">
!  For spatially smoothing the anomalous pbot_t field at each barotropic 
!  time step using a biharmonic operator. May not be necessary when running
!  with barotropic_pred_corr=.true. and pred_corr_gamma > 0.0, since 
!  predictor-corrector has dissipation from pred_corr_gamma > 0.0. 
!  Applicable just for PRESSURE_BASED vertical coordinates. 
!  WARNING: This operator is NOT positive definite, and so can 
!  produce spurious extrema.  It is not recommended just for this 
!  reason.
!  </DATA> 
!  <DATA NAME="smooth_pbot_t_laplacian" TYPE="logical">
!  For spatially smoothing pbot_t-pbot0 on the big time step 
!  using a laplacian operator. For compatibility 
!  and global conservation, must also introduce a mixing 
!  to the thickness weighted tracer concentration in the k=kbot cell. 
!  Applicable just for PRESSURE_BASED vertical coordinates. 
!  </DATA> 
!  <DATA NAME="smooth_pbot_t_biharmonic" TYPE="logical">
!  For spatially smoothing pbot_t-pbot0 on the big time step 
!  by using a biharmonic operator. For compatibility 
!  and global conservation, must also introduce a mixing 
!  to the thickness weighted tracer concentration in the k=kbot cell. 
!  Applicable just for PRESSURE_BASED vertical coordinates. 
!  WARNING: This operator is NOT positive definite, and so can 
!  produce spurious extrema.  It is not recommended just for this 
!  reason.
!  </DATA> 
!  <DATA NAME="pbot_offset" UNITS="Pa" TYPE="real">
!  Uniform offset for use in determining the filter 
!  acting on tracer when smoothing the bottom pressure anomaly.
!  Default pbot_offset=1e-12.
!  </DATA> 
!
!  <DATA NAME="vel_micom_lap" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in the Laplacian smoothing of surface height
!  or anomalous bottom pressure.
!  </DATA>
!  <DATA NAME="vel_micom_bih" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM biharmonic mixing 
!  coefficient used in the biharmonic smoothing of surface height
!  or anomalous bottom pressure.
!  </DATA>
!
!  <DATA NAME="udrho_bt_lap" TYPE="logical">
!  For non-geopotential vertical coordinates, the vertically 
!  integrated horizontal momentum can be noisy. It is therefore
!  useful to add a smoothing operator.  Here, we apply the 
!  laplacian friction as coded in the friction module using 
!  the vertically averaged isotropic viscosity as well as a 
!  background. Do so on each barotropic time step. 
!  Default udrho_bt_lap=.false.
!  </DATA> 
!  <DATA NAME="udrho_bt_bih" TYPE="logical">
!  For non-geopotential vertical coordinates, the vertically 
!  integrated horizontal momentum can be noisy. It is therefore
!  useful to add a smoothing operator.  Here, we apply the 
!  biharmonic friction as coded in the friction module using 
!  the vertically averaged isotropic viscosity as well as a 
!  background. Do so on each barotropic time step.
!  Default udrho_bt_bih=.false.
!  </DATA> 
!  <DATA NAME="udrho_lap" TYPE="logical">
!  For non-geopotential vertical coordinates, the vertically 
!  integrated horizontal momentum can be noisy. It is therefore
!  useful to add a smoothing operator.  Here, we apply the 
!  laplacian friction as coded in the friction module using 
!  the vertically averaged isotropic viscosity as well as a 
!  background. Do so just on the baroclinic time step. 
!  Default udrho_lap=.false.
!  </DATA> 
!  <DATA NAME="udrho_bih" TYPE="logical">
!  For non-geopotential vertical coordinates, the vertically 
!  integrated horizontal momentum can be noisy. It is therefore
!  useful to add a smoothing operator.  Here, we apply the 
!  biharmonic friction as coded in the friction module using 
!  the vertically averaged isotropic viscosity as well as a 
!  background. Do so just on the baroclinic time step.
!  Default udrho_bih=.false.
!  </DATA> 
!
!  <DATA NAME="udrho_lap_vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in the Laplacian smoothing of udrho.
!  Default udrho_lap_vel_micom=.05
!  </DATA>
!  <DATA NAME="udrho_bih_vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM biharmonic mixing 
!  coefficient used in the biharmonic smoothing of udrho.
!  Default udrho_bih_vel_micom=.01 
!  </DATA>
!
!  <DATA NAME="tidal_forcing_m2" TYPE="logical">
!  Forces from lunar M2 tidal constituent. 
!  </DATA> 
!  <DATA NAME="tidal_forcing_8" TYPE="logical">
!  Forces from 8 lunar and solar tidal constituents. 
!  </DATA> 
!  <DATA NAME="tidal_forcing_ideal" TYPE="logical">
!  For ideal tidal forcing, which has a bump configuration.
!  </DATA> 
!
!  <DATA NAME="geoid_forcing" TYPE="logical">
!  For modifying the geoid, implemented as a time independent 
!  tidal forcing. Need to read in a file to obtain the offset
!  geoid profile.  
!  Default geoid_forcing=.false.
!  </DATA> 
!
!  <DATA NAME="barotropic_bmf" TYPE="logical">
!  To apply bottom drag over each barotropic time step. 
!  Note that enhanced drag from unresolved tides is not 
!  facilitated with this approach.  
!  Default barotropic_bmf=.false.
!  </DATA> 
!  <DATA NAME="barotropic_bmf_cdbot" TYPE="real" UNITS="dimensionless">
!  Dimensionless bottom drag coefficient for applying the bottom drag
!  on the barotropic time step. Default barotropic_bmf_cdbot=1e-3.
!  </DATA> 
!
!  <DATA NAME="truncate_eta" TYPE="logical">
!  To truncate the surface height deviation so to ensure positive thickness 
!  within the top cell. This method will not conserve volume or tracer. 
!  It is coded for cases when conservation is not critical but wish to 
!  run GEOPOTENTIAL models w/ large free surface height deviations, such
!  as when running with tides and very fine vertical resolution. 
!  </DATA> 
!  <DATA NAME="verbose_truncate" TYPE="logical">
!  For verbose printout on truncate_eta
!  </DATA> 
!  <DATA NAME="frac_crit_cell_height" UNITS="dimensionless" TYPE="real">
!  When use GEOPOTENTIAL vertical coordinate, the 
!  top model tracer grid cell has thickness dzt(i,j,1) = dzt(1) + eta_t(i,j).
!  0 < frac_crit_cell_height <= 1 sets the fraction of dzt(1) that is allowed
!  prior to bringing the model down due to overly small dzt(i,j,1).
!  </DATA> 
!  <DATA NAME="eta_max" UNITS="meter" TYPE="real">
!  The maximum positive eta_t allowed when truncate_eta is true. 
!  </DATA> 
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  Print out lots of diagnostics of use for debugging.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For brief or full printout on initialization
!  </DATA> 
!  <DATA NAME="diag_step" TYPE="integer">
!  Frequency for output of ascii barotropic diagnostics. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln, grav, radian, c2dbars, pi, seconds_per_day
use diag_manager_mod, only: register_diag_field, register_static_field, send_data, need_data
use fms_mod,          only: write_version_number, read_data, FATAL, WARNING, NOTE
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, file_exist
use fms_io_mod,       only: field_size
use fms_io_mod,       only: register_restart_field, save_restart, restore_state
use fms_io_mod,       only: reset_field_pointer, restart_file_type
use mpp_domains_mod,  only: BGRID_NE, NON_BITWISE_EXACT_SUM, BITWISE_EXACT_SUM
use mpp_domains_mod,  only: mpp_update_domains, mpp_global_field, mpp_global_sum, mpp_global_min
use mpp_domains_mod,  only: mpp_get_domain_components, mpp_get_layout, domain1d, mpp_get_pelist
use mpp_domains_mod,  only: NUPDATE, EUPDATE
use mpp_mod,          only: mpp_chksum, mpp_max, mpp_min, mpp_sum, mpp_root_pe, mpp_pe
use mpp_mod,          only: mpp_error, mpp_broadcast, stdout, stdlog
use mpp_mod,          only: mpp_send, mpp_recv, mpp_sync_self 
use time_manager_mod, only: time_type, increment_time, get_time

use ocean_bih_friction_mod, only: bih_friction_barotropic
use ocean_domains_mod,      only: get_local_indices, get_global_indices, set_ocean_domain
use ocean_domains_mod,      only: get_halo_sizes    
use ocean_lap_friction_mod, only: lap_friction_barotropic
use ocean_obc_mod,          only: ocean_obc_update_boundary, ocean_obc_barotropic, ocean_obc_damp_newton
use ocean_obc_mod,          only: ocean_obc_adjust_forcing_bt, ocean_obc_surface_height
use ocean_obc_mod,          only: ocean_obc_adjust_divud, ocean_obc_ud
use ocean_obc_mod,          only: ocean_obc_check_for_update
use ocean_operators_mod,    only: BDX_ET, BDY_NT, FMX, FMY, FAX, FAY, LAP_T
use ocean_operators_mod,    only: DIV_UD, REMAP_BT_TO_BU, GRAD_BAROTROPIC_P
use ocean_parameters_mod,   only: TWO_LEVEL, THREE_LEVEL
use ocean_parameters_mod,   only: GEOPOTENTIAL, ZSTAR, PRESSURE, PSTAR
use ocean_parameters_mod,   only: DEPTH_BASED, PRESSURE_BASED
use ocean_parameters_mod,   only: missing_value, rho0, rho0r, onehalf 
use ocean_types_mod,        only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,        only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,        only: ocean_external_mode_type, ocean_options_type
use ocean_types_mod,        only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,        only: ocean_adv_vel_type, ocean_prog_tracer_type
use ocean_util_mod,         only: write_timestamp
use ocean_workspace_mod,    only: wrk1_2d, wrk2_2d, wrk3_2d, wrk4_2d, wrk1_v2d, wrk2_v2d

implicit none

public ocean_barotropic_init
public eta_and_pbot_update
public eta_and_pbot_tendency
public eta_and_pbot_diagnose
public update_ocean_barotropic
public ocean_barotropic_forcing
public ocean_mass_forcing 
public ocean_eta_smooth
public ocean_pbot_smooth
public ocean_barotropic_end
public ocean_barotropic_restart

private leap_frog_barotropic_depth
private pred_corr_barotropic_depth
private pred_corr_barotropic_press
private tidal_forcing_init
private geoid_forcing_init
private get_tidal_forcing
private read_barotropic
private eta_smooth_diagnosed 

private eta_check 
private barotropic_energy
private barotropic_integrals
private psi_compute
private eta_terms_diagnose
private eta_truncate
private maximum_convrhoud
private barotropic_chksum
private ideal_initialize_eta 

private

logical  :: zero_tendency=.false.             ! if true, will freeze the barotropic fields 
logical  :: zero_eta_ic=.false.               ! if true, will initialize eta to zero.
logical  :: zero_eta_tendency=.false.         ! if true, will maintain deta_dt at zero.
logical  :: zero_eta_t=.false.                ! if true, will maintain eta_t at zero. 
logical  :: zero_eta_u=.false.                ! if true, will maintain eta_u at zero. 
logical  :: zero_forcing_bt=.false.           ! if true, will set Ext_mode%forcing_bt=0.0.  
logical  :: zero_nonlinear_forcing_bt=.false. ! if true, will only use smf and bmf for Ext_mode%forcing_bt.  

logical  :: ideal_initial_eta=.false.         ! for ideal initialization of eta 
real     :: ideal_initial_eta_amplitude=5.0   ! metre 
real     :: ideal_initial_eta_xwidth=100e3    ! metre
real     :: ideal_initial_eta_ywidth=100e3    ! metre

logical  :: truncate_eta = .false.    ! if true, will keep eta_t small enough to ensure dzt(1) + eta_t > 0.
logical  :: verbose_truncate = .true. ! for full printout of the truncate_eta points. 
real     :: frac_crit_cell_height=.20 ! fraction of dzt(1) that is deemed too thin for upper level dzt(i,j,1)
real     :: eta_max = 5.0             ! max amplitude surface height fluctuation (meter)

logical  :: debug_this_module=.false. ! for debugging--prints out a lot of checksums 
logical  :: verbose_init=.true.       ! for printouts during initialization  
integer  :: diag_step=-1              ! for printing ascii diagnostics 


! for vertical coordinate
integer  :: vert_coordinate
integer  :: vert_coordinate_class
 
! for tidal forcing 
logical  :: tidal_forcing_m2=.false.    ! for tidal forcing by just the M2 constituent 
logical  :: tidal_forcing_8=.false.     ! for tidal forcing by 8 lunar/solar constituents 
logical  :: tidal_forcing_ideal=.false. ! for ideal tidal forcing which has a bump configuration 
logical  :: tidal_forcing=.false.       ! internally set true if any tidal forcing enabled 

! for geoid forcing 
logical  :: geoid_forcing=.false.       ! for implementing a time independent tidal forcing,
                                        ! corresponding to a geoid modification


! to apply bottom drag on barotropic time step 
logical  :: barotropic_bmf       = .false.  
real     :: barotropic_bmf_cdbot = 1e-3

! for smoothing eta/pbot on each small barotropic time step or on the big-time step
logical  :: smooth_eta_t_bt_laplacian=.false.   ! to smooth eta_t_bt with Laplacian operator 
logical  :: smooth_eta_t_bt_biharmonic=.false.  ! to smooth eta_t_bt with biharmonic operator (not recommended)
logical  :: smooth_eta_t_laplacian=.true.       ! to smooth eta_t with Laplacian operator 
logical  :: smooth_eta_t_biharmonic=.false.     ! to smooth eta_t with biharmonic operator (not recommended)
logical  :: smooth_eta_diag_laplacian=.true.    ! to smooth diagnosed eta_t with Laplacian operator 
logical  :: smooth_eta_diag_biharmonic=.false.  ! to smooth diagnosed with biharmonic operator 
logical  :: smooth_anompb_bt_laplacian=.false.  ! to smooth pbot_t_bt-rho0*grav*ht with Laplacian operator 
logical  :: smooth_anompb_bt_biharmonic=.false. ! to smooth pbot_t_bt-rho0*grav*ht with biharmonic operator (not recommended)
logical  :: smooth_pbot_t_laplacian=.true.      ! to smooth pbot_t-pbot0 with Laplacian operator 
logical  :: smooth_pbot_t_biharmonic=.false.    ! to smooth pbot_t-pbot0 with biharmonic operator (not recommended)
real     :: vel_micom_lap=.05                   ! velocity (m/s) to set Laplacian mixing to smooth.
real     :: vel_micom_lap_diag=.2               ! velocity (m/s) to set Laplacian mixing to smooth diagnosed eta.
real     :: vel_micom_bih=.01                   ! velocity (m/s) to set biharmonic mixing to smooth.
real     :: vel_micom_bih_diag=.1               ! velocity (m/s) to set biharmonic mixing to smooth diagnosed eta.
real     :: eta_offset=1e-12                    ! surface height offset (m) for use in smoothing   
real     :: pbot_offset=1e-12                   ! bottom pressure offset (Pa) for use in smoothing 

! for spatially smoothing vertically integrated horizontal
! momentum on barotropic or baroclinic time steps using 
! friction. 
real     :: udrho_lap_vel_micom = .05   
real     :: udrho_bih_vel_micom = .01
logical  :: udrho_bt_lap        = .false. 
logical  :: udrho_bt_bih        = .false. 
logical  :: udrho_lap           = .false. 
logical  :: udrho_bih           = .false. 

! general options for time stepping eta_t and pbot_t 
logical  :: barotropic_time_stepping_mom4p0=.false.
logical  :: barotropic_time_stepping_mom4p1=.false.
logical  :: initsum_with_bar       =.true.
logical  :: initsum_with_bar_mom4p0=.false.
logical  :: initsum_with_bar_mom4p1=.true.

! barotropic leap-frog specific 
logical  :: barotropic_leap_frog=.false. ! if wish to use leap-frog time scheme
logical  :: splitting=.true.             ! internally set false if dtuv=dtbt
real     :: robert_asselin_bt=0.05       ! dimensionless parameter to damp leap-frog splitting mode if dtbt=dtuv  
real     :: twodt=0.0                    ! internally set to 2*dtbt
 
! barotropic predictor-corrector specific 
logical  :: barotropic_pred_corr=.true.  ! if wish to use predictor-corrector scheme 
real     :: pred_corr_gamma=0.2          ! dissipation parameter for predictor-corrector 

! for writing a rstart 
logical  :: write_a_restart=.true. 

! internally set 
logical  :: module_is_initialized=.false. ! to indicate that the module has been properly initialized 
logical  :: have_obc=.false.              ! for running with OBC
logical  :: update_domains_for_obc=.false.! if true call mpp_update_domains for OBC
integer  :: tendency        ! for discretization of time tendency ("threelevel" or "twolevel")
integer  :: nts=1           ! internally set to number of barotropic timesteps per baroclinic dtuv timestep
real     :: dtts            ! tracer time step (secs)
real     :: dtuv            ! baroclinic time step (secs)
real     :: dtbt            ! barotropic time step (secs)
real     :: dteta           ! timestep (secs) determining update of surface height or bottom pressure
real     :: dtime           ! timestep (secs) (dtime=2*dteta for threelevel and dtime=dteta for twolevel)
real     :: dtimer          ! 1/dtime
real     :: area_total_r    ! inverse area (1/m) of k=1 tracer cells 

real     :: p5gravr         ! 1/(2*grav)
real     :: p5grav_rho0r
real     :: grav_r          ! inverse gravitational acceleration 
real     :: grav_rho0r      ! 1.0/(grav*rho0)
real     :: rho0grav        ! rho0*grav
real     :: dtbt_gamma      ! dtbt * pred_corr_gamma 
real     :: dtbt_gamma_rho0r! dtbt * pred_corr_gamma * rho0r
real     :: dtbt_rho0r      ! dtbt * rho0r 

! for specifying transport units
! can either be Sv or mks
character(len=32) :: transport_dims ='Sv (10^9 kg/s)' 
real              :: transport_convert=1.0e-9 


! for psi_compute
integer, allocatable :: pelist_x(:)
integer, allocatable :: pelist_y(:)
integer  :: layout(2)
integer  :: myid_x
integer  :: myid_y
integer  :: size_x
integer  :: size_y
integer  :: this_pe
type(domain1d),save :: Domx
type(domain1d),save :: Domy


character(len=128) :: &
     version='$Id: ocean_barotropic.F90,v 1.1.4.1.2.62.6.12.6.2.2.6.4.3.22.1.4.3 2009/12/10 22:24:26 smg Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

#include <ocean_memory.h>

#ifdef MOM4_STATIC_ARRAYS

real, dimension(isd:ied,jsd:jed,3)   :: eta_t_bt   ! sea surface height(m) on T-cells on barotropic timesteps
real, dimension(isd:ied,jsd:jed)     :: eta_t_init ! sea surface height(m) on T-cells from ideal initial condition
real, dimension(isd:ied,jsd:jed,3)   :: anompb_bt  ! bottom pressure (Pa) on T-cells on barotropic timesteps
real, dimension(isd:ied,jsd:jed,2,3) :: udrho_bt   ! vertical average of depth*velocity on U-cells

real, dimension(isd:ied,jsd:jed)     :: thicku_r   ! inverse thickness of U-cell column (1/metre)

real, dimension(isd:ied,jsd:jed)     :: cori1      ! for the barotropic loop 
real, dimension(isd:ied,jsd:jed)     :: cori2      ! for the barotropic loop 

real, dimension(isd:ied,jsd:jed)     :: smooth_lap      ! coefficient (m^2/s) for lap smoothing of eta_t or pbot-pbot0
real, dimension(isd:ied,jsd:jed)     :: smooth_lap_diag ! coefficient (m^2/s) for lap smoothing of diagnosed eta_t 
real, dimension(isd:ied,jsd:jed)     :: smooth_bih      ! coefficient (m^4/s) for bih smoothing of eta_t or pbot-pbot0
real, dimension(isd:ied,jsd:jed)     :: smooth_bih_diag  ! coefficient (m^4/s) for bih smoothing of diagnosed eta_t
real, dimension(isd:ied,jsd:jed)     :: rho_g           ! rho_surf*grav
real, dimension(isd:ied,jsd:jed)     :: smooth_mask     ! mask array for use in smoothing
real, dimension(isd:ied,jsd:jed)     :: etastar         ! offset eta for use in smoothing
real, dimension(isd:ied,jsd:jed)     :: pbotstar        ! offset pbot for use in smoothing
real, dimension(isd:ied,jsd:jed)     :: tmp             ! temporary array for various calculations 
real, dimension(isd:ied,jsd:jed,2)   :: diffxy          ! difference array for smoothing

real, dimension(isd:ied,jsd:jed,2)   :: friction_lap ! friction operator applied to udrho 
real, dimension(isd:ied,jsd:jed,2)   :: friction_bih ! friction operator applied to udrho 

real, dimension(isd:ied,jsd:jed) :: lap_ceu_back     ! (m^2/sec) background viscosity for barotropic friction
real, dimension(isd:ied,jsd:jed) :: lap_cnu_back     ! (m^2/sec) background viscosity for barotropic friction
real, dimension(isd:ied,jsd:jed) :: bih_ceu_back     ! (m^4/sec) background viscosity for barotropic friction
real, dimension(isd:ied,jsd:jed) :: bih_cnu_back     ! (m^4/sec) background viscosity for barotropic friction
real, dimension(isd:ied,jsd:jed) :: smooth_lap_udrho ! coefficient (m^2/s) for lap smoothing of udrho
real, dimension(isd:ied,jsd:jed) :: smooth_bih_udrho ! coefficient (m^4/s) for bih smoothing of udrho

real, dimension(isd:ied,jsd:jed) :: eta_dynamic_tend   ! contribution to eta_nonbouss from dynamics (m)
real, dimension(isd:ied,jsd:jed) :: eta_water_tend     ! contribution to eta_nonbouss from water forcing (m)
real, dimension(isd:ied,jsd:jed) :: eta_steric_tend    ! contribution to eta_nonbouss from steric effect (m)
real, dimension(isd:ied,jsd:jed) :: eta_nonsteric_tend ! contribution to eta_nonbouss from nonsteric effect (m)
real, dimension(isd:ied,jsd:jed) :: eta_source_tend    ! contribution to eta_nonbouss from sources (m)
real, dimension(isd:ied,jsd:jed) :: eta_smooth_tend    ! contribution to eta_nonbouss from diagnosed eta smoothing (m)

real, dimension(isd:ied,jsd:jed) :: rhodzt          ! vertically integrated density at time tau (m*kg/m^3)
real, dimension(isd:ied,jsd:jed) :: rhodzt_inv      ! inverse of vertically integrated density (m^2/kg)
real, dimension(isd:ied,jsd:jed) :: rhodzt_taup1    ! vertically integrated density at time taup1 (m*kg/m^3)
#else

real, dimension(:,:,:), allocatable   :: eta_t_bt   ! sea surface height(m) on T-cells on barotropic timesteps
real, dimension(:,:), allocatable     :: eta_t_init ! sea surface height(m) on T-cells from ideal initial condition
real, dimension(:,:,:), allocatable   :: anompb_bt  ! bottom pressure (Pa) on T-cells on barotropic timesteps
real, dimension(:,:,:,:), allocatable :: udrho_bt   ! depth averaged velocity * ocean depth on U-cells

real, dimension(:,:), allocatable     :: thicku_r   ! inverse thickness of U-cell column (1/metre)

real, dimension(:,:), allocatable     :: cori1      ! for the barotropic loop
real, dimension(:,:), allocatable     :: cori2      ! for the barotropic loop

real, dimension(:,:), allocatable     :: smooth_lap      ! coefficient (m^2/s) for lap smoothing of eta_t or pbot-pbot0
real, dimension(:,:), allocatable     :: smooth_lap_diag ! coefficient (m^2/s) for lap smoothing of diagnosed eta_t
real, dimension(:,:), allocatable     :: smooth_bih      ! coefficient (m^4/s) for bih smoothing of eta_t or pbot-pbot0
real, dimension(:,:), allocatable     :: smooth_bih_diag ! coefficient (m^4/s) for bih smoothing of diagnosed eta_t
real, dimension(:,:), allocatable     :: rho_g           ! rho_surf*grav 
real, dimension(:,:), allocatable     :: smooth_mask     ! mask array for use in smoothing 
real, dimension(:,:), allocatable     :: etastar         ! offset eta for use in smoothing 
real, dimension(:,:), allocatable     :: pbotstar        ! offset pbot for use in smoothing 
real, dimension(:,:), allocatable     :: tmp             ! temporary array for various calculations 
real, dimension(:,:,:), allocatable   :: diffxy          ! temporary array for smoothing

real, dimension(:,:,:), allocatable   :: friction_lap ! friction operator applied to udrho 
real, dimension(:,:,:), allocatable   :: friction_bih ! friction operator applied to udrho 

real, dimension(:,:), allocatable  :: lap_ceu_back     ! (m^2/sec) background viscosity for barotropic friction
real, dimension(:,:), allocatable  :: lap_cnu_back     ! (m^2/sec) background viscosity for barotropic friction
real, dimension(:,:), allocatable  :: bih_ceu_back     ! (m^4/sec) background viscosity for barotropic friction
real, dimension(:,:), allocatable  :: bih_cnu_back     ! (m^4/sec) background viscosity for barotropic friction
real, dimension(:,:), allocatable  :: smooth_lap_udrho ! coefficient (m^2/s) for lap smoothing of udrho
real, dimension(:,:), allocatable  :: smooth_bih_udrho ! coefficient (m^4/s) for bih smoothing of udrho

real, dimension(:,:), allocatable  :: eta_dynamic_tend   ! contribution to eta_nonbouss from dynamics (m)
real, dimension(:,:), allocatable  :: eta_water_tend     ! contribution to eta_nonbouss from water forcing (m)
real, dimension(:,:), allocatable  :: eta_nonsteric_tend ! contribution to eta_nonbouss from nonsteric effect (m)
real, dimension(:,:), allocatable  :: eta_steric_tend    ! contribution to eta_nonbouss from steric effect (m)
real, dimension(:,:), allocatable  :: eta_source_tend    ! contribution to eta_nonbouss from sources (m)
real, dimension(:,:), allocatable  :: eta_smooth_tend    ! contribution to eta_nonbouss from diagnosed eta smoothing (m)
real, dimension(:,:), allocatable  :: rhodzt            ! vertically integrated density at time tau (m*kg/m^3)
real, dimension(:,:), allocatable  :: rhodzt_inv        ! inverse of vertically integrated density (m^2/kg)
real, dimension(:,:), allocatable  :: rhodzt_taup1      ! vertically integrated density at time taup1 (m*kg/m^3)

#endif

! for diagnostics manager 
logical :: used

integer :: id_eta_t_mod          =-1
integer :: id_eta_t              =-1
integer :: id_eta_t_sq           =-1
integer :: id_eta_t_bt           =-1
integer :: id_deta_dt            =-1
integer :: id_eta_u              =-1
integer :: id_eta_t_bar          =-1
integer :: id_smooth_lap         =-1
integer :: id_smooth_lap_diag    =-1
integer :: id_smooth_bih         =-1
integer :: id_smooth_bih_diag    =-1
integer :: id_patm_for_sea_lev   =-1
integer :: id_sea_lev_for_coupler=-1


integer :: id_smooth_lap_udrho =-1
integer :: id_smooth_bih_udrho =-1
integer :: id_udrho_bt_lap     =-1
integer :: id_udrho_bt_bih     =-1
integer :: id_vdrho_bt_lap     =-1
integer :: id_vdrho_bt_bih     =-1
integer :: id_udrho_lap        =-1
integer :: id_udrho_bih        =-1
integer :: id_vdrho_lap        =-1
integer :: id_vdrho_bih        =-1

integer :: id_pbot_t     =-1
integer :: id_anompb     =-1
integer :: id_anompb_bt  =-1
integer :: id_pbot_u     =-1
integer :: id_dpbot_dt   =-1

integer :: id_forcing_u_bt         =-1
integer :: id_forcing_v_bt         =-1
integer :: id_nonlin_forcing_u_bt  =-1
integer :: id_nonlin_forcing_v_bt  =-1

integer :: id_rhobarz = -1
integer :: id_rhoavg  = -1

integer :: id_sea_level          =-1
integer :: id_sea_level_sq       =-1
integer :: id_eta_nonbouss       = -1
integer :: id_eta_tend           = -1
integer :: id_eta_dynamic_tend   = -1
integer :: id_eta_water_tend     = -1
integer :: id_eta_steric_tend    = -1
integer :: id_eta_steric_tend_A  = -1
integer :: id_eta_steric_tend_B  = -1
integer :: id_eta_nonsteric_tend = -1
integer :: id_eta_source_tend    = -1
integer :: id_eta_smooth_tend    = -1
integer :: id_eta_nonsteric      = -1
integer :: id_eta_steric         = -1
integer :: id_eta_water          = -1
integer :: id_eta_dynamic        = -1
integer :: id_eta_source         = -1
integer :: id_eta_smooth         = -1
integer :: id_eta_surf_temp      = -1
integer :: id_eta_surf_salt      = -1
integer :: id_eta_surf_water     = -1
integer :: id_eta_bott_temp      = -1
integer :: id_eta_stemp_forcing  = -1
integer :: id_eta_btemp_forcing  = -1
integer :: id_eta_salt_forcing   = -1
integer :: id_eta_water_forcing  = -1
integer :: id_eta_dpress_dt      = -1

integer :: id_eta_global            = -1
integer :: id_eta_nonbouss_global   = -1
integer :: id_eta_steric_global     = -1
integer :: id_eta_nonsteric_global  = -1
integer :: id_eta_dynamic_global    = -1
integer :: id_eta_water_global      = -1
integer :: id_eta_source_global     = -1
integer :: id_eta_smooth_global     = -1
integer :: id_eta_surf_temp_global  = -1
integer :: id_eta_surf_salt_global  = -1
integer :: id_eta_surf_water_global = -1
integer :: id_eta_bott_temp_global  = -1

integer :: id_eta_tend_global           = -1
integer :: id_eta_dynamic_tend_global   = -1
integer :: id_eta_water_tend_global     = -1
integer :: id_eta_steric_tend_global    = -1
integer :: id_eta_nonsteric_tend_global = -1
integer :: id_eta_source_tend_global    = -1
integer :: id_eta_smooth_tend_global    = -1

integer :: id_eta_water_forcing_global = -1
integer :: id_eta_stemp_forcing_global = -1
integer :: id_eta_salt_forcing_global  = -1
integer :: id_eta_btemp_forcing_global = -1

integer :: id_patm_t   =-1
integer :: id_patm_u   =-1
integer :: id_dpatm_dt =-1

integer :: id_urhod          =-1
integer :: id_vrhod          =-1
integer :: id_psiu           =-1
integer :: id_psiv           =-1

integer :: id_ps             =-1
integer :: id_psx            =-1
integer :: id_psy            =-1

integer :: id_pb             =-1
integer :: id_grad_anompbx   =-1
integer :: id_grad_anompby   =-1

integer :: id_eta_smoother   =-1
integer :: id_pbot_smoother  =-1
integer :: id_smooth_mask    =-1

integer :: id_conv_rho_ud_t  =-1
integer :: id_eta_eq_tidal   =-1
integer :: id_ke_bt          =-1
integer :: id_pe_bt          =-1
integer :: id_eta_geoid      =-1

integer :: id_pme_velocity    =-1
integer :: id_pme_total       =-1
integer :: id_river_total     =-1
integer :: id_etat_avg        =-1
integer :: id_mass_total      =-1
integer :: id_ext_mode_source =-1

integer :: xhalo, yhalo, halo

! for restart
integer                       :: id_restart(22) = 0
type(restart_file_type), save :: Bar_restart

! for ascii output
integer :: unit=6

! for converting pme to m/sec 
real :: rho_water_r=0.001

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), save    :: Dom_flux

! for geoid forcing  
real, private, dimension(:,:), allocatable :: eta_geoid ! modified geoid (m)

! for tidal forcing 
real :: alphat=1.0  ! ocean loading and self-attraction from tidal forcing (=0.948)
real :: rradian     ! inverse radians
real :: tidal_omega_M2,Love_M2,amp_M2
real :: tidal_omega_S2,Love_S2,amp_S2
real :: tidal_omega_N2,Love_N2,amp_N2
real :: tidal_omega_K2,Love_K2,amp_K2
real :: tidal_omega_K1,Love_K1,amp_K1
real :: tidal_omega_O1,Love_O1,amp_O1
real :: tidal_omega_P1,Love_P1,amp_P1
real :: tidal_omega_Q1,Love_Q1,amp_Q1

! for tidal and geoid forcing  
real, private, dimension(:,:), allocatable :: eta_eq_tidal! equilibrium tidal forcing (m)
real, private, dimension(:,:), allocatable :: cos2lon     ! cosine   of 2*longitude on T-cells
real, private, dimension(:,:), allocatable :: coslon      ! cosine   of   longitude on T-cells
real, private, dimension(:,:), allocatable :: coslat2     ! cosine^2 of   longitude on T-cells
real, private, dimension(:,:), allocatable :: sin2lon     !   sine   of 2*longitude on T-cells
real, private, dimension(:,:), allocatable :: sinlon      !   sine   of   longitude on T-cells
real, private, dimension(:,:), allocatable :: sin2lat     !   sine   of 2*latitude  on T-cells
real, private, dimension(:,:), allocatable :: ideal_amp   ! amplitude for ideal tidal forcing 


namelist /ocean_barotropic_nml/ write_a_restart,                                &
         zero_tendency, zero_eta_ic, zero_eta_t, zero_eta_u,                    &
         zero_nonlinear_forcing_bt, zero_forcing_bt, zero_eta_tendency,         &
         tidal_forcing_m2, tidal_forcing_8, tidal_forcing_ideal, geoid_forcing, &
         barotropic_pred_corr, pred_corr_gamma,                                 &
         barotropic_leap_frog, robert_asselin_bt,                               &
         smooth_eta_t_bt_laplacian, smooth_eta_t_bt_biharmonic,                 &
         smooth_eta_t_laplacian, smooth_eta_t_biharmonic,                       &
         smooth_anompb_bt_laplacian, smooth_anompb_bt_biharmonic,               &
         smooth_pbot_t_laplacian, smooth_pbot_t_biharmonic,                     &
         smooth_eta_diag_laplacian, smooth_eta_diag_biharmonic,                 &
         vel_micom_lap, vel_micom_lap_diag, vel_micom_bih, vel_micom_bih_diag,  &   
         truncate_eta, verbose_truncate, eta_max, frac_crit_cell_height,        &
         verbose_init, debug_this_module, diag_step,                            &
         eta_offset, pbot_offset,                                               &
         barotropic_time_stepping_mom4p0, barotropic_time_stepping_mom4p1,      &
	 initsum_with_bar_mom4p0, initsum_with_bar_mom4p1,                      &
         ideal_initial_eta, ideal_initial_eta_amplitude,                        &
         ideal_initial_eta_xwidth, ideal_initial_eta_ywidth,                    &
         barotropic_bmf, barotropic_bmf_cdbot,                                  &
         udrho_bt_lap, udrho_bt_bih, udrho_lap, udrho_bih,                      &
         udrho_lap_vel_micom, udrho_bih_vel_micom
                                
contains


!#######################################################################
! <SUBROUTINE NAME="ocean_barotropic_init">
!
! <DESCRIPTION>
! Initialize the barotropic module
! </DESCRIPTION>
!
subroutine ocean_barotropic_init(Grid, Domain, Time, Time_steps, Ocean_options, Ext_mode, obc, &
                                 ver_coordinate, ver_coordinate_class, cmip_units, debug)

  type(ocean_grid_type),          intent(in), target   :: Grid
  type(ocean_domain_type),        intent(in), target   :: Domain
  type(ocean_time_type),          intent(in)           :: Time
  type(ocean_time_steps_type),    intent(inout)        :: Time_steps 
  type(ocean_options_type),       intent(inout)        :: Ocean_options
  type(ocean_external_mode_type), intent(inout)        :: Ext_mode
  logical,                        intent(in)           :: obc
  integer,                        intent(in)           :: ver_coordinate
  integer,                        intent(in)           :: ver_coordinate_class
  logical,                        intent(in)           :: cmip_units
  logical,                        intent(in), optional :: debug

  character(len=40) model_type 
  logical :: cfl_error=.false.
  real    :: cfl_grid_factor
  real    :: crit, fudge
  real    :: cgmax, gridmin, cgext, gridsp, dtcg
  real    :: max_dt_for_cgext, max_dt_for_cgext0
  integer :: tau, taup1
  integer :: ioun, io_status, ierr
  integer :: i, j, n
  integer :: icg, jcg
  integer :: smooth_methods=0

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_barotropic_mod (ocean_barotropic_init): module already initialized.')
  endif 

  module_is_initialized = .TRUE.

  if (diag_step == 0) diag_step = 1

  call write_version_number( version, tagname )

  ioun = open_namelist_file()
  read  (ioun, ocean_barotropic_nml,iostat=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_barotropic_nml)  
  write (stdlogunit, ocean_barotropic_nml)
  ierr = check_nml_error(io_status, 'ocean_barotropic_nml')
  call close_file(ioun)

  if (PRESENT(debug) .and. .not. debug_this_module) then
    debug_this_module = debug
  endif 
  if(debug_this_module) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_barotropic with debug_this_module=.true.'  
  endif 

  if(.not. write_a_restart) then 
    write(stdoutunit,'(a)') '==>Note: running ocean_barotropic with write_a_restart=.false.'
    write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
  endif 

  if(pred_corr_gamma > 0.25 .and. barotropic_pred_corr) then 
    write(stdoutunit,'(a)') '==>Warning: ocean_barotropic_mod: pred_corr_gamma > 0.25 may be unstable.'
  endif 

  have_obc = obc
  if (have_obc) update_domains_for_obc = ocean_obc_check_for_update()
  vert_coordinate       = ver_coordinate 
  vert_coordinate_class = ver_coordinate_class
  
  Dom => Domain
  Grd => Grid
  call set_ocean_domain(Dom_flux, Grid, name='horz diff flux', maskmap=Dom%maskmap)

#ifndef MOM4_STATIC_ARRAYS

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Domain, isg, ieg, jsg, jeg)
  nk = Grid%nk
  
  call get_halo_sizes(Domain, xhalo, yhalo)
  if (xhalo /= yhalo .and. have_obc) then 
     call mpp_error(FATAL,&
     '==>Error in ocean_barotropic_init: with static memory and OBC, then xhalo must equal yhalo')
  endif 
  halo = xhalo
  
  allocate (rho_g(isd:ied,jsd:jed))

  allocate (Ext_mode%eta_t(isd:ied,jsd:jed,3))
  allocate (Ext_mode%eta_u(isd:ied,jsd:jed,3))
  allocate (Ext_mode%eta_t_bar(isd:ied,jsd:jed,3))
  allocate (Ext_mode%deta_dt(isd:ied,jsd:jed)) 

  allocate (Ext_mode%pbot_t(isd:ied,jsd:jed,3))
  allocate (Ext_mode%pbot_u(isd:ied,jsd:jed,3))
  allocate (Ext_mode%anompb(isd:ied,jsd:jed,3))
  allocate (Ext_mode%anompb_bar(isd:ied,jsd:jed,3))
  allocate (Ext_mode%dpbot_dt(isd:ied,jsd:jed))

  allocate (Ext_mode%dpatm_dt(isd:ied,jsd:jed))
  allocate (Ext_mode%patm_t(isd:ied,jsd:jed,3))
  allocate (Ext_mode%patm_for_sea_lev(isd:ied,jsd:jed))
  allocate (Ext_mode%patm_u(isd:ied,jsd:jed))

  allocate (Ext_mode%conv_rho_ud_t(isd:ied,jsd:jed,3))
  allocate (Ext_mode%udrho(isd:ied,jsd:jed,2,3))
  allocate (Ext_mode%forcing_bt(isd:ied,jsd:jed,2))
  allocate (Ext_mode%ps(isd:ied,jsd:jed))
  allocate (Ext_mode%grad_ps(isd:ied,jsd:jed,2))
  allocate (Ext_mode%grad_anompb(isd:ied,jsd:jed,2))
  allocate (Ext_mode%press_force(isd:ied,jsd:jed,2))

  allocate (Ext_mode%source(isd:ied,jsd:jed))
  allocate (Ext_mode%eta_smooth(isd:ied,jsd:jed))
  allocate (Ext_mode%pbot_smooth(isd:ied,jsd:jed))

  allocate (Ext_mode%eta_nonbouss(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_nonsteric(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_steric(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_water(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_source(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_dynamic(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_surf_temp(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_surf_salt(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_surf_water(isd:ied,jsd:jed,3)) 
  allocate (Ext_mode%eta_bott_temp(isd:ied,jsd:jed,3)) 

  allocate (eta_t_bt(isd:ied,jsd:jed,3))
  allocate (eta_t_init(isd:ied,jsd:jed))
  allocate (anompb_bt(isd:ied,jsd:jed,3))
  allocate (udrho_bt(isd:ied,jsd:jed,2,3))

  allocate (thicku_r(isd:ied,jsd:jed))
  allocate (cori1(isd:ied,jsd:jed))
  allocate (cori2(isd:ied,jsd:jed))
  
  allocate (smooth_mask(isd:ied,jsd:jed))
  allocate (etastar(isd:ied,jsd:jed))
  allocate (pbotstar(isd:ied,jsd:jed))
  allocate (tmp(isd:ied,jsd:jed))
  allocate (diffxy(isd:ied,jsd:jed,2))

  allocate (eta_dynamic_tend(isd:ied,jsd:jed))
  allocate (eta_water_tend(isd:ied,jsd:jed))
  allocate (eta_steric_tend(isd:ied,jsd:jed))
  allocate (eta_nonsteric_tend(isd:ied,jsd:jed))
  allocate (eta_source_tend(isd:ied,jsd:jed))
  allocate (eta_smooth_tend(isd:ied,jsd:jed))
  allocate (rhodzt(isd:ied,jsd:jed))
  allocate (rhodzt_inv(isd:ied,jsd:jed))
  allocate (rhodzt_taup1(isd:ied,jsd:jed))
  
#endif
  
  tendency = Time_steps%tendency
  dtts     = Time_steps%dtts
  dtuv     = Time_steps%dtuv
  dtbt     = Time_steps%dtbt
  dteta    = Time_steps%dteta
  dtime    = Time_steps%dtime_e
  dtimer   = 1.0/(dtime+epsln)
  
  Ext_mode%eta_t(:,:,:)       = 0.0
  Ext_mode%eta_u(:,:,:)       = 0.0
  Ext_mode%eta_t_bar(:,:,:)   = 0.0
  Ext_mode%deta_dt(:,:)       = 0.0

  Ext_mode%dpbot_dt(:,:)        = 0.0
  Ext_mode%dpatm_dt(:,:)        = 0.0
  Ext_mode%patm_t(:,:,:)        = 0.0
  Ext_mode%patm_for_sea_lev(:,:)= 0.0
  Ext_mode%patm_u(:,:)          = 0.0

  Ext_mode%conv_rho_ud_t(:,:,:)= 0.0
  Ext_mode%udrho(:,:,:,:)      = 0.0
  Ext_mode%forcing_bt(:,:,:)   = 0.0
  Ext_mode%ps(:,:)             = 0.0
  Ext_mode%grad_ps(:,:,:)      = 0.0
  Ext_mode%press_force(:,:,:)  = 0.0
  Ext_mode%grad_anompb(:,:,:)  = 0.0
  Ext_mode%source(:,:)         = 0.0
  Ext_mode%eta_smooth(:,:)     = 0.0
  Ext_mode%pbot_smooth(:,:)    = 0.0

  Ext_mode%eta_nonbouss(:,:,:)   = 0.0
  Ext_mode%eta_nonsteric(:,:,:)  = 0.0
  Ext_mode%eta_steric(:,:,:)     = 0.0
  Ext_mode%eta_water(:,:,:)      = 0.0
  Ext_mode%eta_source(:,:,:)     = 0.0
  Ext_mode%eta_dynamic(:,:,:)    = 0.0
  Ext_mode%eta_surf_temp(:,:,:)  = 0.0
  Ext_mode%eta_surf_salt(:,:,:)  = 0.0
  Ext_mode%eta_surf_water(:,:,:) = 0.0
  Ext_mode%eta_bott_temp(:,:,:)  = 0.0

  eta_t_bt(:,:,:)   = 0.0
  eta_t_init(:,:)   = 0.0
  udrho_bt(:,:,:,:) = 0.0
  thicku_r(:,:)     = 0.0

  ! initialize the bottom pressure.
  ! note that for PRESSURE_BASED models, these 
  ! values are recomputed at time=0 
  ! based on Thickness%pbot0. 
  do n=1,3
    Ext_mode%pbot_t(:,:,n)     = rho0*grav*Grd%ht(:,:)
    Ext_mode%pbot_u(:,:,n)     = rho0*grav*Grd%hu(:,:)
    Ext_mode%anompb(:,:,n)     = 0.0
    Ext_mode%anompb_bar(:,:,n) = 0.0
    anompb_bt(:,:,n)           = 0.0
  enddo

  do j=jsd,jed
    do i=isd,ied
      cori1(i,j) = 0.5*dtbt*Grd%f(i,j) 
      cori2(i,j) = 1.0/(1.0 + cori1(i,j)**2)
    enddo
  enddo

  tmp(:,:)          = 0.0
  rho_g(:,:)        = grav*rho0

  smooth_mask(:,:)  = 0.0
  etastar(:,:)      = 0.0
  pbotstar(:,:)     = 0.0
  diffxy(:,:,:)     = 0.0

  eta_dynamic_tend(:,:)   = 0.0
  eta_water_tend(:,:)     = 0.0
  eta_steric_tend(:,:)    = 0.0
  eta_nonsteric_tend(:,:) = 0.0
  eta_source_tend(:,:)    = 0.0
  eta_smooth_tend(:,:)    = 0.0
  rhodzt(:,:)             = 0.0
  rhodzt_inv(:,:)         = 0.0
  rhodzt_taup1(:,:)       = 0.0

  if(have_obc) then
    area_total_r      = mpp_global_sum(Dom%domain2d,Grd%dat(:,:)*Grd%tmask(:,:,1)* Grd%obc_tmask(:,:), BITWISE_EXACT_SUM)
  else
    area_total_r      = mpp_global_sum(Dom%domain2d,Grd%dat(:,:)*Grd%tmask(:,:,1), BITWISE_EXACT_SUM)
  endif
  area_total_r      = 1.0/(epsln+area_total_r) 
  grav_r            = 1.0/grav
  rho0grav          = rho0*grav
  grav_rho0r        = rho0r/grav
  p5gravr           = 0.5*grav_r
  p5grav_rho0r      = 0.5*grav*rho0r 
  dtbt_gamma        = dtbt*pred_corr_gamma
  dtbt_gamma_rho0r  = dtbt*pred_corr_gamma*rho0r
  dtbt_rho0r        = dtbt*rho0r 


#ifndef MOM4_STATIC_ARRAYS    
  allocate (smooth_lap(isd:ied,jsd:jed))
  allocate (smooth_lap_diag(isd:ied,jsd:jed))
  allocate (smooth_bih(isd:ied,jsd:jed))
  allocate (smooth_bih_diag(isd:ied,jsd:jed))
  allocate (smooth_lap_udrho(isd:ied,jsd:jed))
  allocate (smooth_bih_udrho(isd:ied,jsd:jed))
  allocate (friction_lap(isd:ied,jsd:jed,2))
  allocate (friction_bih(isd:ied,jsd:jed,2))
  allocate (lap_ceu_back(isd:ied,jsd:jed))
  allocate (lap_cnu_back(isd:ied,jsd:jed))
  allocate (bih_ceu_back(isd:ied,jsd:jed))
  allocate (bih_cnu_back(isd:ied,jsd:jed))
#endif    
  friction_lap(:,:,:)= 0.0
  friction_bih(:,:,:)= 0.0
  lap_ceu_back(:,:)  = 0.0
  lap_cnu_back(:,:)  = 0.0
  bih_ceu_back(:,:)  = 0.0
  bih_cnu_back(:,:)  = 0.0

  smooth_lap(:,:) = &
   Grd%tmask(:,:,1)*vel_micom_lap*(2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:)))
  smooth_lap_diag(:,:) = &
   Grd%tmask(:,:,1)*vel_micom_lap_diag*(2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:)))
  smooth_bih(:,:) = &
   Grd%tmask(:,:,1)*vel_micom_bih*(2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:)))**3
  smooth_bih_diag(:,:) = &
   Grd%tmask(:,:,1)*vel_micom_bih_diag*(2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:)))**3

  smooth_lap_udrho(:,:) = &
   Grd%umask(:,:,1)*udrho_lap_vel_micom*(2.0*Grd%dxu(:,:)*Grd%dyu(:,:)/(Grd%dxu(:,:)+Grd%dyu(:,:)))
  smooth_bih_udrho(:,:) = &
   Grd%umask(:,:,1)*udrho_bih_vel_micom*(2.0*Grd%dxu(:,:)*Grd%dyu(:,:)/(Grd%dxu(:,:)+Grd%dyu(:,:)))**3

  ! compute background viscosity on east and north face of U-cells 
  do j=jsc,jec
     do i=isc,iec
        lap_ceu_back(i,j) = 0.5*(smooth_lap_udrho(i,j)+smooth_lap_udrho(i+1,j))
        lap_cnu_back(i,j) = 0.5*(smooth_lap_udrho(i,j)+smooth_lap_udrho(i,j+1))
        bih_ceu_back(i,j) = 0.5*(smooth_bih_udrho(i,j)+smooth_bih_udrho(i+1,j))
        bih_cnu_back(i,j) = 0.5*(smooth_bih_udrho(i,j)+smooth_bih_udrho(i,j+1))
     enddo
  enddo

  call mpp_update_domains (lap_ceu_back, Dom%domain2d)
  call mpp_update_domains (lap_cnu_back, Dom%domain2d)
  call mpp_update_domains (bih_ceu_back, Dom%domain2d)
  call mpp_update_domains (bih_cnu_back, Dom%domain2d)
  bih_ceu_back(:,:) = sqrt(bih_ceu_back(:,:))
  bih_cnu_back(:,:) = sqrt(bih_cnu_back(:,:))

  ! for units placed on diagnosed transport 
  if(cmip_units) then
      transport_convert=1.0
      transport_dims   = 'kg/s'
  else
      transport_convert=1.0e-9 
      transport_dims   = 'Sv (10^9 kg/s)'
  endif


  ! static fields for diagnostic manager 
  id_smooth_lap = register_static_field ('ocean_model', 'smooth_lap', Grd%tracer_axes(1:2), &
                   'laplacian mixing coefficient for smoothing', 'm^2/s',                   &
                   missing_value=missing_value, range=(/0.0,1e12/))
  if (id_smooth_lap > 0) used = send_data(id_smooth_lap, smooth_lap(isc:iec,jsc:jec), Time%model_time)

  id_smooth_lap_diag = register_static_field ('ocean_model', 'smooth_lap_diag', Grd%tracer_axes(1:2), &
                   'laplacian mixing coefficient for smoothing diagnosed eta', 'm^2/s',               &
                   missing_value=missing_value, range=(/0.0,1e12/))
  if (id_smooth_lap_diag > 0) used = send_data(id_smooth_lap_diag, smooth_lap_diag(isc:iec,jsc:jec), Time%model_time)

  id_smooth_bih = register_static_field ('ocean_model', 'smooth_bih', Grd%tracer_axes(1:2), &
                   'biharmonic mixing coefficient for smoothing', 'm^4/s',                  &
                   missing_value=missing_value, range=(/0.0,1e12/))
  if (id_smooth_bih > 0) used = send_data(id_smooth_bih, smooth_bih(isc:iec,jsc:jec), Time%model_time)

  id_smooth_bih_diag = register_static_field ('ocean_model', 'smooth_bih_diag', Grd%tracer_axes(1:2), &
                   'biharmonic mixing coefficient for smoothing diagnosed eta', 'm^4/s',              &
                   missing_value=missing_value, range=(/0.0,1e12/))
  if (id_smooth_bih_diag > 0) used = send_data(id_smooth_bih_diag, smooth_bih_diag(isc:iec,jsc:jec), Time%model_time)

  id_smooth_lap_udrho = register_static_field ('ocean_model', 'smooth_lap_udrho',            &
                   Grd%vel_axes_uv(1:2), 'laplacian mixing coefficient for smoothing udrho', &
                   'm^2/s', missing_value=missing_value, range=(/0.0,1e12/))
  if (id_smooth_lap_udrho > 0) used = send_data(id_smooth_lap_udrho, &
                                      smooth_lap_udrho(isc:iec,jsc:jec), Time%model_time)
  id_smooth_bih_udrho = register_static_field ('ocean_model', 'smooth_bih_udrho',             &
                   Grd%vel_axes_uv(1:2), 'biharmonic mixing coefficient for smoothing udrho', &
                   'm^4/s', missing_value=missing_value, range=(/0.0,1e20/))
  if (id_smooth_bih_udrho > 0) used = send_data(id_smooth_bih_udrho, &
                                      smooth_bih_udrho(isc:iec,jsc:jec), Time%model_time)

  ! convert smooth_bih and smooth_bih_udrho to their sqrt 
  ! form so to apply to both parts of the iterated Laplacian 
  smooth_bih(:,:)       = sqrt(smooth_bih(:,:))
  smooth_bih_udrho(:,:) = sqrt(smooth_bih_udrho(:,:))

  ! register time dependent fields for diagnostic manager 

  if(vert_coordinate_class==DEPTH_BASED) then 
     model_type=' [Boussinesq (volume conserving) model]'
  else 
     model_type=' [non-Boussinesq (mass conserving) model]'
  endif 


  id_sea_level = register_diag_field ('ocean_model', 'sea_level', Grd%tracer_axes(1:2),   &
        Time%model_time, 'effective sea level (eta_t + patm/(rho0*g)) on T cells','meter',&
        missing_value=missing_value, range=(/-1e3,1e3/),                                  &
        standard_name='sea_surface_height_above_geoid' )

  id_sea_level_sq = register_diag_field ('ocean_model', 'sea_level_sq', Grd%tracer_axes(1:2),     &
        Time%model_time, 'square of effective sea level (eta_t + patm/(rho0*g)) on T cells','m^2',&
        missing_value=missing_value, range=(/-1e3,1e3/),                                          &
        standard_name='square_of_sea_surface_height_above_geoid' )

  id_eta_t_mod  = register_diag_field ('ocean_model', 'eta_t_mod', Grd%tracer_axes(1:2), &
              Time%model_time, 'surface height on T cells', 'meter',                     &
              missing_value=missing_value, range=(/-1e3,1e3/))

  id_eta_t  = register_diag_field ('ocean_model', 'eta_t', Grd%tracer_axes(1:2),     &
              Time%model_time, 'surface height on T cells'//trim(model_type),'meter',&
              missing_value=missing_value, range=(/-1e3,1e3/))

  id_eta_t_sq = register_diag_field ('ocean_model', 'eta_t_sq', Grd%tracer_axes(1:2),        &
              Time%model_time, 'square of surface height on T cells'//trim(model_type),'m^2',&
              missing_value=missing_value, range=(/-1e3,1e3/))

  id_eta_t_bar    = register_diag_field ('ocean_model', 'eta_t_bar', Grd%tracer_axes(1:2), &
              Time%model_time, 'surface height on T cells averaged over bar loop', 'meter',&
              missing_value=missing_value, range=(/-1e3,1e3/))

  id_eta_t_bt  = register_diag_field ('ocean_model', 'eta_t_bt', Grd%tracer_axes(1:2),     &
                 Time%model_time, 'surface height on T cells w/i barotropic loop', 'meter',&
                 missing_value=missing_value, range=(/-1e3,1e3/))

  id_deta_dt  = register_diag_field ('ocean_model', 'deta_dt', Grd%tracer_axes(1:2), &
                Time%model_time,'tendency for eta_t', 'm/s',                         &
                missing_value=missing_value, range=(/-1e6,1e6/))

  id_eta_u  = register_diag_field ('ocean_model', 'eta_u', Grd%vel_axes_uv(1:2), &
              Time%model_time,'surface height on U cells', 'meter',              &
              missing_value=missing_value, range=(/-1e3,1e3/))


  id_pbot_u = register_diag_field ('ocean_model', 'pbot_u', Grd%vel_axes_uv(1:2), &
              Time%model_time,'bottom pressure on U cells', 'dbar',               &
              missing_value=missing_value, range=(/-1e6,1e6/))

  id_pbot_t = register_diag_field ('ocean_model', 'pbot_t', Grd%tracer_axes(1:2),   &
              Time%model_time,'bottom pressure on T cells'//trim(model_type),'dbar',&
              missing_value=missing_value, range=(/-1e6,1e6/),                    &
              standard_name='sea_water_pressure_at_sea_floor')

  id_anompb = register_diag_field ('ocean_model', 'anompb', Grd%tracer_axes(1:2), &
              Time%model_time,'T-cell bottom pressure - rho0*grav*ht', 'dbar',    &
              missing_value=missing_value, range=(/-1e6,1e6/))

  id_anompb_bt = register_diag_field ('ocean_model', 'anompb_bt', Grd%tracer_axes(1:2),&
                 Time%model_time,'pbot_t - rho0*grav*ht in barotropic loop', 'dbar',   &
                 missing_value=missing_value, range=(/-1e6,1e6/))

  id_dpbot_dt = register_diag_field ('ocean_model', 'dpbot_dt', Grd%tracer_axes(1:2), &
                Time%model_time,'tendency for pbot on T cells', 'Pa/s',               &
                missing_value=missing_value, range=(/-1e6,1e6/))

  id_patm_t   = register_diag_field ('ocean_model', 'patm_t', Grd%tracer_axes(1:2), &
                Time%model_time,'applied pressure on T cells', 'Pa',                &
                missing_value=missing_value, range=(/-1e6,1e6/),                    &
                standard_name='sea_water_pressure_at_sea_water_surface')

  id_patm_u   = register_diag_field ('ocean_model', 'patm_u', Grd%vel_axes_uv(1:2), &
                Time%model_time,'applied pressure on U cells', 'Pa',                &
                missing_value=missing_value, range=(/-1e6,1e6/))

  id_patm_for_sea_lev = register_diag_field('ocean_model','patm_for_sea_lev', Grd%tracer_axes(1:2),&
       Time%model_time, 'applied pressure used to compute sea level for sea ice dynamics',         &
       'Pa', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_sea_lev_for_coupler = register_diag_field('ocean_model','sea_lev_for_coupler', Grd%tracer_axes(1:2),&
       Time%model_time, 'sea level passed to coupler for use in computing sea ice dynamics',             &
       'm', missing_value=missing_value,range=(/-1.e10,1.e10/))   

  id_dpatm_dt = register_diag_field ('ocean_model', 'dpatm_dt', Grd%tracer_axes(1:2), &
                Time%model_time,'tendency for patm on T cells', 'Pa/s',               &
                missing_value=missing_value, range=(/-1e6,1e6/))

  id_ps   = register_diag_field ('ocean_model', 'ps', Grd%tracer_axes(1:2),  &
            Time%model_time,'surf+atmos pressure', 'Pascal',                 &
            missing_value=missing_value, range=(/-1e5,1e5/))

  id_psx  = register_diag_field ('ocean_model', 'psx', Grd%vel_axes_uv(1:2), &
            Time%model_time,'zonal surf pressure gradient', 'Pa/m',          &
            missing_value=missing_value, range=(/-1e9,1e9/))

  id_psy  = register_diag_field ('ocean_model', 'psy', Grd%vel_axes_uv(1:2), &
            Time%model_time, 'meridional surf pressure gradient', 'Pa/m',    &
            missing_value=missing_value, range=(/-1e9,1e9/))

  id_pb   = register_diag_field ('ocean_model', 'pb', Grd%tracer_axes(1:2), &
            Time%model_time,'effective bottom pressure ', 'Pa',             &
            missing_value=missing_value, range=(/-1e9,1e9/))

  id_grad_anompbx  = register_diag_field ('ocean_model', 'grad_anompbx', Grd%vel_axes_uv(1:2), &
            Time%model_time,'zonal bottom pressure force', 'Pa/m',                             &
            missing_value=missing_value, range=(/-1e9,1e9/))

  id_grad_anompby  = register_diag_field ('ocean_model', 'grad_anompby', Grd%vel_axes_uv(1:2), &
            Time%model_time, 'meridional bottom pressure force', 'Pa/m',                       &
            missing_value=missing_value, range=(/-1e9,1e9/))

  id_urhod = register_diag_field ('ocean_model', 'urhod', Grd%vel_axes_uv(1:2), &
             Time%model_time,'depth and density weighted u', '(kg/m^3)*m^2/s',  &
             missing_value=missing_value, range=(/-1e6,1e6/))

  id_vrhod   = register_diag_field ('ocean_model', 'vrhod', Grd%vel_axes_uv(1:2), &
               Time%model_time,'depth and density weighted v', '(kg/m^3)*m^2/s',  &
               missing_value=missing_value, range=(/-1e6,1e6/))

  id_psiu   = register_diag_field ('ocean_model', 'psiu', Grd%tracer_axes_flux_x(1:2),   &
              Time%model_time,'quasi-barotropic strmfcn psiu (compatible with tx_trans)',&
              trim(transport_dims), missing_value=missing_value, range=(/-1e10,1e10/),   &
              standard_name='ocean_barotropic_mass_streamfunction')

  id_psiv   = register_diag_field ('ocean_model', 'psiv', Grd%tracer_axes_flux_y(1:2),   &
              Time%model_time,'quasi-barotropic strmfcn psiv (compatible with ty_trans)',&
              trim(transport_dims), missing_value=missing_value, range=(/-1e10,1e10/))

  id_conv_rho_ud_t = register_diag_field ('ocean_model', 'conv_rho_ud_t', Grd%tracer_axes(1:2), &
                     Time%model_time,'convergence rho*ud on T cells', '(kg/m^3)*(m/s)',         &
                     missing_value=missing_value, range=(/-1e3,1e3/))

  id_eta_smoother = register_diag_field ('ocean_model', 'eta_smoother', Grd%tracer_axes(1:2), &
                  Time%model_time,'surface smoother applied to eta', 'm/s',                   &
                  missing_value=missing_value, range=(/-1e6,1e6/))  

  id_pbot_smoother  = register_diag_field ('ocean_model', 'pbot_smooth',  Grd%tracer_axes(1:2), &
                    Time%model_time,'bottom smoother applied to bottom pressure', 'dbar/s',     &
                    missing_value=missing_value, range=(/-1e6,1e6/))  

  id_smooth_mask  = register_diag_field ('ocean_model', 'smooth_mask', Grd%tracer_axes(1:2), &
                    Time%model_time,'mask for eta/pbot smoothing', 'dimensionless',          &
                    missing_value=missing_value, range=(/-1e1,1e1/))  

  id_ke_bt = register_diag_field ('ocean_model', 'ke_bt',      &
             Time%model_time,'barotropic ke', '10^-15 joules', &
             missing_value=missing_value, range=(/-10.0,1000.0/))

  id_pe_bt = register_diag_field ('ocean_model', 'pe_bt',      &
             Time%model_time,'barotropic pe', '10^-15 joules', &
             missing_value=missing_value, range=(/-10.0,1000.0/))

  id_pme_total = register_diag_field ('ocean_model', 'pme_total',       &
                 Time%model_time,'mass pme added over dtuv time', 'kg', &
                 missing_value=missing_value, range=(/-1e18,1e18/))

  id_river_total = register_diag_field ('ocean_model', 'river_total',       &
                   Time%model_time,'mass river added over dtuv time', 'kg', &
                   missing_value=missing_value, range=(/-10.0,1e18/))

  id_etat_avg  = register_diag_field ('ocean_model', 'etat_avg', &
                 Time%model_time, 'area averaged eta_t','m',     &
                 missing_value=missing_value, range=(/-1e3,1e3/))

  id_mass_total = register_diag_field ('ocean_model', 'mass_total', &
                  Time%model_time, 'total ocean mass','kg',         &
                  missing_value=missing_value, range=(/0.0,1e12/))

  id_forcing_u_bt = register_diag_field ('ocean_model', 'forcing_u_bt', Grd%vel_axes_uv(1:2), &
                  Time%model_time,'i-forcing of barotropic system, including smf and bmf',    &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e9,1e9/))

  id_forcing_v_bt = register_diag_field ('ocean_model', 'forcing_v_bt', Grd%vel_axes_uv(1:2), &
                  Time%model_time,'j-forcing of barotropic system, including smf and bmf',    &
                  '(kg/m^3)*(m^2/s^2)', missing_value=missing_value, range=(/-1e9,1e9/))

  id_nonlin_forcing_u_bt = register_diag_field ('ocean_model', 'nonlin_forcing_u_bt',            &
                  Grd%vel_axes_uv(1:2), Time%model_time,                                         &
                  'i-forcing of barotropic system, excluding smf and bmf', '(kg/m^3)*(m^2/s^2)', &
                  missing_value=missing_value, range=(/-1e9,1e9/))
  id_nonlin_forcing_v_bt = register_diag_field ('ocean_model', 'nonlin_forcing_v_bt',            &
                  Grd%vel_axes_uv(1:2), Time%model_time,                                         &
                  'j-forcing of barotropic system, excluding smf and bmf', '(kg/m^3)*(m^2/s^2)', &
                  missing_value=missing_value, range=(/-1e9,1e9/))

  id_ext_mode_source = register_diag_field ('ocean_model', 'ext_mode_source', Grd%tracer_axes(1:2),      &
                       Time%model_time,'vertically integrated mass source on T cells', '(kg/m^3)*(m/s)', &
                       missing_value=missing_value, range=(/-1e9,1e9/))

  id_eta_nonbouss = register_diag_field ('ocean_model', 'eta_nonbouss', Grd%tracer_axes(1:2), &
                    Time%model_time, 'surface height including steric contribution', 'meter', &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_nonsteric = register_diag_field ('ocean_model', 'eta_nonsteric', Grd%tracer_axes(1:2), &
                   Time%model_time, 'eta_nonbouss arising from dynamics+water+source', 'meter', &
                   missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_steric = register_diag_field ('ocean_model', 'eta_steric', Grd%tracer_axes(1:2), &
                  Time%model_time, 'eta_nonbouss arising from steric effects', 'meter',   &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_water  = register_diag_field ('ocean_model', 'eta_water', Grd%tracer_axes(1:2), &
                  Time%model_time, 'eta arising from boundary fluxes of water', 'meter', &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_dynamic = register_diag_field ('ocean_model', 'eta_dynamic', Grd%tracer_axes(1:2),            &
                  Time%model_time, 'eta from convergence of mass(nonBouss) or volume(Bouss)', 'meter', &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_source = register_diag_field ('ocean_model', 'eta_source', Grd%tracer_axes(1:2), &
                  Time%model_time, 'eta arising from source term', 'meter',               &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_surf_temp = register_diag_field ('ocean_model', 'eta_surf_temp', Grd%tracer_axes(1:2), &
                    Time%model_time, 'eta arising from surface temp forcing', 'meter',          &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_surf_salt = register_diag_field ('ocean_model', 'eta_surf_salt', Grd%tracer_axes(1:2), &
                    Time%model_time, 'eta arising from surface salt forcing', 'meter',          &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_surf_water = register_diag_field ('ocean_model', 'eta_surf_water', Grd%tracer_axes(1:2), &
                    Time%model_time, 'eta arising from surface water forcing', 'meter',           &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_bott_temp = register_diag_field ('ocean_model', 'eta_bott_temp', Grd%tracer_axes(1:2), &
                    Time%model_time, 'eta arising from bottom temp forcing', 'meter',           &
                    missing_value=missing_value, range=(/-1e4,1e4/))

  id_eta_stemp_forcing = register_diag_field ('ocean_model', 'eta_stemp_forcing', Grd%tracer_axes(1:2), &
                    Time%model_time, 'surface temperature forcing for eta', 'meter',                    &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_btemp_forcing = register_diag_field ('ocean_model', 'eta_btemp_forcing', Grd%tracer_axes(1:2), &
                    Time%model_time, 'bottom temperature forcing for eta', 'meter',                     &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_salt_forcing = register_diag_field ('ocean_model', 'eta_salt_forcing', Grd%tracer_axes(1:2), &
                    Time%model_time, 'surface salt forcing for eta', 'meter',                         &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_water_forcing = register_diag_field ('ocean_model', 'eta_water_forcing', Grd%tracer_axes(1:2), &
                    Time%model_time, 'surface water forcing for eta', 'meter',                          &
                    missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_dpress_dt = register_diag_field ('ocean_model', 'eta_dpress_dt', Grd%tracer_axes(1:2), &
                    Time%model_time, 'vertical integral of grav*w*drho_dpress', 'm/s',          &
                    missing_value=missing_value, range=(/-1e4,1e4/))

  id_eta_stemp_forcing_global = register_diag_field ('ocean_model', 'eta_stemp_forcing_global', &
                  Time%model_time, 'global ave surface temperature forcing for eta', 'meter',   & 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_btemp_forcing_global = register_diag_field ('ocean_model', 'eta_btemp_forcing_global', &
                  Time%model_time, 'global ave bottom temperature forcing for eta', 'meter',   & 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_salt_forcing_global = register_diag_field ('ocean_model', 'eta_salt_forcing_global', &
                  Time%model_time, 'global ave surface salt forcing for eta', 'meter',        & 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_water_forcing_global = register_diag_field ('ocean_model', 'eta_water_forcing_global', &
                  Time%model_time, 'global ave surface water forcing for eta', 'meter',     & 
                  missing_value=missing_value, range=(/-1e4,1e4/))


  id_eta_global = register_diag_field ('ocean_model', 'eta_global',                 &
                  Time%model_time, 'global ave eta_t plus patm_t/(g*rho0)', 'meter',& 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_nonbouss_global = register_diag_field ('ocean_model', 'eta_nonbouss_global', &
                  Time%model_time, 'global ave eta_nonbouss', 'meter',                & 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_steric_global = register_diag_field ('ocean_model', 'eta_steric_global', &
                  Time%model_time, 'global ave eta_steric', 'meter',              &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_nonsteric_global = register_diag_field ('ocean_model', 'eta_nonsteric_global', &
                  Time%model_time, 'global ave eta_nonsteric', 'meter',                 &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_dynamic_global = register_diag_field ('ocean_model', 'eta_dynamic_global', &
                  Time%model_time, 'global ave eta_dynamic', 'meter',               &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_water_global = register_diag_field ('ocean_model', 'eta_water_global', &
                  Time%model_time, 'global ave eta_water', 'meter',             &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_source_global = register_diag_field ('ocean_model', 'eta_source_global', &
                  Time%model_time, 'global ave eta_source', 'meter',              &
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_surf_temp_global = register_diag_field ('ocean_model', 'eta_surf_temp_global',    &
                  Time%model_time, 'global ave eta forced by surface temp forcing', 'meter',& 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_surf_salt_global = register_diag_field ('ocean_model', 'eta_surf_salt_global',     &
                  Time%model_time, 'global ave eta forced by surface salt forcing', 'meter',& 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_surf_water_global = register_diag_field ('ocean_model', 'eta_surf_water_global',    &
                  Time%model_time, 'global ave eta forced by surface water forcing', 'meter',& 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_bott_temp_global = register_diag_field ('ocean_model', 'eta_bott_temp_global',    &
                  Time%model_time, 'global ave eta forced by bottom temp forcing', 'meter',& 
                  missing_value=missing_value, range=(/-1e6,1e6/))

  id_eta_tend  = register_diag_field ('ocean_model', 'eta_tend', Grd%tracer_axes(1:2), &
                 Time%model_time, 'change in eta_nonbouss over a time step', 'meter',  &
                 missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_dynamic_tend  = register_diag_field ('ocean_model', 'eta_dynamic_tend', Grd%tracer_axes(1:2), &
              Time%model_time, 'dynamical contributions to eta_nonbouss over a time step', 'meter',    &
              missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_water_tend = register_diag_field ('ocean_model', 'eta_water_tend', Grd%tracer_axes(1:2),&
              Time%model_time, 'water contributions to eta_nonbouss over a time step', 'meter',  &
              missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_steric_tend = register_diag_field ('ocean_model', 'eta_steric_tend', Grd%tracer_axes(1:2), &
              Time%model_time, 'steric contributions to eta_nonbouss over a time step', 'meter',    &
              missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_steric_tend_A = register_diag_field ('ocean_model', 'eta_steric_tend_A', Grd%tracer_axes(1:2),&
              Time%model_time, 'thickness ratio for steric contributions', 'dimensionless',            &
              missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_steric_tend_B = register_diag_field ('ocean_model', 'eta_steric_tend_B', Grd%tracer_axes(1:2),&
              Time%model_time, 'density ratio for steric contributions', 'dimensionless',              &
              missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_nonsteric_tend = register_diag_field ('ocean_model', 'eta_nonsteric_tend', Grd%tracer_axes(1:2), &
              Time%model_time, 'nonsteric contributions to eta_nonbouss over a time step', 'meter',       &
              missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_source_tend = register_diag_field ('ocean_model', 'eta_source_tend', Grd%tracer_axes(1:2), &
              Time%model_time, 'source contributions to eta_nonbouss over a time step', 'meter',    &
              missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_smooth_tend = register_diag_field ('ocean_model', 'eta_smooth_tend', Grd%tracer_axes(1:2), &
              Time%model_time, 'smoothing  contributions to on diagnosed eta used for eta_nonbouss',&
              'meter', missing_value=missing_value, range=(/-1e4,1e4/))

  id_eta_tend_global = register_diag_field ('ocean_model', 'eta_tend_global',                    &
                  Time%model_time, 'global ave change in eta_nonbouss over a time step', 'meter',& 
                  missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_dynamic_tend_global = register_diag_field ('ocean_model', 'eta_dynamic_tend_global',            &
                  Time%model_time, 'global ave dynamical contributions to eta_nonbouss over a time step',&
                  'meter', missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_water_tend_global = register_diag_field ('ocean_model', 'eta_water_tend_global',               &
                  Time%model_time, 'global ave eustatic contributions to eta_nonbouss over a time step',&
                  'meter', missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_steric_tend_global = register_diag_field ('ocean_model', 'eta_steric_tend_global',           &
                  Time%model_time, 'global ave steric contributions to eta_nonbouss over a time step',&
                  'meter', missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_source_tend_global = register_diag_field ('ocean_model', 'eta_source_tend_global',           &
                  Time%model_time, 'global ave source contributions to eta_nonbouss over a time step',&
                  'meter', missing_value=missing_value, range=(/-1e4,1e4/))
  id_eta_smooth_tend_global = register_diag_field ('ocean_model', 'eta_smooth_tend_global',             &
                  Time%model_time, 'global ave smoother contributions to eta_nonbouss over a time step',&
                  'meter', missing_value=missing_value, range=(/-1e4,1e4/))


  id_udrho_bt_lap = register_diag_field ('ocean_model', 'udrho_bt_lap',                 &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'laplacian friction  to i-component udrho on barotropic time step', &
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_udrho_bt_bih = register_diag_field ('ocean_model', 'udrho_bt_bih',                 &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'biharmonic friction  to i-component udrho on barotropic time step',&
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_udrho_lap    = register_diag_field ('ocean_model', 'udrho_lap',                    &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'laplacian friction  to i-component udrho on baroclinic time step', &
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_udrho_bih    = register_diag_field ('ocean_model', 'udrho_bih',                    &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'biharmonic friction  to i-component udrho on baroclinic time step',&
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))

  id_vdrho_bt_lap = register_diag_field ('ocean_model', 'vdrho_bt_lap',                 &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'laplacian friction  to j-component vdrho on barotropic time step', &
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_vdrho_bt_bih = register_diag_field ('ocean_model', 'vdrho_bt_bih',                 &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'biharmonic friction  to j-component vdrho on barotropic time step',&
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_vdrho_lap    = register_diag_field ('ocean_model', 'vdrho_lap',                    &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'laplacian friction  to j-component vdrho on baroclinic time step', &
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))
  id_vdrho_bih    = register_diag_field ('ocean_model', 'vdrho_bih',                    &
                    Grd%vel_axes_uv(1:2), Time%model_time,                              &
                    'biharmonic friction  to j-component vdrho on baroclinic time step',&
                    '(m^2/s^2)*(kg/m^3)', missing_value=missing_value, range=(/-1e10,1e10/))

  id_rhobarz = register_diag_field ('ocean_model', 'rhobarz', Grd%tracer_axes(1:2), &
              Time%model_time, 'vertically averaged in-situ density', 'kg/m^3',     &
              missing_value=missing_value, range=(/-1e1,1e6/))

  id_rhoavg  = register_diag_field ('ocean_model', 'rhoavg',                                &
               Time%model_time, 'global averaged in-situ density from ocean_barotropic_mod',&
               'kg/m^3', missing_value=missing_value, range=(/-1e1,1e6/))

  if(vert_coordinate_class==PRESSURE_BASED) then 
      rho_water_r=0.001
      id_pme_velocity = register_diag_field('ocean_model','pme_velocity', Grd%tracer_axes(1:2),        &
           Time%model_time, 'net (precip-evap)(kg/(m2*sec) into ocean, divided by 1000kg/m3', 'm/sec', &
           missing_value=missing_value,range=(/-1e6,1e6/))   
  else 
      rho_water_r=rho0r
      id_pme_velocity = register_diag_field('ocean_model','pme_velocity', Grd%tracer_axes(1:2),       &
           Time%model_time, 'net (precip-evap)(kg/(m2*sec) into ocean, divided by 1035kg/m3', 'm/sec',&
           missing_value=missing_value,range=(/-1e6,1e6/))   
  endif


  if (.NOT. file_exist('INPUT/ocean_barotropic.res.nc')) then

      if (.NOT.Time%init) then
           call mpp_error(FATAL,'Expecting file INPUT/ocean_barotropic.res.nc to exist.&
           &This file was not found and Time%init=.false.')
      endif 

      if( file_exist('INPUT/eta_t_ic.nc') ) then
          call read_data('INPUT/eta_t_ic.nc', 'eta_t',Ext_mode%eta_t(:,:,1), &
               Dom%domain2d,timelevel=1)
          Ext_mode%eta_t(:,:,1)=Grd%tmask(:,:,1)*Ext_mode%eta_t(:,:,1)
          call mpp_update_domains (Ext_mode%eta_t(:,:,1), Dom%domain2d)
          Ext_mode%eta_t(:,:,2)        = Ext_mode%eta_t(:,:,1)
          Ext_mode%eta_t(:,:,3)        = Ext_mode%eta_t(:,:,1)
          do n=1,3
             eta_t_bt(:,:,n)           = Ext_mode%eta_t(:,:,n)
             Ext_mode%eta_t_bar(:,:,n) = Ext_mode%eta_t(:,:,n)
             Ext_mode%eta_u(:,:,n)     = Grd%umask(:,:,1)*REMAP_BT_TO_BU_LOCAL(Ext_mode%eta_t(:,:,n))
             call mpp_update_domains (Ext_mode%eta_u(:,:,n), Dom%domain2d)
          enddo
      endif
      if( file_exist('INPUT/pbot_t_ic.nc') ) then
          call read_data('INPUT/pbot_t_ic.nc', 'pbot_t',Ext_mode%pbot_t(:,:,1), &
               Dom%domain2d,timelevel=1)
          call mpp_update_domains (Ext_mode%pbot_t(:,:,1), Dom%domain2d)
          Ext_mode%pbot_t(:,:,2)     = Ext_mode%pbot_t(:,:,1)
          Ext_mode%pbot_t(:,:,3)     = Ext_mode%pbot_t(:,:,1)
          do n=1,3
             anompb_bt(:,:,n)           = Ext_mode%pbot_t(:,:,n) - rho0grav*Grd%ht(:,:)
             Ext_mode%anompb_bar(:,:,n) = Ext_mode%pbot_t(:,:,n) - rho0grav*Grd%ht(:,:)
             Ext_mode%anompb(:,:,n)     = Ext_mode%pbot_t(:,:,n) - rho0grav*Grd%ht(:,:)
             Ext_mode%patm_t(:,:,n)     = 0.0
             Ext_mode%pbot_u(:,:,n)     = Grd%umask(:,:,1)*REMAP_BT_TO_BU_LOCAL(Ext_mode%pbot_t(:,:,n))
             call mpp_update_domains (Ext_mode%pbot_u(:,:,n), Dom%domain2d)
          enddo
      endif
      if( file_exist('INPUT/udrho_ic.nc') ) then
          call read_data('INPUT/udrho_ic.nc', 'udrho',Ext_mode%udrho(:,:,1,1), &
               Dom%domain2d,timelevel=1)
          call read_data('INPUT/udrho_ic.nc', 'vdrho',Ext_mode%udrho(:,:,2,1), &
               Dom%domain2d,timelevel=1)
          call mpp_update_domains(Ext_mode%udrho(:,:,1,1), Ext_mode%udrho(:,:,2,1), &
                                  Dom%domain2d, gridtype=BGRID_NE)
          if(have_obc) then
             call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,1),'M','n')
             call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,1),'M','t')
             call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,1),'Z','t')
             call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,1),'Z','n')
          endif

          Ext_mode%udrho(:,:,:,2) = Ext_mode%udrho(:,:,:,1)
          Ext_mode%udrho(:,:,:,3) = Ext_mode%udrho(:,:,:,1)
          udrho_bt(:,:,:,1)       = Ext_mode%udrho(:,:,:,1)
          udrho_bt(:,:,:,2)       = Ext_mode%udrho(:,:,:,1)
          udrho_bt(:,:,:,3)       = Ext_mode%udrho(:,:,:,1)
      endif
      if(have_obc) then
         call ocean_obc_update_boundary(Ext_mode%eta_t(:,:,:), 'T')    
         call ocean_obc_update_boundary(Ext_mode%eta_u(:,:,:), 'C')    
         call ocean_obc_update_boundary(Ext_mode%eta_t_bar(:,:,:), 'T')    
         call ocean_obc_update_boundary(eta_t_bt(:,:,:), 'T')    
      endif
  endif
  call read_barotropic(Time, Ext_mode)

  if(Time%init .and. zero_eta_ic) then 
      call mpp_error(NOTE, &
      '==>Zeroeing initial barotropic fields. This overrides any initial or restart files.')
      Ext_mode%eta_t(:,:,:)         = 0.0
      Ext_mode%conv_rho_ud_t(:,:,:) = 0.0
      Ext_mode%eta_t_bar(:,:,:)     = 0.0
      Ext_mode%eta_u(:,:,:)         = 0.0
      Ext_mode%udrho(:,:,:,:)       = 0.0
      do n=1,3
        Ext_mode%pbot_t(:,:,n)      = rho0*grav*Grd%ht(:,:)
        Ext_mode%pbot_u(:,:,n)      = rho0*grav*Grd%hu(:,:)
        Ext_mode%anompb(:,:,n)      = 0.0
        Ext_mode%anompb_bar(:,:,n)  = 0.0
        Ext_mode%patm_t(:,:,n)      = 0.0
      enddo
  endif

  if(Time%init .and. ideal_initial_eta) then 
      call mpp_error(NOTE, &
      '==>Initializing eta_t to internally generated ideal profile.')
      call ideal_initialize_eta
      do n=1,3
         Ext_mode%eta_t(:,:,n) = eta_t_init(:,:)
         Ext_mode%eta_u(:,:,n) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%eta_t(:,:,n))
      enddo
  endif

  if (zero_tendency) then
     call mpp_error(NOTE,&
     'Using zero_tendency in ocean_barotropic_mod.  Will freeze external mode fields')
      Ocean_options%barotropic_tendency = 'Did not time step the barotropic fields.'   
  else
      Ocean_options%barotropic_tendency = 'Time stepped the barotropic fields.'   
  endif 
  if (zero_forcing_bt) then
     call mpp_error(NOTE,&
     'Using zero_forcing_bt in ocean_barotropic_mod.  Will set Ext_mode%forcing_bt=0.0.')
  endif 
  if (zero_nonlinear_forcing_bt) then
     call mpp_error(NOTE,&
     'Using zero_nonlinear_forcing_bt in ocean_barotropic_mod: only force with smf and bmf')
  endif 

  if(udrho_bt_lap) then 
       call mpp_error(NOTE,&
       'udrho_bt_lap in ocean_barotropic_mod: smoothing udrho on each barotropic time step.')
       smooth_methods=smooth_methods+1
  endif
  if(udrho_bt_bih) then 
       call mpp_error(NOTE,&
       'udrho_bt_bih in ocean_barotropic_mod: smoothing udrho on each barotropic time step.')
       smooth_methods=smooth_methods+1
  endif 
  if(udrho_lap) then 
       call mpp_error(NOTE,&
       'udrho_lap in ocean_barotropic_mod: smoothing udrho on each baroclinic time step.')
       smooth_methods=smooth_methods+1
  endif 
  if(udrho_bih) then 
       call mpp_error(NOTE,&
       'udrho_bih in ocean_barotropic_mod: smoothing udrho on each baroclinic time step.')
       smooth_methods=smooth_methods+1
  endif

  if(barotropic_bmf) then 
     call mpp_error(NOTE,&
     '==>ocean_barotropic_mod: barotropic_bmf=.true. computes bottom drag on each barotropic step.')
  endif


  if (barotropic_leap_frog) then
     if(vert_coordinate_class==PRESSURE_BASED) then 
       call mpp_error(NOTE,&
       'barotropic_leap_frog w/ press-based coord is not implemented. Use barotropic_pred_corr.')
     endif  
     call mpp_error(NOTE,&
     'Using barotropic_leap_frog=.true. for integrating barotropic dynamics.')
     write (stdoutunit,'(/a,f6.2)')&
     ' Robert-Asselin filter on barotropic dynamics has value= ',robert_asselin_bt
     Time_steps%barotropic_scheme = 'barotropic_leap_frog'
  endif 
  if (barotropic_pred_corr) then
     call mpp_error(NOTE,&
     'Using barotropic_pred_corr=.true. for integrating barotropic dynamics.')
     write (stdoutunit,'(/a,f6.2)')&
     ' Predictor-Corrector time filter on barotropic dynamics has value= ',pred_corr_gamma
     Time_steps%barotropic_scheme = 'barotropic_pred_corr'
  endif 
  if (barotropic_pred_corr .and. barotropic_leap_frog) then
   call mpp_error(FATAL,&
   '==>ocean_barotropic_mod: barotropic_pred_corr & barotropic_leap_frog cannot both be true.')
  endif 
  if (.not. barotropic_pred_corr .and. .not. barotropic_leap_frog) then
   call mpp_error(FATAL,&
   '==>ocean_barotropic_mod: barotropic_pred_corr & barotropic_leap_frog cannot both be false.')
  endif 

  if(barotropic_time_stepping_mom4p0) then 
      write (stdoutunit,'(/a)') &
      ' Updating eta_t or pbot_t using a big time step as in mom4p0. '
      initsum_with_bar = initsum_with_bar_mom4p0
      if(initsum_with_bar_mom4p0) then
        write (stdoutunit,'(/a)') &
      ' Initialise sum of barotropic sea level with last average. This is NOT the default. '
      else
        write (stdoutunit,'(/a)') &
      ' Initialise sum of barotropic sea level with eta_t or pbot_t. This is the default. '
      endif
  endif
  if(barotropic_time_stepping_mom4p1) then  
      write (stdoutunit,'(/a)') &
      ' Updating eta_t or pbot_t with a time average assuming barotropic_pred_corr=.true. '
      if(.not. barotropic_pred_corr) then 
          call mpp_error(FATAL,&
          '==>ocean_barotropic_mod: barotropic_pred_corr must be true if using barotropic_time_stepping_mom4p1.')
      endif
      initsum_with_bar = initsum_with_bar_mom4p1
      if(initsum_with_bar_mom4p1) then
        write (stdoutunit,'(/a)') &
      ' Initialise sum of barotropic sea level with last average. This is the default. '
      else
        write (stdoutunit,'(/a)') &
      ' Initialise sum of barotropic sea level with eta_t or pbot_t. This is NOT the default. '
      endif
  endif
  if(barotropic_time_stepping_mom4p0 .and. barotropic_time_stepping_mom4p1) then 
     call mpp_error(FATAL,&
     '==>ocean_barotropic_mod: barotropic_time_stepping_mom4p0 and barotropic_time_stepping_mom4p1 cannot both be true.')
  endif 
  if(.not. barotropic_time_stepping_mom4p0 .and. .not. barotropic_time_stepping_mom4p1) then 
     call mpp_error(FATAL,&
     '==>ocean_barotropic_mod: barotropic_time_stepping_mom4p0 and barotropic_time_stepping_mom4p1 cannot both be false.')
  endif 


  ! timestep consistancy checks
  if (nint(dtuv) > 0) then

    if (nint(dtbt) > nint(dtuv)) then
      call mpp_error(FATAL, &
      '==>Error from ocean_barotropic_mod: dtbt > dtuv is not allowed.')

    elseif (nint(dtbt) == 0) then
      call mpp_error(FATAL,&
      '==>Error from ocean_barotropic_mod: dtbt=0 when dtuv > 0 is not implemented.')

    elseif (nint(dtuv) == nint(dtbt)) then

      call mpp_error(NOTE,&
      '==>From ocean_barotropic_mod: running an unsplit system with dtbt=dtuv.')

      splitting = .false.

      nts = 1
      twodt = (2.0*dtuv)
      write (stdoutunit,'(/a)')' From ocean_barotropic_mod: un-split momentum equations &
                               &imply one barotropic timestep per baroclinic step.    &
                               &External mode solver for udrho is NOT used. '
    else

      splitting = .true.

      nts  = nint((2.0*dtuv)/dtbt)
      crit = 1.e-6
      if (abs((2.0*dtuv)/nts - dtbt) > crit) then
        write (unit,'(/1x,a/)') &
        '==>Error: Splitting momentum equations requires mod(2*dtuv,dtbt)=0; it is not.'
        write (unit,*)      '          Current setting: '
        write (unit,*)      '             dtuv =  ',dtuv,' sec.'
        write (unit,*)      '             dtbt =  ',dtbt,' sec.'
        write (unit,*)      '          New setting should be:  '
        write (unit,*)      '             dtbt =  ',(2.0*dtuv)/nts,' sec.'
        call mpp_error(FATAL,'==>Error from ocean_barotropic_mod: badly chosen time steps')
      endif

      if(barotropic_leap_frog) twodt = 2.0*(2.0*dtuv)/nts

      write (stdoutunit,'(/1x,a,i4,a)')' ==> Note: The barotropic dynamics integrate ',&
                                      nts, ' timesteps for every one baroclinic timestep.'
    endif
  endif

  ! Determine maximum timestep allowable to ensure barotropic 
  ! time stepping is CFL stable with external gravity waves

  cgmax            = 1.0
  icg              = isc
  jcg              = jsc
  gridmin          = 1.0e20
  max_dt_for_cgext = 1.e6

  if(barotropic_leap_frog) cfl_grid_factor=0.5
  if(barotropic_pred_corr) cfl_grid_factor=1.0
  do j=jsc,jec
     do i=isc,iec
        if (Grd%kmu(i,j) > 0) then
            cgext  = sqrt(grav*(Grd%zw(Grd%kmu(i,j))))
            gridsp = 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j))
            gridsp = sqrt(1.0/gridsp) 
            dtcg   = cfl_grid_factor*gridsp/cgext
            if (dtcg < max_dt_for_cgext) then
                gridmin = gridsp; cgmax = cgext
                max_dt_for_cgext = dtcg; icg = i; jcg = j
            endif
        endif
     enddo
  enddo

  if (nint(dtbt) > 0) then

    cgmax = nint(cgmax)
    max_dt_for_cgext = int(max_dt_for_cgext)

    fudge = 1 + 1.e-12*mpp_pe()  ! to distinguish redundancies between processors
    max_dt_for_cgext  = max_dt_for_cgext*fudge
    max_dt_for_cgext0 = max_dt_for_cgext
    call mpp_min (max_dt_for_cgext)
    call mpp_max (cgmax)

    ! show the most unstable location for barotropic gravity waves
    if (max_dt_for_cgext == max_dt_for_cgext0) then

      if (dtbt <= max_dt_for_cgext) then
        write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')                      &
        ' Barotropic stability most nearly violated at U-cell (i,j) = (',&
        icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
        write(unit,'(a,i6)')    '         The number of kmu-levels at this point is ',Grd%kmu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dxu grid spacing (m) at this point is ',Grd%dxu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dyu grid spacing (m) at this point is ',Grd%dyu(icg,jcg)
        cfl_error=.false.
      else
        write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')                   &
        '==>Error: Barotropic stability violated at U-cell (i,j) = (',&
        icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
        write(unit,'(a,i6)')    '         The number of kmu-levels at this point is ',Grd%kmu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dxu grid spacing (m) at this point is ',Grd%dxu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dyu grid spacing (m) at this point is ',Grd%dyu(icg,jcg) 
        cfl_error=.true.  
      endif

      write (unit,'(a,f5.1,a/a,f8.3,a,f8.3,a)')                             &
      '         where the barotropic gravity wave speed is ~',cgmax,' m/s.',&
      '         "dtbt" must be less than ',max_dt_for_cgext/fudge,' sec.   dtbt = ',dtbt,' sec.'

    endif

    if(cfl_error) then 
        call mpp_error(FATAL,'==>Error from ocean_barotropic_mod: time step instability &
                              &detected for external gravity waves. Time step too large.' )
    endif 

    max_dt_for_cgext = nint(max_dt_for_cgext)

  endif

  if(vert_coordinate_class==DEPTH_BASED) then 

      if(smooth_eta_t_bt_laplacian) then 
          write(stdoutunit,'(/1x,a)') &
          '==>smooth_eta_t_bt_laplacian smooths eta_t_bt on each barotropic time step (dtbt).'
      endif
      if(smooth_eta_t_bt_biharmonic) then 
          write(stdoutunit,'(/1x,a)') &
          '==>smooth_eta_t_bt_biharmonic acts on eta_t_bt each dtbt time step. Can produce extrema. Not recommended.'
      endif
      if(smooth_eta_t_bt_laplacian .and. smooth_eta_t_bt_biharmonic) then 
          call mpp_error(FATAL, &
          '==>ocean_barotropic_mod: smooth_eta_t_bt_laplacian & smooth_eta_t_bt_biharmonic cannot both be true.')
      endif

      if(smooth_eta_t_laplacian) then 
          write(stdoutunit,'(/1x,a)') &
          '==>Using smooth_eta_t_laplacian to smooth eta_t.'
      endif
      if(smooth_eta_t_biharmonic) then 
          write(stdoutunit,'(/1x,a)') &
          '==>Using smooth_eta_t_biharmonic to smooth eta_t. Can produce extrema. Not recommended.'
      endif
      if(smooth_eta_t_laplacian .and. smooth_eta_t_biharmonic) then 
          call mpp_error(FATAL,&
          '==>ocean_barotropic_mod: smooth_eta_t_laplacian & smooth_eta_t_biharmonic cannot both be true.')
      endif

  endif

  if(vert_coordinate_class==PRESSURE_BASED) then 

      if(smooth_anompb_bt_laplacian) then 
          write(stdoutunit,'(/1x,a)') &
               '==>smooth_anompb_bt_laplacian smooths anompb_bt on each barotropic time step (dtbt).'
      endif
      if(smooth_anompb_bt_biharmonic) then 
          write(stdoutunit,'(/1x,a)') &
          '==>smooth_anompb_bt_biharmonic acts on anompb_bt each dtbt time step. Can produce extrema. Not recommended.'
      endif
      if(smooth_anompb_bt_laplacian .and. smooth_anompb_bt_biharmonic) then 
          call mpp_error(FATAL, &
          '==>ocean_barotropic_mod: smooth_anompb_bt_laplacian & smooth_anompb_bt_biharmonic cannot both be true.')
      endif

      if(smooth_pbot_t_laplacian) then 
          write(stdoutunit,'(/1x,a)') &
          '==>Using smooth_pbot_t_laplacian to smooth pbot_t-pbot0 on surface height time step (dteta).'
      endif
      if(smooth_pbot_t_biharmonic) then 
          write(stdoutunit,'(/1x,a)') &
          '==>Using smooth_pbot_t_biharmonic to smooth pbot_t-pbot0 on time step (dteta). Can produce extrema. Not recommended.'
      endif
      if(smooth_pbot_t_laplacian .and. smooth_pbot_t_biharmonic) then 
          call mpp_error(FATAL,&
          '==>ocean_barotropic_mod: smooth_pbot_t_laplacian & smooth_pbot_t_biharmonic cannot both be true.')
      endif
      
      if(have_obc) then 
          write(stdoutunit,'(/1x,a)') &
          '==>ocean_barotropic_mod: OBCs are completely untested for pressure based co-ordinates.'
      endif

  endif

  if(truncate_eta) then 
    write(stdoutunit,'(/a)')'==>Warning: truncate_eta=.true.' 
    write(stdoutunit,'(a,f10.4,a,f10.4)') &
   '   Will enforce dzt(1)+eta_t (m) > ',Grd%dzt(1)*frac_crit_cell_height, ' and eta_t (m) < ', eta_max
  endif 

  call tidal_forcing_init(Time)
  call geoid_forcing_init(Time)

  ! get pelist for psi_compute
  call mpp_get_layout (Dom%domain2d, layout)
  allocate ( pelist_x(layout(1)), pelist_y(layout(2)))
  call mpp_get_domain_components (Dom%domain2d, Domx, Domy)
  call mpp_get_pelist( Domx, pelist_x )
  call mpp_get_pelist( Domy, pelist_y )
  size_x = size(pelist_x(:))
  size_y = size(pelist_y(:))
  this_pe = mpp_pe()
  myid_y = 0    
  do i = 1,size_y
     if (this_pe == pelist_y(i)) then
        myid_y = i
        exit
     endif
  enddo
  if(myid_y == 0) call mpp_error(FATAL,'==>Error in ocean_barotropic psi_compute myid_y')
  myid_x = 0
  do i = 1,size_x
     if (this_pe == pelist_x(i)) then
        myid_x = i
        exit
     endif
  enddo
  if(myid_x == 0) call mpp_error(FATAL,'==>Error in ocean_barotropic psi_compute myid_x')


  if(debug_this_module) then 
      tau   = Time%tau
      taup1 = Time%taup1
      if(tendency==THREE_LEVEL) then 
          write(stdoutunit,*) ' '
          write(stdoutunit,*) 'From ocean_barotropic_mod: initial external mode chksums at tau'
          call write_timestamp(Time%model_time)
          call barotropic_chksum(Ext_mode, tau)
      endif
      write(stdoutunit,*) ' '
      write(stdoutunit,*) 'From ocean_barotropic_mod: initial external mode chksums at taup1'
      call write_timestamp(Time%model_time)
      call barotropic_chksum(Ext_mode, taup1)
  endif

end subroutine ocean_barotropic_init
! </SUBROUTINE> NAME="ocean_barotropic_init"


!#######################################################################
! <SUBROUTINE NAME="eta_and_pbot_update">
!
! <DESCRIPTION>
! Time step the surface height or pbot using a "big" time step. 
! These fields are coincident in time with tracers.  
!
! NOTE: For pbot_t updates, we time step anompb for accuracy, then 
! add rho0grav*ht to get pbot_t. 
!
! </DESCRIPTION>
!
subroutine eta_and_pbot_update (Time, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer  :: i, j
  integer  :: taum1, tau, taup1

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1


  ! update eta_t with a "big leap-frog" or "big forward"
  if(vert_coordinate_class==DEPTH_BASED) then 

      ! use eta_t_bar(taum1) to suppress splitting between eta_u and eta_t 
      if(have_obc .and. barotropic_time_stepping_mom4p0) then 
          do j=jsc,jec
             do i=isc,iec
                Ext_mode%eta_t(i,j,taup1) = Ext_mode%eta_t_bar(i,j,taum1) &
                                            + dtime*Ext_mode%deta_dt(i,j)
             enddo
          enddo
      else 
          do j=jsc,jec
             do i=isc,iec
                Ext_mode%eta_t(i,j,taup1) = Ext_mode%eta_t(i,j,taum1) &
                                            + dtime*Ext_mode%deta_dt(i,j)
             enddo
          enddo
      endif

      ! for debugging 
      if(zero_eta_t) then 
         Ext_mode%eta_t(:,:,taup1) = 0.0
      endif 

      if(have_obc) call ocean_obc_surface_height(Time, Ext_mode, dtime)  
      if(truncate_eta) call eta_truncate(Time, Ext_mode)
      call mpp_update_domains (Ext_mode%eta_t(:,:,taup1), Dom%domain2d)
      if(have_obc) call ocean_obc_update_boundary(Ext_mode%eta_t(:,:,taup1),'T')

  ! update pbot_t with a "big leap-frog" or "big forward"
  elseif(vert_coordinate_class==PRESSURE_BASED) then 

      ! use anompb_bar(taum1) to keep from splitting between pbot_t and pbot_u
      if(have_obc .and. barotropic_time_stepping_mom4p0) then 
          do j=jsc,jec
             do i=isc,iec
                Ext_mode%anompb(i,j,taup1) = Ext_mode%anompb_bar(i,j,taum1) &
                                             + dtime*Ext_mode%dpbot_dt(i,j) 
             enddo
          enddo
      else
          do j=jsc,jec
             do i=isc,iec
                Ext_mode%anompb(i,j,taup1) = Ext_mode%anompb(i,j,taum1) &
                                             + dtime*Ext_mode%dpbot_dt(i,j)
             enddo
          enddo
      endif
      call mpp_update_domains (Ext_mode%anompb(:,:,taup1), Dom%domain2d)
      do j=jsd,jed
         do i=isd,ied
            Ext_mode%pbot_t(i,j,taup1) = Ext_mode%anompb(i,j,taup1) + rho0grav*Grd%ht(i,j)
         enddo
      enddo

      if(have_obc) then 
          call ocean_obc_surface_height(Time, Ext_mode, dtime)   
          call ocean_obc_update_boundary(Ext_mode%pbot_t(:,:,taup1),'T')
      endif

  endif

  ! when dtuv=dtbt, we do not compute the time averaged eta_t_bar and anompb_bar,
  ! but we still use these fields.  So set them equal to eta_t and anompb.
  if(.not. splitting) then 
      if(vert_coordinate_class==DEPTH_BASED) then 
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t(i,j,taup1)
             enddo
          enddo
      elseif(vert_coordinate_class==PRESSURE_BASED) then
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%anompb_bar(i,j,taup1) = Ext_mode%anompb(i,j,taup1)
             enddo
          enddo
      endif
  endif


end subroutine eta_and_pbot_update
! </SUBROUTINE> NAME="eta_and_pbot_update"


!#######################################################################
! <SUBROUTINE NAME="eta_and_pbot_diagnose">
!
! <DESCRIPTION>
! Diagnose surface height or pbot depending on the vertical coordinate.
!
! Note that dzt has been updated to taup1 before this routine is called
! since we have already called update_tcell_thickness.
!
! Also, compute geodepth_zt in this routine.  It is necessary to do 
! do this step here, since for pressure coordinate models we do not
! know eta_t(taup1) until this routine...
!
! </DESCRIPTION>
!
subroutine eta_and_pbot_diagnose (Time, Dens, Thickness, patm, pme, river, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in)    :: patm
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river 

  integer  :: tau, taup1
  integer  :: i,j,k,kbot
  real     :: geopotential_anomaly(0:nk)

  tau   = Time%tau
  taup1 = Time%taup1

  ! diagnose bottom pressure 
  if(vert_coordinate_class==DEPTH_BASED) then       
      Ext_mode%pbot_t(:,:,taup1) = 0.0
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%pbot_t(i,j,taup1) = Ext_mode%pbot_t(i,j,taup1) &
                   +Grd%tmask(i,j,k)*grav*Dens%rho(i,j,k,taup1)*Thickness%dzt(i,j,k)
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            Ext_mode%pbot_t(i,j,taup1) = Ext_mode%pbot_t(i,j,taup1) + Ext_mode%patm_t(i,j,taup1)
         enddo
      enddo
      Ext_mode%anompb(:,:,taup1) = Ext_mode%pbot_t(:,:,taup1) - rho0grav*Grd%ht(:,:)
      Ext_mode%pbot_u(:,:,taup1) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%pbot_t(:,:,taup1))
      call mpp_update_domains (Ext_mode%pbot_u(:,:,taup1), Dom%domain2d)
      if(have_obc) call ocean_obc_update_boundary(Ext_mode%pbot_u(:,:,taup1),'T')        

  ! diagnose surface height...
  ! use three approaches, which are identical in continuum
  elseif(vert_coordinate_class==PRESSURE_BASED) then

      wrk1_2d(:,:)              = 0.0
      Ext_mode%eta_t(:,:,taup1) = 0.0

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%eta_t(i,j,taup1) = Ext_mode%eta_t(i,j,taup1) &
                     + Grd%tmask(i,j,k)*Thickness%dzt(i,j,k)
               wrk1_2d(i,j)              = wrk1_2d(i,j) &
                     -Grd%tmask(i,j,k)*Thickness%dzt(i,j,k)*(rho0r*Dens%rho(i,j,k,taup1)-1.0)
            enddo
         enddo
      enddo
      do j=jsd,jed
         do i=isd,ied
            Ext_mode%eta_t(i,j,taup1) = Ext_mode%eta_t(i,j,taup1) - Grd%ht(i,j)
            wrk1_2d(i,j) = wrk1_2d(i,j) - Grd%ht(i,j) &
             +grav_r*rho0r*(Ext_mode%pbot_t(i,j,taup1) - Ext_mode%patm_t(i,j,taup1))
         enddo
      enddo

      call eta_smooth_diagnosed(Time, Ext_mode%eta_t(:,:,taup1))

      Ext_mode%eta_u(:,:,taup1) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%eta_t(:,:,taup1))
      call mpp_update_domains (Ext_mode%eta_u(:,:,taup1), Dom%domain2d)
      if(have_obc) call ocean_obc_update_boundary(Ext_mode%eta_u(:,:,taup1),'T')        

  endif

  ! set the sea level field to be transferred to coupler in ocean_sbc.F90.
  ! note that the if-test is for backward compatibility...bits change but
  ! gradients of Thickness%sea_lev should remain the same.
  if(vert_coordinate == GEOPOTENTIAL) then
      do j=jsd,jed
         do i=isd,ied
            Thickness%sea_lev(i,j) = Thickness%dzt(i,j,1) + Ext_mode%patm_for_sea_lev(i,j)*grav_rho0r
         enddo
      enddo
  else 
      do j=jsd,jed
         do i=isd,ied
            Thickness%sea_lev(i,j) = Ext_mode%eta_t(i,j,taup1) + Ext_mode%patm_for_sea_lev(i,j)*grav_rho0r
         enddo
      enddo
  endif

  ! compute depth from z=0 to T-point, for use in computing geopotential. 
  do j=jsd,jed
     do i=isd,ied
        Thickness%geodepth_zt(i,j,1) = Thickness%depth_zt(i,j,1) - Ext_mode%eta_t(i,j,taup1)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%geodepth_zt(i,j,k) = Thickness%geodepth_zt(i,j,k-1) + Thickness%dzwt(i,j,k-1)
        enddo
     enddo
  enddo

  ! find max and min of some quantities 
  if (diag_step > 0 .and. mod(Time%itt, diag_step) == 0) then 
     if(vert_coordinate == GEOPOTENTIAL) then
        call eta_check(Time, Ext_mode)
     endif 
     call maximum_convrhoud(Time, Ext_mode)
  endif 

  ! global integrals for diagnostics 
  call barotropic_integrals (Time, Ext_mode, patm, pme, river)

  ! send to diag_manager     

  if (id_eta_t > 0)    used = send_data (id_eta_t, Ext_mode%eta_t(:,:,tau), &
                              Time%model_time, rmask=Grd%tmask(:,:,1),      &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_eta_t_sq > 0) used = send_data (id_eta_t_sq, Ext_mode%eta_t(:,:,tau)*Ext_mode%eta_t(:,:,tau), &
                              Time%model_time, rmask=Grd%tmask(:,:,1),                                 &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_eta_t_mod > 0) used = send_data (id_eta_t_mod, wrk1_2d(:,:),     &
                              Time%model_time, rmask=Grd%tmask(:,:,1),    &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_eta_u > 0)    used = send_data (id_eta_u, Ext_mode%eta_u(:,:,tau), &
                              Time%model_time, rmask=Grd%umask(:,:,1),      &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_pbot_t > 0)   used = send_data (id_pbot_t, c2dbars*Ext_mode%pbot_t(:,:,tau), &
                              Time%model_time, rmask=Grd%tmask(:,:,1),                &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_anompb > 0)   used = send_data (id_anompb, c2dbars*Ext_mode%anompb(:,:,tau), &
                              Time%model_time, rmask=Grd%tmask(:,:,1),                &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_pbot_u > 0)   used = send_data (id_pbot_u, c2dbars*Ext_mode%pbot_u(:,:,tau), &
                              Time%model_time, rmask=Grd%umask(:,:,1),                &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! remove inverse barometer from atmos and sea ice to get effective sea level 
  if (id_sea_level > 0) then 
       used = send_data (id_sea_level,                              &
       Ext_mode%eta_t(:,:,tau)+grav_rho0r*Ext_mode%patm_t(:,:,tau), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                     &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  if (id_sea_level_sq > 0) then 
       wrk1_2d(:,:)= 0.0       
       do j=jsc,jec
          do i=isc,iec
             wrk1_2d(i,j) = (Ext_mode%eta_t(i,j,tau)+grav_rho0r*Ext_mode%patm_t(i,j,tau))**2
          enddo
       enddo
       used = send_data (id_sea_level_sq, wrk1_2d(:,:),&
       Time%model_time, rmask=Grd%tmask(:,:,1),        &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


  if(id_patm_for_sea_lev > 0) then 
      used = send_data(id_patm_for_sea_lev,Ext_mode%patm_for_sea_lev(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),                     &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  if(id_sea_lev_for_coupler > 0) then 
      used = send_data(id_sea_lev_for_coupler,Thickness%sea_lev(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),                &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 



end subroutine eta_and_pbot_diagnose
! </SUBROUTINE> NAME="eta_and_pbot_diagnose"



!#######################################################################
! <SUBROUTINE NAME="eta_and_pbot_tendency">
!
! <DESCRIPTION>
!  Compute tendency for surface height or bottom pressure. 
! </DESCRIPTION>
!
subroutine eta_and_pbot_tendency(Time, pme, river, Ext_mode)
  
  type(ocean_time_type),          intent(in)    :: Time
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: i, j, tau
  tau = Time%tau

  Ext_mode%deta_dt(:,:)  = 0.0
  Ext_mode%dpbot_dt(:,:) = 0.0

  if(vert_coordinate_class==DEPTH_BASED) then 
      do j=jsd,jed
         do i=isd,ied       
            Ext_mode%deta_dt(i,j) = Grd%tmask(i,j,1)                              &
                 *rho0r*( Ext_mode%conv_rho_ud_t(i,j,tau) + pme(i,j) + river(i,j) &
                 +Ext_mode%source(i,j) )
         enddo
      enddo

  elseif(vert_coordinate_class==PRESSURE_BASED) then 
      do j=jsd,jed
         do i=isd,ied       
            Ext_mode%dpbot_dt(i,j) = Ext_mode%dpatm_dt(i,j) + grav*Grd%tmask(i,j,1) &
                 *( Ext_mode%conv_rho_ud_t(i,j,tau) + pme(i,j) + river(i,j)         &
                 +Ext_mode%source(i,j) ) 
         enddo
      enddo

  endif

  ! over-write the above for the zero_tendency test
  if (zero_tendency .or. zero_eta_tendency) then 
     Ext_mode%deta_dt(:,:)  = 0.0
     Ext_mode%dpbot_dt(:,:) = 0.0
  endif 

  if (id_deta_dt > 0)  used = send_data (id_deta_dt, Ext_mode%deta_dt(:,:), &
                              Time%model_time, rmask=Grd%tmask(:,:,1),      &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_dpbot_dt > 0) used = send_data (id_dpbot_dt, Ext_mode%dpbot_dt(:,:), &
                              Time%model_time, rmask=Grd%tmask(:,:,1),        &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! mass flux per area from pme divided by density factor to convert to m/sec
  if (id_pme_velocity > 0) then
      used = send_data(id_pme_velocity, pme(:,:)*rho_water_r, &
             Time%model_time, rmask=Grd%tmask(:,:,1),         &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


end subroutine eta_and_pbot_tendency
! </SUBROUTINE> NAME="eta_and_pbot_tendency"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_barotropic">
!
! <DESCRIPTION>
! Time step the external mode fields using either a leap-frog scheme
! or a predictor-corrector scheme.  Time average these fields to 
! update the vertically integrated density weighted velocity 
! (urhod,vrhod) and the time averaged surface height eta_t_bar 
! or bottom pressure anompb_bar. 
!
! NOTE:   surface pressure gradient and gradient of anomalous 
! bottom pressure are needed for the energy analysis. 
! 
! Also, if splitting=false or update_velocity_via_uprime=.false.,
! use these for velocity update in ocean_velocity_mod. 
!
! Include the tidal forcing if present.  
!
! Include geoid perturbation if present.
! 
! Use time averaged eta and pbot to ensure the most stable pressure 
! gradient for use with splitting=false or update_velocity_via_uprime=.false. 
! 
! </DESCRIPTION>
!
subroutine update_ocean_barotropic (Time, Dens, Thickness, Velocity, Adv_vel, &
                                    Theta, Salinity, Ext_mode, patm, pme, river)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_thickness_type),     intent(in)    :: Thickness 
  type(ocean_velocity_type),      intent(in)    :: Velocity
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_prog_tracer_type),   intent(in)    :: Theta
  type(ocean_prog_tracer_type),   intent(in)    :: Salinity 
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in)    :: patm
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river

  type(time_type)                  :: next_time 
  type(time_type)                  :: time_bt 
  real, dimension(isd:ied,jsd:jed) :: psiu
  real, dimension(isd:ied,jsd:jed) :: psiv
  real                             :: rnts, rntsp1
  integer                          :: tau, taup1, itime 
  integer                          :: i, j, k, n
  integer                          :: day, sec 
  real                             :: dayr    

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1

  ! for setting surface pressure and bottom pressure anomaly 
  if(tendency==TWO_LEVEL) then 
     itime=taup1
  else 
     itime=tau
  endif

  ! for time averaging to get eta_t_bar and anompb_bar 
  rntsp1 = 1.0/(float(nts)+1.0)

  ! for time averaging to get udrho 
  if(barotropic_time_stepping_mom4p0) then 
      rnts = 1.0/float(nts)
  else 
      rnts = 2.0/(float(nts)*float(nts+1))
  endif 

  ! compute inverse thickness of U-cell column
  do j=jsd,jed
     do i=isd,ied
        thicku_r(i,j) = 0.0
        k = Grd%kmu(i,j)
        if(k > 1) then 
            k = Grd%kmu(i,j)
            thicku_r(i,j) = Grd%umask(i,j,1)/(epsln+Thickness%thicku(i,j,tau))  
        endif
     enddo
  enddo

  if(vert_coordinate==GEOPOTENTIAL) then 
      if(tendency==THREE_LEVEL) then 
          do j=jsd,jed
             do i=isd,ied
                rho_g(i,j) = grav*Dens%rho(i,j,1,tau)
             enddo
          enddo
      elseif(tendency==TWO_LEVEL) then 
          do j=jsd,jed
             do i=isd,ied
                rho_g(i,j) = grav*Dens%rho(i,j,1,taup1)
             enddo
          enddo
      endif
  endif

  ! only call external mode time schemes if splitting barotropic from baroclinic
  if(splitting) then 
      
      if(vert_coordinate_class==DEPTH_BASED) then 
        if(barotropic_leap_frog) call leap_frog_barotropic_depth (Time, Thickness, Ext_mode, patm, pme, river)
        if(barotropic_pred_corr) call pred_corr_barotropic_depth (Time, Thickness, Ext_mode, patm, pme, river)
      elseif(vert_coordinate_class==PRESSURE_BASED) then 
        if(barotropic_pred_corr) call pred_corr_barotropic_press (Time, Thickness, Ext_mode, pme, river)
      endif  

      ! time averaged eta_t_bar
      if(vert_coordinate_class==DEPTH_BASED) then
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%eta_t_bar(i,j,taup1) = rntsp1*Ext_mode%eta_t_bar(i,j,taup1) 
             enddo
          enddo
          call mpp_update_domains(Ext_mode%eta_t_bar(:,:,taup1), Dom%domain2d)

      ! time averaged anompb_bar
      elseif(vert_coordinate_class==PRESSURE_BASED) then
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%anompb_bar(i,j,taup1) = rntsp1*Ext_mode%anompb_bar(i,j,taup1) 
             enddo
          enddo
          call mpp_update_domains(Ext_mode%anompb_bar(:,:,taup1), Dom%domain2d)

      endif

      ! time averaged vertically integrated momentum 
      do n=1,2
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%udrho(i,j,n,taup1) = rnts*Ext_mode%udrho(i,j,n,taup1)
            enddo
         enddo
      enddo
      call mpp_update_domains(Ext_mode%udrho(:,:,1,taup1),Ext_mode%udrho(:,:,2,taup1),&
                              Dom%domain2d,gridtype=BGRID_NE) 

      ! convergence of vertically integrated momentum 
      Ext_mode%conv_rho_ud_t(:,:,taup1)  = -DIV_UD(Ext_mode%udrho(:,:,:,taup1))
      if(have_obc) call ocean_obc_adjust_divud(Ext_mode%conv_rho_ud_t(:,:,taup1))        
      call mpp_update_domains (Ext_mode%conv_rho_ud_t(:,:,taup1), Dom%domain2d)


  endif  ! endif for splitting
         ! if not splitting, then full velocity
         ! update is performed in ocean_velocity_mod
  
  
  ! closure for eta_u 
  if(vert_coordinate_class==DEPTH_BASED) then
      Ext_mode%eta_u(:,:,taup1) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%eta_t_bar(:,:,taup1))
      call mpp_update_domains (Ext_mode%eta_u(:,:,taup1), Dom%domain2d)
      if(have_obc) call ocean_obc_update_boundary(Ext_mode%eta_u(:,:,taup1),'T')        

      ! for debugging 
      if(zero_eta_u) then 
         Ext_mode%eta_u(:,:,taup1) = 0.0
      endif 


  ! closure for pbot_u
  elseif(vert_coordinate_class==PRESSURE_BASED) then
      Ext_mode%pbot_u(:,:,taup1) = &
      Grd%umask(:,:,1)*REMAP_BT_TO_BU(rho0grav*Grd%ht(:,:)+Ext_mode%anompb_bar(:,:,taup1))
      call mpp_update_domains (Ext_mode%pbot_u(:,:,taup1), Dom%domain2d)
      if(have_obc) call ocean_obc_update_boundary(Ext_mode%pbot_u(:,:,taup1),'T')        
  endif


  ! overwrite the above for debugging cases run with zero tendency 
  if (nint(dtbt)==0 .or. zero_tendency) then 
      do n=1,2
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%udrho(i,j,n,taup1) = Ext_mode%udrho(i,j,n,tau)
            enddo
         enddo
      enddo
      if(vert_coordinate_class==DEPTH_BASED) then 
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%eta_t(i,j,taup1) = Ext_mode%eta_t(i,j,tau)
                Ext_mode%eta_u(i,j,taup1) = Ext_mode%eta_u(i,j,tau)
             enddo
          enddo
      elseif(vert_coordinate_class==PRESSURE_BASED) then 
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%anompb(i,j,taup1) = Ext_mode%anompb(i,j,tau)
                Ext_mode%pbot_t(i,j,taup1) = Ext_mode%pbot_t(i,j,tau)
                Ext_mode%pbot_u(i,j,taup1) = Ext_mode%pbot_u(i,j,tau)
             enddo
          enddo
      endif
      do j=jsd,jed
         do i=isd,ied
            Ext_mode%conv_rho_ud_t(i,j,taup1) = Ext_mode%conv_rho_ud_t(i,j,tau)
         enddo
      enddo
  endif 

  ! surface pressure gradient and gradient of anomalous bottom pressure.
  ! these fields are needed for the energy analysis. See comments at top
  ! of this subroutine for more details. 
  wrk1_2d(:,:)= 0.0
  if(tidal_forcing) then 
      time_bt = increment_time(Time%model_time, int(dteta),0) 
      call get_time(time_bt, sec, day)                                 
      dayr = day + sec/86400.0  
      call get_tidal_forcing(Time, dayr)
      do j=jsd,jed
         do i=isd,ied
            Ext_mode%ps(i,j) = rho_g(i,j) & 
             *(alphat*Ext_mode%eta_t_bar(i,j,itime)-eta_eq_tidal(i,j)) + patm(i,j)  
            wrk1_2d(i,j)     = Ext_mode%anompb_bar(i,j,itime) & 
                               +rho0grav*(alphat*Ext_mode%eta_t_bar(i,j,taup1)-eta_eq_tidal(i,j))
         enddo
      enddo
  else 
      do j=jsd,jed
         do i=isd,ied
            Ext_mode%ps(i,j) = rho_g(i,j)*Ext_mode%eta_t_bar(i,j,itime) + patm(i,j)  
            wrk1_2d(i,j)     = Ext_mode%anompb_bar(i,j,itime)
         enddo
      enddo
  endif
  if(geoid_forcing) then 
      do j=jsd,jed
         do i=isd,ied
            Ext_mode%ps(i,j) = Ext_mode%ps(i,j) - rho_g(i,j)*eta_geoid(i,j)  
            wrk1_2d(i,j)     = Ext_mode%anompb_bar(i,j,itime) - rho0grav*eta_geoid(i,j)
         enddo
      enddo
  endif

  Ext_mode%grad_ps(:,:,:)     = GRAD_BAROTROPIC_P(Ext_mode%ps(:,:))
  Ext_mode%grad_anompb(:,:,:) = GRAD_BAROTROPIC_P(wrk1_2d(:,:))

  ! diagnose contributions to sea level
  call eta_terms_diagnose(Time, Dens, Thickness, Theta, Salinity, Adv_vel, Ext_mode, pme, river)


  ! send diagnostics to diag_manager     
  if (id_eta_t_bar > 0) used = send_data (id_eta_t_bar, Ext_mode%eta_t_bar(:,:,tau), &
                              Time%model_time, rmask=Grd%tmask(:,:,1),               &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_urhod > 0)     used = send_data (id_urhod, Ext_mode%udrho(:,:,1,tau), &
                              Time%model_time, rmask=Grd%umask(:,:,1),         &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_vrhod > 0)     used = send_data (id_vrhod, Ext_mode%udrho(:,:,2,tau), &
                              Time%model_time, rmask=Grd%umask(:,:,1),         &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_conv_rho_ud_t > 0) used = send_data (id_conv_rho_ud_t, Ext_mode%conv_rho_ud_t(:,:,tau), &
                                   Time%model_time, rmask=Grd%tmask(:,:,1),                      &
                                   is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_ps > 0)       used = send_data (id_ps, Ext_mode%ps(:,:),      &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_psx > 0)      used = send_data (id_psx, Ext_mode%grad_ps(:,:,1), &
                              Time%model_time, rmask=Grd%umask(:,:,1),    &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_psy > 0)      used = send_data (id_psy, Ext_mode%grad_ps(:,:,2), &
                              Time%model_time, rmask=Grd%umask(:,:,1),    &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_pb > 0)       used = send_data (id_pb, wrk1_2d(:,:),          &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_grad_anompbx > 0)  used = send_data (id_grad_anompbx, Ext_mode%grad_anompb(:,:,1), &
                            Time%model_time, rmask=Grd%umask(:,:,1),                        &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_grad_anompby > 0)  used = send_data (id_grad_anompby, Ext_mode%grad_anompb(:,:,2), &
                            Time%model_time, rmask=Grd%umask(:,:,1),                        &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! compute barotropic energetics  
  next_time = increment_time(Time%model_time, int(dtts), 0)
  if (need_data(id_ke_bt,next_time) .or. need_data(id_pe_bt,next_time)) then
    call barotropic_energy (Time, Ext_mode)
  endif 

  ! diagnose quasi-streamfunctions for vertically integrated transport
  psiu = 0.0
  psiv = 0.0
  if (need_data(id_psiu, next_time) .or. need_data(id_psiv, next_time)) then 
      call psi_compute(Time, Adv_vel, psiu, psiv)
      if(id_psiu > 0) used = send_data (id_psiu, psiu(:,:),           &
                             Time%model_time, rmask=Grd%tmask(:,:,1), &
                             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if(id_psiv > 0) used = send_data (id_psiv, psiv(:,:),           &
                             Time%model_time, rmask=Grd%tmask(:,:,1), &
                             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  if(have_obc) then
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,taup1),'M','n')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,taup1),'M','t')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,taup1),'Z','t')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,taup1),'Z','n')
  endif


  ! compute checksums for debugging 
  if (debug_this_module) then
     write(stdoutunit,*) ' '
     write(stdoutunit,*) 'From ocean_barotropic_mod: chksums at taup1'
     call write_timestamp(Time%model_time)
     call barotropic_chksum(Ext_mode, taup1)
  endif

end subroutine update_ocean_barotropic
! </SUBROUTINE> NAME="update_ocean_barotropic"


!#######################################################################
! <SUBROUTINE NAME="ocean_barotropic_forcing">
!
! <DESCRIPTION>
! Construct the vertically integrated forcing. This forcing is to be  
! held constant over the barotropic timesteps. At the time of calling
! this routine, accel has the contributions from explicit-time 
! forcing, except for the following:
!
! 1. Coriolis force is updated on the barotropic time steps when 
!    integrating the barotropic dynamics.  So it should not 
!    be included in forcing_bt.  
!
! 2. Contributions from smf and bmf are added to forcing_bt to allow 
!    for simpler handling of vertical friction implicitly in time. 
!
! 3. The accel array is already thickness and density weighted, so 
!    a vertical density weighted integral is here just a vertical sum.
!
! </DESCRIPTION>
!
subroutine ocean_barotropic_forcing(Time, Thickness, Velocity, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time 
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_velocity_type),      intent(inout) :: Velocity
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: i, j, k, n
  integer :: tau, taup1
  
  tau   = Time%tau
  taup1 = Time%taup1

  ! initialize the forcing to zero 
  do n=1,2
     do j=jsd,jed
        do i=isd,ied
           Ext_mode%forcing_bt(i,j,n) = 0.0
           wrk1_v2d(i,j,n)            = 0.0
        enddo
     enddo
  enddo

  ! vertical integral of nonlinear forcing
  if(.not. zero_nonlinear_forcing_bt) then 
      do n=1,2
         do k=1,nk 
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v2d(i,j,n) = wrk1_v2d(i,j,n) + Grd%umask(i,j,k)*Velocity%accel(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
  endif

  ! external forcing of barotropic velocity from smf and bmf. 
  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           Ext_mode%forcing_bt(i,j,n) = wrk1_v2d(i,j,n)  &
           + Grd%umask(i,j,1)*(Velocity%smf(i,j,n)-Velocity%bmf(i,j,n))
        enddo
     enddo
  enddo

  ! remove bottom drag from forcing_bt when it is instead
  ! applied on each barotropic step 
  if(barotropic_bmf) then 
      do n=1,2
         do j=jsc,jec
            do i=isc,iec
               Ext_mode%forcing_bt(i,j,n) = Ext_mode%forcing_bt(i,j,n) &
                                            + Grd%umask(i,j,1)*Velocity%bmf(i,j,n)
            enddo
         enddo
      enddo
  endif

  if (have_obc) call ocean_obc_adjust_forcing_bt(Ext_mode)

  call mpp_update_domains(Ext_mode%forcing_bt(:,:,1),Ext_mode%forcing_bt(:,:,2), &
                          Dom%domain2d, gridtype=BGRID_NE)



  ! horizontal friction due to just the time averaged barotropic velocity 
  friction_lap(:,:,:) = 0.0 
  friction_bih(:,:,:) = 0.0 
  if(udrho_lap) then 
      call lap_friction_barotropic(lap_ceu_back, lap_cnu_back, &
                                   Ext_mode%udrho(:,:,:,tau), friction_lap)
  endif
  if(udrho_bih) then 
      call bih_friction_barotropic(bih_ceu_back, bih_cnu_back, &
                                   Ext_mode%udrho(:,:,:,tau), friction_bih)
  endif
  do n=1,2
     do j=jsd,jed
        do i=isd,ied
           Ext_mode%forcing_bt(i,j,n) = Ext_mode%forcing_bt(i,j,n) &
                                      + friction_lap(i,j,n) + friction_bih(i,j,n)
           Velocity%lap_friction_bt(i,j,n) = friction_lap(i,j,n)
           Velocity%bih_friction_bt(i,j,n) = friction_bih(i,j,n)
        enddo
     enddo
  enddo


  ! no external forcing applied to barotropic velocity 
  if(zero_forcing_bt) then 
      do n=1,2
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%forcing_bt(i,j,n) = 0.0
            enddo
         enddo
      enddo
  endif

  if (id_nonlin_forcing_u_bt > 0) then 
       used = send_data (id_nonlin_forcing_u_bt, wrk1_v2d(:,:,1), &
                         Time%model_time, rmask=Grd%umask(:,:,1), &
                         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if (id_nonlin_forcing_v_bt > 0) then 
       used = send_data (id_nonlin_forcing_v_bt, wrk1_v2d(:,:,2), &
                         Time%model_time, rmask=Grd%umask(:,:,1), &
                         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if (id_forcing_u_bt > 0) then 
       used = send_data (id_forcing_u_bt, Ext_mode%forcing_bt(:,:,1), &
                         Time%model_time, rmask=Grd%umask(:,:,1),     &
                         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if (id_forcing_v_bt > 0) then 
       used = send_data (id_forcing_v_bt, Ext_mode%forcing_bt(:,:,2), &
                         Time%model_time, rmask=Grd%umask(:,:,1),     &
                         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 
  if (id_udrho_lap > 0) then
      used = send_data (id_udrho_lap, friction_lap(:,:,1), &
             Time%model_time, rmask=Grd%umask(:,:,1),      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_udrho_bih > 0) then
      used = send_data (id_udrho_bih, friction_bih(:,:,1), &
             Time%model_time, rmask=Grd%umask(:,:,1),      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_vdrho_lap > 0) then
      used = send_data (id_vdrho_lap, friction_lap(:,:,2), &
             Time%model_time, rmask=Grd%umask(:,:,1),      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_vdrho_bih > 0) then
      used = send_data (id_vdrho_bih, friction_bih(:,:,2), &
             Time%model_time, rmask=Grd%umask(:,:,1),      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif



end subroutine ocean_barotropic_forcing
! </SUBROUTINE> NAME="ocean_barotropic_forcing"


!#######################################################################
! <SUBROUTINE NAME="ocean_mass_forcing">
!
! <DESCRIPTION>
! Construct the vertically integrated mass source for 
! evolution of the free surface and/or bottom pressure.  
!
! Also construct the time tendency for the atmospheric pressure.
!
! </DESCRIPTION>
!
subroutine ocean_mass_forcing(Time, Thickness, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time 
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: i, j, k
  integer :: taum1, tau, taup1

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  ! vertically integrated mass source (kg/m^3)*(m/sec)  
  ! needed on halos in order for deta_dt and dpbot_dt to 
  ! also be on halos. 
  do j=jsd,jed
     do i=isd,ied
        Ext_mode%source(i,j) = 0.0
     enddo
  enddo
  do k=1,nk 
     do j=jsd,jed
        do i=isd,ied
           Ext_mode%source(i,j) = Ext_mode%source(i,j) + Thickness%mass_source(i,j,k)
        enddo
     enddo
  enddo
 
  ! time tendency for atmospheric pressure (Pa/sec)  
  ! values needed in halos in order for Thickness%rho_dzt_tendency
  ! to be defined well on halos. 
  do j=jsd,jed
     do i=isd,ied
        Ext_mode%dpatm_dt(i,j) = &
        Grd%tmask(i,j,1)*(Ext_mode%patm_t(i,j,taup1)-Ext_mode%patm_t(i,j,taum1))*dtimer
     enddo
  enddo

  ! atmospheric pressure on velocity cell point 
  Ext_mode%patm_u(:,:) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%patm_t(:,:,taup1))
 
  if (id_ext_mode_source > 0) used = send_data (id_ext_mode_source, Ext_mode%source(:,:), &
                                     Time%model_time, rmask=Grd%tmask(:,:,1), &
                                     is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

  if (id_patm_t > 0) used   = send_data (id_patm_t, Ext_mode%patm_t(:,:,taup1), &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_patm_u > 0) used   = send_data (id_patm_u, Ext_mode%patm_u(:,:), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_dpatm_dt > 0) used = send_data (id_dpatm_dt, Ext_mode%dpatm_dt(:,:), &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


end subroutine ocean_mass_forcing
! </SUBROUTINE> NAME="ocean_mass_forcing"


!#######################################################################
! <SUBROUTINE NAME="leap_frog_barotropic_depth">
!
! <DESCRIPTION>
! Integrate barotropic dynamics for "nts" timesteps using leapfrog.
! Assume a depth-like vertical coordinate, so solve for surface height. 
! </DESCRIPTION>
!
subroutine leap_frog_barotropic_depth (Time, Thickness, Ext_mode, patm, pme, river)

  type(ocean_time_type),          intent(in)    :: Time 
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in)    :: patm
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river

  type(time_type)                  :: time_bt 
  real, dimension(isd:ied,jsd:jed) :: steady_forcing

  integer :: fstaum1, fstau, fstaup1 
  integer :: itime, i, j, n
  integer :: tau, taum1, taup1
  integer :: day, sec 
  real    :: dayr    
  real    :: dtime

  if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      'Error from ocean_barotropic_mod (leap_frog_barotropic): module must be initialized')
  endif

  eta_eq_tidal = 0.0

  taum1 = Time%taum1
  tau   = Time%tau 
  taup1 = Time%taup1

  itime   = 1
  fstaum1 = mod(itime-1,3) + 1  
  fstau   = mod(itime+0,3) + 1

  if(tidal_forcing) then 
      time_bt = increment_time(Time%model_time, int(dteta),0) 
      call get_time(time_bt, sec, day)                                 
      dayr = day + sec/86400.0     
  endif

  do j=jsd,jed
     do i=isd,ied

        eta_t_bt(i,j,fstaum1)         = Ext_mode%eta_t_bar(i,j,tau)
        eta_t_bt(i,j,fstau)           = Ext_mode%eta_t_bar(i,j,tau)
        Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t(i,j,tau)     ! for accumulation of eta_t_bt

        ! Compute that part of the eta_t_bt forcing which is constant
        ! over barotropic integration. units (kg/m^3)*(m/s).
        ! Remove Ext_mode%eta_smooth from Ext_mode%source, 
        ! since Ext_mode%eta_smooth is only to be used for 
        ! smoothing eta_t, not for smoothing eta_t_bt.
        steady_forcing(i,j) = Grd%tmask(i,j,1)*( pme(i,j) + river(i,j) &
                             +Ext_mode%source(i,j) - Ext_mode%eta_smooth(i,j))

     enddo
  enddo

  do n=1,2
     do j=jsd,jed
        do i=isd,ied
           udrho_bt(i,j,n,fstaum1)     = Ext_mode%udrho(i,j,n,tau)
           udrho_bt(i,j,n,fstau)       = Ext_mode%udrho(i,j,n,tau)
           Ext_mode%udrho(i,j,n,taup1) = Ext_mode%udrho(i,j,n,tau)  ! to accumulate udrho_bt  
        enddo
     enddo
  enddo

  ! step over the nts barotropic time steps 
  do itime=1,nts

     if (tidal_forcing) then
         dayr = dayr + dtbt/86400.0 
         call get_tidal_forcing(Time, dayr)
     endif

     ! set leap-frog time level indices
     fstaum1 = mod(itime-1,3) + 1 
     fstau   = mod(itime+0,3) + 1
     fstaup1 = mod(itime+1,3) + 1

     if (itime == 1) then
         dtime = 0.5*twodt ! first timestep is euler forward
     else
         dtime = twodt ! remaining timesteps are leapfrog
     endif

     ! update free surface height on barotropic time step
     tmp(:,:) = -DIV_UD(udrho_bt(:,:,:,fstau)) * Grd%tmask(:,:,1)
     ! set DIV_UD to zero at open boundaries if needed
     if(have_obc) call ocean_obc_adjust_divud(tmp)        

     do j=jsd,jed
        do i=isd,ied      
           eta_t_bt(i,j,fstaup1) = eta_t_bt(i,j,fstaum1) &
           + dtime*rho0r*(tmp(i,j) + steady_forcing(i,j))
        enddo
     enddo

     ! open boundary condition at barotropic time step
     ! This message passing is needed only, if the boundary scheme needs 
     ! values from two domains. A check should be introduced, which 
     ! switches off this call, if not needed, since message passing is  
     ! done here at high frequency.
     if(have_obc) call ocean_obc_barotropic(eta_t_bt, fstaum1, fstau, fstaup1, dtime)

     ! smooth eta (relative to eta_t_bt=0)
     if(smooth_eta_t_bt_laplacian) then
         eta_t_bt(:,:,fstaup1) = eta_t_bt(:,:,fstaup1) &
           + dtime*LAP_T(eta_t_bt(:,:,fstaum1),smooth_lap(:,:))
     endif
     if(smooth_eta_t_bt_biharmonic) then
         tmp(:,:) = -LAP_T(eta_t_bt(:,:,fstaum1), smooth_bih(:,:))
         call mpp_update_domains (tmp(:,:), Dom%domain2d)
         eta_t_bt(:,:,fstaup1) = eta_t_bt(:,:,fstaup1) &
         + dtime*LAP_T(tmp(:,:), smooth_bih(:,:))
     endif

     call mpp_update_domains (eta_t_bt(:,:,fstaup1), Dom%domain2d)   

     if (tidal_forcing) then
         do j=jsd,jed
            do i=isd, ied
               Ext_mode%ps(i,j) = rho_g(i,j) &
               *(alphat*eta_t_bt(i,j,fstau)-eta_eq_tidal(i,j)) + patm(i,j)
            enddo
         enddo
     else 
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%ps(i,j)=rho_g(i,j)*eta_t_bt(i,j,fstau) + patm(i,j)
            enddo
         enddo
     endif
     if(geoid_forcing) then 
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%ps(i,j) = Ext_mode%ps(i,j) - rho_g(i,j)*eta_geoid(i,j)
            enddo
         enddo
     endif 

     ! compute the pressure force per area (N/m^2) 
     Ext_mode%grad_ps(:,:,:) = GRAD_BAROTROPIC_P(Ext_mode%ps(:,:))
     do j=jsd,jed
        do i=isd,ied
           Ext_mode%press_force(i,j,1) = -Thickness%thicku(i,j,tau)*Ext_mode%grad_ps(i,j,1)
           Ext_mode%press_force(i,j,2) = -Thickness%thicku(i,j,tau)*Ext_mode%grad_ps(i,j,2)
        enddo
     enddo
     if (Grd%tripolar) then 
         call mpp_update_domains (Ext_mode%press_force(:,:,1),Ext_mode%press_force(:,:,2), &
                                  Dom%domain2d, gridtype=BGRID_NE)
     endif


     ! apply bottom drag to the vertically integrated momentum 
     wrk1_v2d(:,:,:) = 0.0
     if(barotropic_bmf) then 
         do n=1,2
            do j=jsd,jed
               do i=isd,ied
                  wrk1_v2d(i,j,n) = -rho0r*barotropic_bmf_cdbot &
                  *thicku_r(i,j)*udrho_bt(i,j,n,fstau)*thicku_r(i,j)*udrho_bt(i,j,n,fstau)
               enddo
            enddo
         enddo
     endif

     ! horizontal friction
     friction_lap(:,:,:) = 0.0 
     friction_bih(:,:,:) = 0.0 
     if(udrho_bt_lap .or. udrho_bt_bih) then 
         call mpp_update_domains (udrho_bt(:,:,1,fstau),udrho_bt(:,:,2,fstau), &
                                  Dom%domain2d, gridtype=BGRID_NE)
     endif 
     if(udrho_bt_lap) then 
         call lap_friction_barotropic(lap_ceu_back, lap_cnu_back, &
                                      udrho_bt(:,:,:,fstau), friction_lap)
     endif 
     if(udrho_bt_bih) then 
         call bih_friction_barotropic(bih_ceu_back, bih_cnu_back, &
                                      udrho_bt(:,:,:,fstau), friction_bih)
     endif

     ! update vertically integrated horizontal velocity using leap-frog time step 
     do j=jsd,jed
        do i=isd,ied
           udrho_bt(i,j,1,fstaup1) = udrho_bt(i,j,1,fstaum1) &
           + dtime*( Grd%f(i,j)*udrho_bt(i,j,2,fstau)        &
           + Ext_mode%press_force(i,j,1)                     &
           + Ext_mode%forcing_bt(i,j,1)                      &
           + wrk1_v2d(i,j,1)                                 &
           + friction_lap(i,j,1) + friction_bih(i,j,1) )

           udrho_bt(i,j,2,fstaup1) = udrho_bt(i,j,2,fstaum1) &
           + dtime*(-Grd%f(i,j)*udrho_bt(i,j,1,fstau)        &
           + Ext_mode%press_force(i,j,2)                     &
           + Ext_mode%forcing_bt(i,j,2)                      &
           + wrk1_v2d(i,j,2)                                 &
           + friction_lap(i,j,2) + friction_bih(i,j,2) )
        enddo
     enddo

     ! update barotropic velocity on the global halo 
     ! gradient accross boundary = 0 of cross boundary flux 
     ! minimize +/- structure for along boundary flux
     if(have_obc) then
         call ocean_obc_update_boundary(udrho_bt(:,:,1,fstaup1),'M','s')
         call ocean_obc_update_boundary(udrho_bt(:,:,2,fstaup1),'M','i')
         call ocean_obc_update_boundary(udrho_bt(:,:,1,fstaup1),'Z','i')
         call ocean_obc_update_boundary(udrho_bt(:,:,2,fstaup1),'Z','s')
     endif

     ! suppress leap-frog splitting mode using Robert-Asselin time filter 
     if(robert_asselin_bt > 0.0) then 
         do j=jsd,jed
            do i=isd,ied
               eta_t_bt(i,j,fstau)   = eta_t_bt(i,j,fstau) + &
               robert_asselin_bt*(0.5*(eta_t_bt(i,j,fstaup1)+eta_t_bt(i,j,fstaum1))-eta_t_bt(i,j,fstau))
               udrho_bt(i,j,:,fstau) = udrho_bt(i,j,:,fstau)  + &
               robert_asselin_bt*(0.5*(udrho_bt(i,j,:,fstaup1)+udrho_bt(i,j,:,fstaum1))-udrho_bt(i,j,:,fstau))
            enddo
         enddo
     endif 

     ! accumulate to build time average   
     do j=jsd,jed
        do i=isd,ied
           Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t_bar(i,j,taup1) + eta_t_bt(i,j,fstaup1)
           Ext_mode%udrho(i,j,1,taup1)   = Ext_mode%udrho(i,j,1,taup1)   + udrho_bt(i,j,1,fstaup1)
           Ext_mode%udrho(i,j,2,taup1)   = Ext_mode%udrho(i,j,2,taup1)   + udrho_bt(i,j,2,fstaup1)
        enddo
     enddo

     ! take a sample from the centre of the barotropic loop 
     if(itime==int(nts/2)) then 
         if (id_eta_t_bt > 0) then
             used = send_data (id_eta_t_bt, eta_t_bt(:,:,tau), &
                    Time%model_time, rmask=Grd%tmask(:,:,1),   &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_udrho_bt_lap > 0) then
             used = send_data (id_udrho_bt_lap, friction_lap(:,:,1), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_udrho_bt_bih > 0) then
             used = send_data (id_udrho_bt_bih, friction_bih(:,:,1), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_vdrho_bt_lap > 0) then
             used = send_data (id_vdrho_bt_lap, friction_lap(:,:,2), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_vdrho_bt_bih > 0) then
             used = send_data (id_vdrho_bt_bih, friction_bih(:,:,2), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
     endif


  enddo ! end of barotropic integration itime=1,nts 

  if(have_obc) call ocean_obc_update_boundary(Ext_mode%eta_t_bar(:,:,taup1),'T')

end subroutine leap_frog_barotropic_depth
! </SUBROUTINE> NAME="leap_frog_barotropic_depth"



!#######################################################################
! <SUBROUTINE NAME="pred_corr_barotropic_depth">
!
! <DESCRIPTION>
! Integrate barotropic dynamics for "nts" timesteps using 
! predictor-corrector. Assume depth-like vertical coordinate
! so solve for surface height.   
!
! This scheme is more stable than the leap_frog since it can run with 
! longer time steps to resolve external mode gravity waves.  It also  
! provides some extra smoothing when pred_corr_gamma > 0 and so the options 
! smooth_eta_t_bt_laplacian and smooth_eta_t_bt_biharmonic may not be
! needed.  
!
! Note that OBC has not been tested for use with predictor-corrector.
! 
! </DESCRIPTION>
!
subroutine pred_corr_barotropic_depth (Time, Thickness, Ext_mode, patm, pme, river)

  type(ocean_time_type),          intent(in)    :: Time 
  type(ocean_thickness_type),     intent(in)    :: Thickness 
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in)    :: patm
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river

  type(time_type)                  :: time_bt 
  real, dimension(isd:ied,jsd:jed) :: steady_forcing

  integer :: itime, i, j, n
  integer :: tau, taup1
  integer :: fstau, fstaup1 
  integer :: day, sec
  real    :: dayr 
  real    :: urhod_tmp1, urhod_tmp2

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     'Error from ocean_barotropic_mod (pred_corr_barotropic_depth): module must be initialized')
  endif
 
  eta_eq_tidal = 0.0

  tau   = Time%tau 
  taup1 = Time%taup1

  itime   = 1
  fstau   = mod(itime+0,2) + 1
  fstaup1 = mod(itime+1,2) + 1

  if(tidal_forcing) then 
      time_bt = increment_time(Time%model_time, int(dteta),0) 
      call get_time(time_bt, sec, day)                                 
      dayr = day + sec/86400.0     
  endif

  ! initialize for accumulating to take time average. eta_t_bar includes
  ! endpoints, whereas udrho does not. hence the different initialization.
  if(initsum_with_bar) then 
    do j=jsd,jed
       do i=isd,ied
    	  Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t_bar(i,j,tau) 
    	  Ext_mode%udrho(i,j,1,taup1)	= 0.0
    	  Ext_mode%udrho(i,j,2,taup1)	= 0.0
       enddo
    enddo
  else
    do j=jsd,jed
       do i=isd,ied
    	  Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t(i,j,tau) 
    	  Ext_mode%udrho(i,j,1,taup1)	= 0.0
    	  Ext_mode%udrho(i,j,2,taup1)	= 0.0
       enddo
    enddo
  endif
  do j=jsd,jed
     do i=isd,ied

        eta_t_bt(i,j,fstau)   = Ext_mode%eta_t_bar(i,j,tau)
        eta_t_bt(i,j,fstaup1) = Ext_mode%eta_t_bar(i,j,tau)
        udrho_bt(i,j,1,fstau) = Ext_mode%udrho(i,j,1,tau)
        udrho_bt(i,j,2,fstau) = Ext_mode%udrho(i,j,2,tau)

        ! Compute that part of the eta_t_bt forcing which is constant
        ! over barotropic integration. units (kg/m^3)*(m/s).
        ! Remove Ext_mode%eta_smooth from Ext_mode%source, 
        ! since Ext_mode%eta_smooth is only to be used for 
        ! smoothing eta_t, not for smoothing eta_t_bt.
        steady_forcing(i,j) = Grd%tmask(i,j,1)*( pme(i,j) + river(i,j) &
                             +Ext_mode%source(i,j) - Ext_mode%eta_smooth(i,j))

     enddo
  enddo

  ! step over the nts barotropic time steps 
  do itime=1,nts

     ! set predictor-corrector time level indices
     fstau   = mod(itime+0,2) + 1
     fstaup1 = mod(itime+1,2) + 1

     ! predictor step for eta
     if(pred_corr_gamma > 0.0) then 

         tmp(:,:) = -DIV_UD(udrho_bt(:,:,:,fstau)) * Grd%tmask(:,:,1)          
         if(have_obc) call ocean_obc_adjust_divud(tmp)        
         do j=jsd,jed
            do i=isd,ied
               eta_t_bt(i,j,fstaup1) = eta_t_bt(i,j,fstau) &
               + dtbt_gamma_rho0r*(tmp(i,j) + steady_forcing(i,j)) 
            enddo
         enddo
         if(have_obc) call ocean_obc_barotropic(eta_t_bt, fstau, fstau, fstaup1, dtbt_gamma)

     endif

     call mpp_update_domains (eta_t_bt(:,:,fstaup1), Dom%domain2d)      

     if (tidal_forcing) then
         dayr = dayr + dtbt/86400.0 
         call get_tidal_forcing(Time, dayr)
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%ps(i,j) = rho_g(i,j) &
               *(alphat*eta_t_bt(i,j,fstaup1)-eta_eq_tidal(i,j)) + patm(i,j)
            enddo
         enddo
     else 
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%ps(i,j) = rho_g(i,j)*eta_t_bt(i,j,fstaup1) + patm(i,j)
            enddo
         enddo
     endif
     if(geoid_forcing) then 
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%ps(i,j) = Ext_mode%ps(i,j) - rho_g(i,j)*eta_geoid(i,j)
            enddo
         enddo
     endif 


     ! compute the pressure force per area (N/m^2) 
     Ext_mode%grad_ps(:,:,:) = GRAD_BAROTROPIC_P(Ext_mode%ps(:,:))
     do j=jsd,jed
        do i=isd,ied
           Ext_mode%press_force(i,j,1) = -Thickness%thicku(i,j,tau)*Ext_mode%grad_ps(i,j,1)
           Ext_mode%press_force(i,j,2) = -Thickness%thicku(i,j,tau)*Ext_mode%grad_ps(i,j,2)
        enddo
     enddo
     if (Grd%tripolar) then 
         call mpp_update_domains (Ext_mode%press_force(:,:,1),Ext_mode%press_force(:,:,2), &
                                  Dom%domain2d, gridtype=BGRID_NE)
     endif

     ! apply bottom drag to the vertically integrated momentum 
     wrk1_v2d(:,:,:) = 0.0

     if(have_obc) call ocean_obc_damp_newton(udrho_bt(:,:,:,fstau),wrk1_v2d)
 
     if(barotropic_bmf) then 
         do n=1,2
            do j=jsd,jed
               do i=isd,ied
                  wrk1_v2d(i,j,n) = wrk1_v2d(i,j,n) -rho0r*barotropic_bmf_cdbot &
                  *thicku_r(i,j)*udrho_bt(i,j,n,fstau)*thicku_r(i,j)*udrho_bt(i,j,n,fstau)
               enddo
            enddo
         enddo
     endif

     ! horizontal friction
     friction_lap(:,:,:) = 0.0 
     friction_bih(:,:,:) = 0.0 

     if(udrho_bt_lap .or. udrho_bt_bih) then 

        call mpp_update_domains (udrho_bt(:,:,1,fstau),udrho_bt(:,:,2,fstau), &
                                  Dom%domain2d, gridtype=BGRID_NE)
        if(udrho_bt_lap) then 
            call lap_friction_barotropic(lap_ceu_back, lap_cnu_back, &
                                      udrho_bt(:,:,:,fstau), friction_lap)
        endif 
        if(udrho_bt_bih) then 
            call bih_friction_barotropic(bih_ceu_back, bih_cnu_back, &
                                      udrho_bt(:,:,:,fstau), friction_bih)
        endif

        ! update with explicit time pieces and implicit Coriolis 
        do j=jsd,jed
           do i=isd,ied

              urhod_tmp1 =  udrho_bt(i,j,1,fstau)	    &
              + dtbt*( 0.5*Grd%f(i,j)*udrho_bt(i,j,2,fstau) &
              + Ext_mode%press_force(i,j,1)		    &
              + Ext_mode%forcing_bt(i,j,1)		    &
              + wrk1_v2d(i,j,1) 			    &
              + friction_lap(i,j,1) + friction_bih(i,j,1) )

              urhod_tmp2 =  udrho_bt(i,j,2,fstau)	    &
              + dtbt*(-0.5*Grd%f(i,j)*udrho_bt(i,j,1,fstau) &
              + Ext_mode%press_force(i,j,2)		    &
              + Ext_mode%forcing_bt(i,j,2)		    &
              + wrk1_v2d(i,j,2) 			    &
              + friction_lap(i,j,2) + friction_bih(i,j,2) )

              udrho_bt(i,j,1,fstaup1) = cori2(i,j)*(urhod_tmp1 + cori1(i,j)*urhod_tmp2)
              udrho_bt(i,j,2,fstaup1) = cori2(i,j)*(urhod_tmp2 - cori1(i,j)*urhod_tmp1)
           enddo
        enddo

     else

        do j=jsd,jed
           do i=isd,ied

              urhod_tmp1 =  udrho_bt(i,j,1,fstau)	    &
              + dtbt*( 0.5*Grd%f(i,j)*udrho_bt(i,j,2,fstau) &
              + Ext_mode%press_force(i,j,1)		    &
              + Ext_mode%forcing_bt(i,j,1)		    &
              + wrk1_v2d(i,j,1) )

              urhod_tmp2 =  udrho_bt(i,j,2,fstau)	    &
              + dtbt*(-0.5*Grd%f(i,j)*udrho_bt(i,j,1,fstau) &
              + Ext_mode%press_force(i,j,2)		    &
              + Ext_mode%forcing_bt(i,j,2)		    &
              + wrk1_v2d(i,j,2) )

              udrho_bt(i,j,1,fstaup1) = cori2(i,j)*(urhod_tmp1 + cori1(i,j)*urhod_tmp2)
              udrho_bt(i,j,2,fstaup1) = cori2(i,j)*(urhod_tmp2 - cori1(i,j)*urhod_tmp1)

           enddo
        enddo

     endif ! endif for if(udrho_bt_lap .or. udrho_bt_bih) then 


     ! update barotropic velocity on the global halo 
     ! gradient accross boundary = 0 of cross boundary flux 
     ! minimize +/- structure for along boundary flux
     if (have_obc) call ocean_obc_ud(eta_t_bt(:,:,fstaup1),udrho_bt(:,:,:,fstaup1))        


     ! corrector step for eta 
     tmp(:,:) = -DIV_UD(udrho_bt(:,:,:,fstaup1)) * Grd%tmask(:,:,1)
     if(have_obc) call ocean_obc_adjust_divud(tmp)        
     do j=jsd,jed
        do i=isd,ied          
           eta_t_bt(i,j,fstaup1) = eta_t_bt(i,j,fstau) &
           + dtbt_rho0r*(tmp(i,j) + steady_forcing(i,j))
        enddo
     enddo
     if(have_obc) call ocean_obc_barotropic(eta_t_bt, fstau, fstau, fstaup1, dtbt)

     ! smooth eta_t_bt (relative to eta_t_bt=0).
     ! time consuming due to mpp_update_domain calls
     if(smooth_eta_t_bt_laplacian) then
         eta_t_bt(:,:,fstaup1) = eta_t_bt(:,:,fstaup1) &
         + dtbt*LAP_T(eta_t_bt(:,:,fstau), smooth_lap(:,:))
     endif
     if(smooth_eta_t_bt_biharmonic) then
         tmp(:,:) = -LAP_T(eta_t_bt(:,:,fstau), smooth_bih(:,:))
         call mpp_update_domains (tmp(:,:), Dom%domain2d)
         eta_t_bt(:,:,fstaup1) = eta_t_bt(:,:,fstaup1) &
         + dtbt*LAP_T(tmp(:,:),smooth_bih(:,:))
     endif
     if (update_domains_for_obc.or.smooth_eta_t_bt_laplacian.or.smooth_eta_t_bt_biharmonic) &
	 call mpp_update_domains (eta_t_bt(:,:,fstaup1), Dom%domain2d)

     ! accumulate for time average 
     if(barotropic_time_stepping_mom4p0) then 
         do j=jsd,jed
            do i=isd,ied   
               Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t_bar(i,j,taup1) &
                    + eta_t_bt(i,j,fstaup1)
               Ext_mode%udrho(i,j,1,taup1)   = Ext_mode%udrho(i,j,1,taup1)   &
                    + udrho_bt(i,j,1,fstaup1)
               Ext_mode%udrho(i,j,2,taup1)   = Ext_mode%udrho(i,j,2,taup1)   &
                    + udrho_bt(i,j,2,fstaup1)
            enddo
         enddo
     else 
         do j=jsd,jed
            do i=isd,ied   
               Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t_bar(i,j,taup1) &
                    + eta_t_bt(i,j,fstaup1)
               Ext_mode%udrho(i,j,1,taup1)   = Ext_mode%udrho(i,j,1,taup1)   &
                    + float(nts-itime+1)*udrho_bt(i,j,1,fstaup1)
               Ext_mode%udrho(i,j,2,taup1)   = Ext_mode%udrho(i,j,2,taup1)   &
                    + float(nts-itime+1)*udrho_bt(i,j,2,fstaup1)
            enddo
         enddo
     endif


     ! take a sample from the centre of the barotropic loop 
     if(itime==int(nts/2)) then 
         if (id_eta_t_bt > 0) then
             used = send_data (id_eta_t_bt, eta_t_bt(:,:,tau), &
                    Time%model_time, rmask=Grd%tmask(:,:,1),   &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_udrho_bt_lap > 0) then
             used = send_data (id_udrho_bt_lap, friction_lap(:,:,1), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_udrho_bt_bih > 0) then
             used = send_data (id_udrho_bt_bih, friction_bih(:,:,1), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_vdrho_bt_lap > 0) then
             used = send_data (id_vdrho_bt_lap, friction_lap(:,:,2), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_vdrho_bt_bih > 0) then
             used = send_data (id_vdrho_bt_bih, friction_bih(:,:,2), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
     endif


  enddo ! end of barotropic integration itime=1,nts 

  if(have_obc) call ocean_obc_update_boundary(Ext_mode%eta_t_bar(:,:,taup1),'T')

end subroutine pred_corr_barotropic_depth
! </SUBROUTINE> NAME="pred_corr_barotropic_depth"



!#######################################################################
! <SUBROUTINE NAME="pred_corr_barotropic_press">
!
! <DESCRIPTION>
! Integrate barotropic dynamics for "nts" timesteps using 
! predictor-corrector. Assume pressure-like vertical coordinate
! so solve here for the bottom pressure. 
!
! This scheme provides some smoothing of small spatial scales
! when pred_corr_gamma > 0.  
!
! NOTE: the pressure gradient force is based on gradients of 
! (pbot_t_bt - rho0*grav*ht) rather than gradients of pbot_t_bt.
! This approach aims to improve accuracy of the pressure force. 
!
! </DESCRIPTION>
!
subroutine pred_corr_barotropic_press (Time, Thickness, Ext_mode, pme, river)

  type(ocean_time_type),          intent(in)    :: Time 
  type(ocean_thickness_type),     intent(in)    :: Thickness 
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river

  type(time_type)                    :: time_bt 
  real, dimension(isd:ied,jsd:jed)   :: steady_forcing
  real, dimension(isd:ied,jsd:jed,3) :: anompb_obc

  integer :: itime, i, j, n
  integer :: tau, taup1
  integer :: fstau, fstaup1 
  integer :: day, sec
  real    :: dayr 
  real    :: urhod_tmp1, urhod_tmp2

  if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      'Error from ocean_barotropic_mod (pred_corr_barotropic_press): module not initialized')
  endif
 
  tau   = Time%tau 
  taup1 = Time%taup1

  eta_eq_tidal = 0.0 

  itime   = 1
  fstau   = mod(itime+0,2) + 1
  fstaup1 = mod(itime+1,2) + 1

  if(tidal_forcing) then 
      time_bt = increment_time(Time%model_time, int(dteta),0) 
      call get_time(time_bt, sec, day)                                 
      dayr = day + sec/86400.0     
  endif

  anompb_obc(:,:,:) = 0.0

  ! initialize for accumulating to take time average. anompb_bar includes
  ! endpoints, whereas udrho does not. hence the different initialization.
  if(initsum_with_bar) then 
     do j=jsd,jed
       do i=isd,ied
          Ext_mode%anompb_bar(i,j,taup1) = Ext_mode%anompb_bar(i,j,tau)
          Ext_mode%udrho(i,j,1,taup1)    = 0.0
          Ext_mode%udrho(i,j,2,taup1)    = 0.0
      enddo
     enddo
  else
     do j=jsd,jed
       do i=isd,ied
          Ext_mode%anompb_bar(i,j,taup1) = Ext_mode%anompb(i,j,tau)
          Ext_mode%udrho(i,j,1,taup1)    = 0.0
          Ext_mode%udrho(i,j,2,taup1)    = 0.0
       enddo
     enddo
  endif
  do j=jsd,jed
     do i=isd,ied
 
        anompb_bt(i,j,fstau)   = Ext_mode%anompb_bar(i,j,tau)
        anompb_bt(i,j,fstaup1) = Ext_mode%anompb_bar(i,j,tau)
        udrho_bt(i,j,1,fstau)  = Ext_mode%udrho(i,j,1,tau)
        udrho_bt(i,j,2,fstau)  = Ext_mode%udrho(i,j,2,tau)

        ! compute that part of the anompb_bt forcing which is constant 
        ! over barotropic integration. 
        ! units (m/s^2)*(kg/m^3)*(m/s) = kg/(m*s^3) = Pa/s
        ! Remove Ext_mode%pbot_smooth from Ext_mode%source, 
        ! since Ext_mode%pbot_smooth is only to be used for 
        ! smoothing pbot, not for smoothing anompb_bt.
        steady_forcing(i,j) = Grd%tmask(i,j,1)                                            &
        *(grav*(pme(i,j) + river(i,j) + Ext_mode%source(i,j) - Ext_mode%pbot_smooth(i,j)) &
        + Ext_mode%dpatm_dt(i,j))
 
     enddo
  enddo

  ! step over the nts barotropic time steps 
  do itime=1,nts

     ! set predictor-corrector time level indices
     fstau   = mod(itime+0,2) + 1
     fstaup1 = mod(itime+1,2) + 1

     ! predictor step for eta
     if(pred_corr_gamma > 0.0) then 

         tmp(:,:) = -grav*DIV_UD(udrho_bt(:,:,:,fstau))          
         do j=jsd,jed
            do i=isd,ied
               anompb_bt(i,j,fstaup1) = anompb_bt(i,j,fstau) &
               + dtbt_gamma*(tmp(i,j) + steady_forcing(i,j)) 
            enddo
         enddo

         if(have_obc) then 
         ! open boundary condition at barotropic time step
         ! This message passing is needed only, if the boundary scheme needs 
         ! values from two domains. A check should be introduced, which 
         ! switches off this call, if not needed, since message passing is  
         ! done here at high frequency.
             do j=jsd,jed
                do i=isd,ied
                   anompb_obc(i,j,fstaup1) = anompb_bt(i,j,fstaup1)
                enddo
             enddo
             call ocean_obc_barotropic(anompb_obc, fstau, fstau, fstaup1, pred_corr_gamma*dtbt)
         endif

     endif

     call mpp_update_domains (anompb_bt(:,:,fstaup1), Dom%domain2d)      

     wrk1_2d(:,:) = 0.0
     if (tidal_forcing) then
         dayr = dayr + dtbt/86400.0 
         call get_tidal_forcing(Time, dayr)
         do j=jsd,jed
            do i=isd,ied
               wrk1_2d(i,j) = Grd%tmask(i,j,1)*(anompb_bt(i,j,fstaup1)                     &
                              +rho0grav*(alphat*Ext_mode%eta_t(i,j,tau)-eta_eq_tidal(i,j)) &
                               )
            enddo
         enddo
     else
         do j=jsd,jed
            do i=isd,ied
               wrk1_2d(i,j) = Grd%tmask(i,j,1)*anompb_bt(i,j,fstaup1)
            enddo
         enddo
     endif
     if(geoid_forcing) then 
         do j=jsd,jed
            do i=isd,ied
               wrk1_2d(i,j) = wrk1_2d(i,j) - Grd%tmask(i,j,1)*rho0grav*eta_geoid(i,j)
            enddo
         enddo
     endif 


     ! compute the pressure force per area (N/m^2) 
     Ext_mode%grad_anompb(:,:,:) = GRAD_BAROTROPIC_P(wrk1_2d(:,:))
     do j=jsd,jed
        do i=isd,ied
           Ext_mode%press_force(i,j,1) = -rho0r*Thickness%mass_u(i,j,tau)*Ext_mode%grad_anompb(i,j,1)
           Ext_mode%press_force(i,j,2) = -rho0r*Thickness%mass_u(i,j,tau)*Ext_mode%grad_anompb(i,j,2)
        enddo
     enddo
     if (Grd%tripolar) then 
         call mpp_update_domains (Ext_mode%press_force(:,:,1),Ext_mode%press_force(:,:,2), &
                                  Dom%domain2d, gridtype=BGRID_NE)
     endif

     ! apply bottom drag to the vertically integrated momentum 
     wrk1_v2d(:,:,:) = 0.0
     if(barotropic_bmf) then 
         do n=1,2
            do j=jsd,jed
               do i=isd,ied
                  wrk1_v2d(i,j,n) = -rho0r*barotropic_bmf_cdbot &
                  *thicku_r(i,j)*udrho_bt(i,j,n,fstau)*thicku_r(i,j)*udrho_bt(i,j,n,fstau)
               enddo
            enddo
         enddo
     endif

     ! horizontal friction
     friction_lap(:,:,:) = 0.0 
     friction_bih(:,:,:) = 0.0 

     if(udrho_bt_lap .or. udrho_bt_bih) then 

         call mpp_update_domains (udrho_bt(:,:,1,fstau),udrho_bt(:,:,2,fstau), &
              Dom%domain2d, gridtype=BGRID_NE)
         if(udrho_bt_lap) then 
             call lap_friction_barotropic(lap_ceu_back, lap_cnu_back, &
                  udrho_bt(:,:,:,fstau), friction_lap)
         endif
         if(udrho_bt_bih) then 
             call bih_friction_barotropic(bih_ceu_back, bih_cnu_back, &
                  udrho_bt(:,:,:,fstau), friction_bih)
         endif


         ! update with explicit time pieces and implicit Coriolis 
         do j=jsd,jed
            do i=isd,ied

               urhod_tmp1 = udrho_bt(i,j,1,fstau)                  &
                    + dtbt*( 0.5*Grd%f(i,j)*udrho_bt(i,j,2,fstau)  &
                    + Ext_mode%press_force(i,j,1)                  &
                    + Ext_mode%forcing_bt(i,j,1)                   &
                    + wrk1_v2d(i,j,1)                              &
                    + friction_lap(i,j,1) + friction_bih(i,j,1) )

               urhod_tmp2 = udrho_bt(i,j,2,fstau)                  &
                    + dtbt*(-0.5*Grd%f(i,j)*udrho_bt(i,j,1,fstau)  &
                    + Ext_mode%press_force(i,j,2)                  & 
                    + Ext_mode%forcing_bt(i,j,2)                   &
                    + wrk1_v2d(i,j,2)                              &
                    + friction_lap(i,j,2) + friction_bih(i,j,2) )

               udrho_bt(i,j,1,fstaup1) = cori2(i,j)*(urhod_tmp1 + cori1(i,j)*urhod_tmp2)
               udrho_bt(i,j,2,fstaup1) = cori2(i,j)*(urhod_tmp2 - cori1(i,j)*urhod_tmp1)

            enddo
         enddo

     else 

         do j=jsd,jed
            do i=isd,ied

               urhod_tmp1 = udrho_bt(i,j,1,fstau)                 &
                    + dtbt*( 0.5*Grd%f(i,j)*udrho_bt(i,j,2,fstau) &
                    + Ext_mode%press_force(i,j,1)                 &
                    + Ext_mode%forcing_bt(i,j,1)                  &
                    + wrk1_v2d(i,j,1) )                                     

               urhod_tmp2 = udrho_bt(i,j,2,fstau)                  &
                    + dtbt*(-0.5*Grd%f(i,j)*udrho_bt(i,j,1,fstau)  &
                    + Ext_mode%press_force(i,j,2)                  & 
                    + Ext_mode%forcing_bt(i,j,2)                   &
                    + wrk1_v2d(i,j,2) )                                    

               udrho_bt(i,j,1,fstaup1) = cori2(i,j)*(urhod_tmp1 + cori1(i,j)*urhod_tmp2)
               udrho_bt(i,j,2,fstaup1) = cori2(i,j)*(urhod_tmp2 - cori1(i,j)*urhod_tmp1)

            enddo
         enddo

     endif ! endif for if(udrho_bt_lap .or. udrho_bt_bih) then 


     if(have_obc) then
     ! update barotropic velocity on the global halo 
     ! gradient accross boundary = 0 of cross boundary flux 
     ! minimize +/- structure for along boundary flux
         call ocean_obc_update_boundary(udrho_bt(:,:,1,fstaup1),'M','n')
         call ocean_obc_update_boundary(udrho_bt(:,:,2,fstaup1),'M','t')
         call ocean_obc_update_boundary(udrho_bt(:,:,1,fstaup1),'Z','t')
         call ocean_obc_update_boundary(udrho_bt(:,:,2,fstaup1),'Z','n')
     endif

     ! corrector step for anompb_bt
     tmp(:,:) = -grav*DIV_UD(udrho_bt(:,:,:,fstaup1))
     do j=jsd,jed
        do i=isd,ied          
           anompb_bt(i,j,fstaup1) = anompb_bt(i,j,fstau) + dtbt*(tmp(i,j) + steady_forcing(i,j))
        enddo
     enddo

     ! open boundary condition at barotropic time step
     ! This message passing is needed only, if the boundary scheme needs 
     ! values from two domains. A check should be introduced, which 
     ! switches off this call, if not needed, since message passing is  
     ! done here at high frequency.
     if(have_obc) then 
         do j=jsd,jed
            do i=isd,ied
               anompb_obc(i,j,fstaup1) = anompb_bt(i,j,fstaup1) 
            enddo
         enddo
         call ocean_obc_barotropic(anompb_obc, fstau, fstau, fstaup1, dtbt)
     endif

     ! smooth anompb_bt.
     ! time consuming due to mpp_update_domain calls
     if(smooth_anompb_bt_laplacian) then
         tmp(:,:) = anompb_bt(:,:,fstau)
         call mpp_update_domains (tmp(:,:), Dom%domain2d)
         anompb_bt(:,:,fstaup1) = anompb_bt(:,:,fstaup1) &
         + dtbt*LAP_T(tmp(:,:), smooth_lap(:,:))
     endif
     if(smooth_anompb_bt_biharmonic) then
         tmp(:,:) = anompb_bt(:,:,fstau) 
         call mpp_update_domains (tmp(:,:), Dom%domain2d)
         tmp(:,:) = -LAP_T(tmp(:,:), smooth_bih(:,:))
         call mpp_update_domains (tmp(:,:), Dom%domain2d)
         anompb_bt(:,:,fstaup1) = anompb_bt(:,:,fstaup1) &
         + dtbt*LAP_T(tmp(:,:),smooth_bih(:,:))
     endif

     ! accumulate for time average 
     if(barotropic_time_stepping_mom4p0) then 
         do j=jsd,jed
            do i=isd,ied   
               Ext_mode%anompb_bar(i,j,taup1) = Ext_mode%anompb_bar(i,j,taup1) &
                                                + anompb_bt(i,j,fstaup1)
               Ext_mode%udrho(i,j,1,taup1)    = Ext_mode%udrho(i,j,1,taup1)    &
                                                + udrho_bt(i,j,1,fstaup1)
               Ext_mode%udrho(i,j,2,taup1)    = Ext_mode%udrho(i,j,2,taup1)    &
                                                + udrho_bt(i,j,2,fstaup1)
            enddo
         enddo
     else 
         do j=jsd,jed
            do i=isd,ied   
               Ext_mode%anompb_bar(i,j,taup1) = Ext_mode%anompb_bar(i,j,taup1) &
                                                + anompb_bt(i,j,fstaup1)
               Ext_mode%udrho(i,j,1,taup1)    = Ext_mode%udrho(i,j,1,taup1)    &
                                                + float(nts-itime+1)*udrho_bt(i,j,1,fstaup1)
               Ext_mode%udrho(i,j,2,taup1)    = Ext_mode%udrho(i,j,2,taup1)    &
                                                + float(nts-itime+1)*udrho_bt(i,j,2,fstaup1)
            enddo
         enddo
     endif

     ! take a sample from the centre of the barotropic loop 
     if(itime==int(nts/2)) then 
         if (id_eta_t_bt > 0) then
             used = send_data (id_eta_t_bt, eta_t_bt(:,:,fstau), &
                    Time%model_time, rmask=Grd%tmask(:,:,1),     &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_anompb_bt > 0) then
             used = send_data (id_anompb_bt, anompb_bt(:,:,fstau), &
                    Time%model_time, rmask=Grd%tmask(:,:,1),       &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_udrho_bt_lap > 0) then
             used = send_data (id_udrho_bt_lap, friction_lap(:,:,1), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_udrho_bt_bih > 0) then
             used = send_data (id_udrho_bt_bih, friction_bih(:,:,1), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_vdrho_bt_lap > 0) then
             used = send_data (id_vdrho_bt_lap, friction_lap(:,:,2), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
         if (id_vdrho_bt_bih > 0) then
             used = send_data (id_vdrho_bt_bih, friction_bih(:,:,2), &
                    Time%model_time, rmask=Grd%umask(:,:,1),         &
                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
         endif
     endif


  enddo ! end of barotropic integration itime=1,nts 


  if(have_obc) call ocean_obc_update_boundary(Ext_mode%anompb_bar(:,:,taup1),'T')

end subroutine pred_corr_barotropic_press
! </SUBROUTINE> NAME="pred_corr_barotropic_press"


!#######################################################################
! <SUBROUTINE NAME="eta_smooth_diagnosed">
!
! <DESCRIPTION>
!  Smooth eta_t for case when running with pressure models and wish 
!  to have a diagnostic eta field smoothed. 
! </DESCRIPTION>
!
subroutine eta_smooth_diagnosed(Time, eta)

  type(ocean_time_type),      intent(in)    :: Time
  real, dimension(isd:,jsd:), intent(inout) :: eta

  integer :: i,j
  tmp(:,:) = 0.0

  ! tmp has dimensions (m/s) 
  if(smooth_eta_diag_laplacian) then 
    tmp(:,:) =  LAP_T(eta(:,:),smooth_lap_diag(:,:), update=.true.)
  elseif(smooth_eta_diag_biharmonic) then 
    tmp(:,:) = -LAP_T(eta(:,:),smooth_bih_diag(:,:), update=.true.)
    call mpp_update_domains (tmp, Dom%domain2d)
    tmp(:,:) =  LAP_T(tmp(:,:),smooth_bih_diag(:,:), update=.true.)
  endif 
  call mpp_update_domains (tmp, Dom%domain2d)

  do j=jsd,jed
     do i=isd,ied
        eta(i,j)             = eta(i,j) + tmp(i,j)*dtime
        eta_smooth_tend(i,j) = Grd%tmask(i,j,1)*dtime*tmp(i,j)
     enddo
  enddo  

  if (id_eta_smoother > 0) then 
     used = send_data (id_eta_smoother, tmp(:,:)*dtimer,& 
     Time%model_time, rmask=Grd%tmask(:,:,1),           &
     is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
  endif 

end subroutine eta_smooth_diagnosed
! </SUBROUTINE> NAME="eta_smooth_diagnosed"


!#######################################################################
! <SUBROUTINE NAME="ocean_eta_smooth">
!
! <DESCRIPTION>
!  Compute tendency for smoothing eta and tracer. 
!
!  Use either a laplacian or a biharmonic smoothing operator. 
!  Recommend against using the biharmonic, since it is not a 
!  positive definite operator and so can lead to extrema. 
! </DESCRIPTION>
!
subroutine ocean_eta_smooth(Time, Thickness, Ext_mode, T_prog)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)

  integer :: i, j, n, nprog, taum1
  real    :: eta_min

  nprog = size(T_prog(:))  

  ! initialise some fields for later use with OBC and return 
  if(.not. smooth_eta_t_laplacian .and. .not. smooth_eta_t_biharmonic) then
      if(have_obc) then
          Ext_mode%eta_smooth  = 0.0
          do n=1,nprog
             T_prog(n)%eta_smooth = 0.
          enddo
      endif
      return 
  endif

  taum1            = Time%taum1
  smooth_mask(:,:) = 0.0
  etastar(:,:)     = 0.0
  tmp(:,:)         = 0.0
  diffxy(:,:,:)    = 0.0  

  do n=1,nprog
     do j=jsd,jed
        do i=isd,ied
           T_prog(n)%eta_smooth(i,j) = 0.0
        enddo
     enddo
  enddo

  ! define etastar to be > 0 in order to have postive tracer diffusion 
  eta_min = mpp_global_min(Dom%domain2d,Ext_mode%eta_t(:,:,taum1))
  eta_min = abs(eta_min)
  do j=jsd,jed
     do i=isd,ied
        etastar(i,j) = Ext_mode%eta_t(i,j,taum1) &
                       + Grd%tmask(i,j,1)*(eta_min + eta_offset) 
     enddo
  enddo

  do j=jsd,jec
     do i=isd,iec
        diffxy(i,j,1) = etastar(i+1,j)-etastar(i,j)
        diffxy(i,j,2) = etastar(i,j+1)-etastar(i,j)
     enddo
  enddo

  do j=jsc,jec
     do i=isc,iec
        if(diffxy(i-1,j,1) /= 0.0 .or. diffxy(i,j,1) /= 0.0 .or. &
           diffxy(i,j-1,2) /= 0.0 .or. diffxy(i,j,2) /= 0.0 ) then 
           smooth_mask(i,j) = Grd%tmask(i,j,1)
        endif
     enddo
  enddo
  call mpp_update_domains (smooth_mask, Dom%domain2d)
   

  ! tmp has dimensions (m/s) 
  if(smooth_eta_t_laplacian) then 
    tmp(:,:) =  LAP_T(etastar(:,:),smooth_lap(:,:), update=.true.)
  else 
    tmp(:,:) = -LAP_T(etastar(:,:),smooth_bih(:,:), update=.true.)
    call mpp_update_domains (tmp, Dom%domain2d)
    tmp(:,:) =  LAP_T(tmp(:,:),smooth_bih(:,:), update=.true.)
  endif 

  ! need mass_source on full data domain in order to compute Adv_vel%wrho_bu
  call mpp_update_domains (tmp, Dom%domain2d)
  do j=jsd,jed
     do i=isd,ied
        tmp(i,j)                     = rho0*tmp(i,j)*smooth_mask(i,j)
        Ext_mode%eta_smooth(i,j)     = tmp(i,j)              
        Thickness%mass_source(i,j,1) = Thickness%mass_source(i,j,1) + tmp(i,j)
     enddo
  enddo
  
  if (id_eta_smoother > 0) used = send_data (id_eta_smoother, tmp(:,:)*rho0r,  & 
                                  Time%model_time, rmask=Grd%tmask(:,:,1),     &
                                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

  ! T_prog%eta_smooth has dimensions tracer concentration * (kg/m^3)*(m/s).
  ! note that tracer filter is zero when eta_t is zero, as we wish since in 
  ! this case there is no need to apply a filter.
  if(smooth_eta_t_laplacian) then 
      do n=1,nprog
         do j=jsd,jed
            do i=isd,ied
               tmp(i,j) = rho0*etastar(i,j)*T_prog(n)%field(i,j,1,taum1)
            enddo
         enddo
         tmp(:,:) = LAP_T(tmp(:,:),smooth_lap(:,:), update=.true.)
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%eta_smooth(i,j) = tmp(i,j)*smooth_mask(i,j)
            enddo
         enddo
      enddo
  else 
      do n=1,nprog
         do j=jsd,jed
            do i=isd,ied
               tmp(i,j) = rho0*etastar(i,j)*T_prog(n)%field(i,j,1,taum1)
            enddo
         enddo
         tmp(:,:) = -LAP_T(tmp(:,:),smooth_bih(:,:), update=.true.)
         call mpp_update_domains (tmp, Dom%domain2d)
         tmp(:,:) =  LAP_T(tmp(:,:),smooth_bih(:,:), update=.true.)
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%eta_smooth(i,j) = tmp(i,j)*smooth_mask(i,j)
            enddo
         enddo
      enddo
  endif

  if (id_smooth_mask > 0) used  = send_data (id_smooth_mask, smooth_mask(:,:), & 
                                  Time%model_time, rmask=Grd%tmask(:,:,1),     &
                                  is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

end subroutine ocean_eta_smooth
! </SUBROUTINE> NAME="ocean_eta_smooth"


!#######################################################################
! <SUBROUTINE NAME="ocean_pbot_smooth">
!
! <DESCRIPTION>
!  Compute tendency for diffusion of (pbot_t-pbot0) in both pbot_t 
!  and tracer. Need to consider tracer in order to maintain compability
!  between tracer and mass conservation equations.  
!
!  Use either a laplacian or a biharmonic smoother. 
!
!  Recommend against using the biharmonic, since it is 
!  NOT a positive definite operator and so can lead to extrema.
! </DESCRIPTION>
!
subroutine ocean_pbot_smooth(Time, Thickness, Ext_mode, T_prog)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(inout) :: Thickness
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)

  integer :: i, j, kb, n
  integer :: nprog, taum1
  real    :: press_min

  nprog = size(T_prog(:))  

  ! initialise some fields for later use with OBC and return 
  if(.not. smooth_pbot_t_laplacian .and. .not. smooth_pbot_t_biharmonic) then 
      if(have_obc) then
          Ext_mode%pbot_smooth = 0.0
          do n=1,nprog
             T_prog(n)%pbot_smooth = 0.0
          enddo
      endif
      return 
  endif

  taum1 = Time%taum1
  smooth_mask(:,:) = 0.0
  pbotstar(:,:)    = 0.0
  tmp(:,:)         = 0.0
  diffxy(:,:,:)    = 0.0

  do n=1,nprog
     do j=jsd,jed
        do i=isd,ied
           T_prog(n)%pbot_smooth(i,j) = 0.0
        enddo
     enddo
  enddo

  ! define pbotstar to be > 0 in order to have positive tracer diffusion 
  do j=jsd,jed
     do i=isd,ied
        tmp(i,j) = Grd%tmask(i,j,1)*(Ext_mode%pbot_t(i,j,taum1)-Thickness%pbot0(i,j))
     enddo
  enddo
  press_min = mpp_global_min(Dom%domain2d,tmp(:,:))
  press_min = abs(press_min)
  do j=jsd,jed
     do i=isd,ied
        pbotstar(i,j) = tmp(i,j) + Grd%tmask(i,j,1)*(press_min + pbot_offset) 
     enddo
  enddo

  do j=jsd,jec
     do i=isd,iec
        diffxy(i,j,1) = Grd%tmask(i+1,j,1)*Grd%tmask(i,j,1)*(pbotstar(i+1,j)-pbotstar(i,j))
        diffxy(i,j,2) = Grd%tmask(i,j+1,1)*Grd%tmask(i,j,1)*(pbotstar(i,j+1)-pbotstar(i,j))
     enddo
  enddo

  do j=jsc,jec
     do i=isc,iec
        if(diffxy(i-1,j,1) /= 0.0 .or. diffxy(i,j,1) /= 0.0 .or. &
           diffxy(i,j-1,2) /= 0.0 .or. diffxy(i,j,2) /= 0.0 ) then 
           smooth_mask(i,j) = Grd%tmask(i,j,1)
        endif
     enddo
  enddo
  call mpp_update_domains (smooth_mask, Dom%domain2d)

  if(smooth_pbot_t_laplacian) then 
    tmp(:,:) =  LAP_T(pbotstar(:,:),smooth_lap(:,:), update=.true.)
  else 
    tmp(:,:) = -LAP_T(tmp(:,:),smooth_bih(:,:), update=.true.)
    call mpp_update_domains (tmp, Dom%domain2d)
    tmp(:,:) =  LAP_T(tmp(:,:),smooth_bih(:,:), update=.true.)
  endif 

  ! need mass_source on full data domain in order to compute Adv_vel%wrho_bu.
  ! pbot_smooth has dimensions (kg/m^3)*(m/s)
  call mpp_update_domains (tmp, Dom%domain2d)

  do j=jsd,jed
     do i=isd,ied
        kb=Grd%kmt(i,j)
        if(kb>0) then 
            Ext_mode%pbot_smooth(i,j)     = tmp(i,j)*smooth_mask(i,j)*grav_r
            Thickness%mass_source(i,j,kb) = Thickness%mass_source(i,j,kb)  &
                                            + Ext_mode%pbot_smooth(i,j)
        endif 
     enddo
  enddo

  if (id_pbot_smoother > 0) then
   used  = send_data (id_pbot_smoother, grav*c2dbars*Ext_mode%pbot_smooth(:,:), & 
                      Time%model_time, rmask=Grd%tmask(:,:,1),                  &
                      is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
  endif 


  ! T_prog%pbot_smooth has dimensions tracer concentration * (kg/m^3)*(m/s) 
  if(smooth_pbot_t_laplacian) then 
      do n=1,nprog
         tmp(:,:)=0.0
         do j=jsd,jed
            do i=isd,ied
               kb=Grd%kmt(i,j)
               if(kb>0) then 
                   tmp(i,j) = Grd%tmask(i,j,1)*grav_r*pbotstar(i,j)*T_prog(n)%field(i,j,kb,taum1)
               endif
            enddo
         enddo
         tmp(:,:) = LAP_T(tmp(:,:),smooth_lap(:,:), update=.true.)
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%pbot_smooth(i,j) = tmp(i,j)*smooth_mask(i,j)
            enddo
         enddo
      enddo
  else 
      do n=1,nprog
         tmp(:,:)=0.0
         do j=jsd,jed
            do i=isd,ied
               kb=Grd%kmt(i,j)
               if(kb>0) then 
                   tmp(i,j) = Grd%tmask(i,j,1)*grav_r*pbotstar(i,j)*T_prog(n)%field(i,j,kb,taum1)
               endif
            enddo
         enddo
         tmp(:,:) = -LAP_T(tmp(:,:),smooth_bih(:,:), update=.true.)
         call mpp_update_domains (tmp, Dom%domain2d)
         tmp(:,:) =  LAP_T(tmp(:,:),smooth_bih(:,:), update=.true.)
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%pbot_smooth(i,j) = tmp(i,j)*smooth_mask(i,j)
            enddo
         enddo
      enddo
  endif

  if (id_smooth_mask > 0) then 
    used  = send_data (id_smooth_mask, smooth_mask(:,:), & 
    Time%model_time, rmask=Grd%tmask(:,:,1),             &
    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
  endif 

end subroutine ocean_pbot_smooth
! </SUBROUTINE> NAME="ocean_pbot_smooth"


!#######################################################################
! <SUBROUTINE NAME="barotropic_integrals">
!
! <DESCRIPTION>
! Compute area averaged fresh water and surface height and ocean mass.
! </DESCRIPTION>
!
subroutine barotropic_integrals (Time, Ext_mode, patm, pme, river)
  
  type(ocean_time_type),          intent(in) :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  real, dimension(isd:,jsd:),     intent(in) :: patm
  real, dimension(isd:,jsd:),     intent(in) :: pme
  real, dimension(isd:,jsd:),     intent(in) :: river 

  real    :: pme_total, river_total, etat_avg, mass_total
  integer :: i, j, tau
  type(time_type) :: next_time

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_barotropic_mod (barotropic_integrals): module must be initialized')
  endif 

  next_time = increment_time(Time%model_time, int(dtts), 0)  
  
  if (diag_step > 0 .and. mod(Time%itt, diag_step) == 0&
       .or. need_data(id_pme_total,next_time) .or. need_data(id_river_total,next_time)&
       .or. need_data(id_etat_avg,next_time).or. need_data(id_mass_total,next_time)) then  

      etat_avg    = 0.0 
      pme_total   = 0.0       
      river_total = 0.0 
      mass_total  = 0.0

      tau = Time%tau

      if(have_obc) then
        do j=jsc,jec
          do i=isc,iec
            etat_avg    = etat_avg    + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)*Grd%obc_tmask(i,j)
            pme_total   = pme_total   + Grd%dat(i,j)*Grd%tmask(i,j,1)*pme(i,j)*Grd%obc_tmask(i,j)
            river_total = river_total + Grd%dat(i,j)*Grd%tmask(i,j,1)*river(i,j)*Grd%obc_tmask(i,j)
            mass_total  = mass_total  + Grd%dat(i,j)*Grd%tmask(i,j,1) &
                                       *(Ext_mode%pbot_t(i,j,tau)-patm(i,j))*Grd%obc_tmask(i,j)
          enddo
        enddo
      else
        do j=jsc,jec
          do i=isc,iec
            etat_avg    = etat_avg    + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)
            pme_total   = pme_total   + Grd%dat(i,j)*Grd%tmask(i,j,1)*pme(i,j)
            river_total = river_total + Grd%dat(i,j)*Grd%tmask(i,j,1)*river(i,j)
            mass_total  = mass_total  + Grd%dat(i,j)*Grd%tmask(i,j,1) &
                                       *(Ext_mode%pbot_t(i,j,tau)-patm(i,j))
          enddo
        enddo
      endif
      call mpp_sum (etat_avg)
      call mpp_sum (pme_total)
      call mpp_sum (river_total)
      call mpp_sum (mass_total)

      etat_avg    = etat_avg/Grd%tcella(1)
      pme_total   = pme_total*dtuv
      river_total = river_total*dtuv

      if (id_etat_avg    > 0) used = send_data (id_etat_avg,    etat_avg,    Time%model_time)
      if (id_river_total > 0) used = send_data (id_river_total, river_total, Time%model_time)
      if (id_pme_total   > 0) used = send_data (id_pme_total,   pme_total,   Time%model_time)
      if (id_mass_total  > 0) used = send_data (id_mass_total,  mass_total,  Time%model_time)

  endif
  
end subroutine barotropic_integrals
! </SUBROUTINE> NAME="barotropic_integrals"



!#######################################################################
! <SUBROUTINE NAME="barotropic_energy">
!
! <DESCRIPTION>
!  Compute energetics of vertically integrated flow. 
! </DESCRIPTION>
!
subroutine barotropic_energy (Time, Ext_mode)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode

  real :: ke_bt, pe_bt, depth
  integer :: i, j, tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_barotropic_mod (barotropic_energy): module must be initialized')
  endif 

  tau = Time%tau

  ke_bt = 0.0 
  pe_bt = 0.0
  if(have_obc) then
    do j=jsc,jec
      do i=isc,iec
        depth = (Grd%hu(i,j) + Ext_mode%eta_u(i,j,tau))*Grd%obc_umask(i,j) + epsln
        ke_bt = ke_bt + Grd%dau(i,j)  &
                *((Ext_mode%udrho(i,j,1,tau))**2 + (Ext_mode%udrho(i,j,2,tau))**2)*Grd%obc_umask(i,j)/depth
        pe_bt = pe_bt + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)**2 * Grd%obc_tmask(i,j)
      enddo
    enddo
  else
    do j=jsc,jec
      do i=isc,iec
        depth = Grd%hu(i,j) + Ext_mode%eta_u(i,j,tau) + epsln
        ke_bt = ke_bt + Grd%dau(i,j)  &
                *((Ext_mode%udrho(i,j,1,tau))**2 + (Ext_mode%udrho(i,j,2,tau))**2)/depth
        pe_bt = pe_bt + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)**2
      enddo
    enddo
  endif

  call mpp_sum (ke_bt)
  call mpp_sum (pe_bt)

  ke_bt = ke_bt*0.5*rho0r
  pe_bt = pe_bt*grav*rho0

  if (id_ke_bt > 0) used = send_data (id_ke_bt, ke_bt/1.e15, Time%model_time)
  if (id_pe_bt > 0) used = send_data (id_pe_bt, pe_bt/1.e15, Time%model_time)

end subroutine barotropic_energy
! </SUBROUTINE> NAME="barotropic_energy"


!#######################################################################
! <SUBROUTINE NAME="read_barotropic">
!
! <DESCRIPTION>
!  Read in external mode fields from restart file.
! </DESCRIPTION>
subroutine read_barotropic(Time, Ext_mode)

  type(ocean_time_type),          intent(in)           :: Time
  type(ocean_external_mode_type), intent(inout)        :: Ext_mode

  character*128 file_name
  integer               :: tau, taum1, taup1
  integer, dimension(4) :: siz 

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  file_name = 'ocean_barotropic.res.nc'
  if(tendency==THREE_LEVEL) then 
      id_restart(1) = register_restart_field(Bar_restart, file_name, 'eta_t', Ext_mode%eta_t(:,:,tau),&
               Ext_mode%eta_t(:,:,taup1), domain=Dom%domain2d)
      id_restart(2) = register_restart_field(Bar_restart, file_name, 'anompb',Ext_mode%anompb(:,:,tau), &
               Ext_mode%anompb(:,:,taup1), domain=Dom%domain2d)       
      id_restart(3) = register_restart_field(Bar_restart, file_name, 'conv_rho_ud_t',Ext_mode%conv_rho_ud_t(:,:,tau), &
               Ext_mode%conv_rho_ud_t(:,:,taup1), domain=Dom%domain2d)  
      id_restart(4) = register_restart_field(Bar_restart, file_name, 'eta_t_bar',Ext_mode%eta_t_bar(:,:,tau), &
               Ext_mode%eta_t_bar(:,:,taup1), domain=Dom%domain2d)  
      id_restart(5) = register_restart_field(Bar_restart, file_name, 'anompb_bar',Ext_mode%anompb_bar(:,:,tau), &
               Ext_mode%anompb_bar(:,:,taup1), domain=Dom%domain2d)  
      id_restart(6) = register_restart_field(Bar_restart, file_name, 'eta_u',Ext_mode%eta_u(:,:,tau), &
               Ext_mode%eta_u(:,:,taup1), domain=Dom%domain2d)  
      id_restart(7) = register_restart_field(Bar_restart, file_name, 'pbot_u',Ext_mode%pbot_u(:,:,tau), &
               Ext_mode%pbot_u(:,:,taup1), domain=Dom%domain2d)  
      id_restart(8) = register_restart_field(Bar_restart, file_name, 'patm_t',Ext_mode%patm_t(:,:,tau), &
               Ext_mode%patm_t(:,:,taup1), domain=Dom%domain2d)  
      id_restart(9) = register_restart_field(Bar_restart, file_name, 'udrho',Ext_mode%udrho(:,:,1,tau), &
               Ext_mode%udrho(:,:,1,taup1), domain=Dom%domain2d) 
      id_restart(10)= register_restart_field(Bar_restart, file_name, 'vdrho',Ext_mode%udrho(:,:,2,tau), &
               Ext_mode%udrho(:,:,2,taup1), domain=Dom%domain2d) 
      id_restart(11)= register_restart_field(Bar_restart, file_name, 'eta_nonbouss',Ext_mode%eta_nonbouss(:,:,tau),&
               Ext_mode%eta_nonbouss(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(12)= register_restart_field(Bar_restart, file_name, 'eta_nonsteric',Ext_mode%eta_nonsteric(:,:,tau),&
               Ext_mode%eta_nonsteric(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(13)= register_restart_field(Bar_restart, file_name, 'eta_steric',Ext_mode%eta_steric(:,:,tau),&
               Ext_mode%eta_steric(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(14)= register_restart_field(Bar_restart, file_name, 'eta_dynamic',Ext_mode%eta_dynamic(:,:,tau),&
               Ext_mode%eta_dynamic(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(15)= register_restart_field(Bar_restart, file_name, 'eta_water',Ext_mode%eta_water(:,:,tau),&
               Ext_mode%eta_water(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(16)= register_restart_field(Bar_restart, file_name, 'eta_source',Ext_mode%eta_source(:,:,tau),&
               Ext_mode%eta_source(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(17)= register_restart_field(Bar_restart, file_name, 'eta_surf_temp',Ext_mode%eta_surf_temp(:,:,tau),&
               Ext_mode%eta_surf_temp(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(18)= register_restart_field(Bar_restart, file_name, 'eta_surf_salt',Ext_mode%eta_surf_salt(:,:,tau),&
               Ext_mode%eta_surf_salt(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(19)= register_restart_field(Bar_restart, file_name, 'eta_surf_water',Ext_mode%eta_surf_water(:,:,tau),&
               Ext_mode%eta_surf_water(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(20)= register_restart_field(Bar_restart, file_name, 'eta_bott_temp',Ext_mode%eta_bott_temp(:,:,tau),&
               Ext_mode%eta_bott_temp(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
  elseif(tendency==TWO_LEVEL) then    
      id_restart(1) = register_restart_field(Bar_restart, file_name, 'eta_t', &
               Ext_mode%eta_t(:,:,taup1), domain=Dom%domain2d)
      id_restart(2) = register_restart_field(Bar_restart, file_name, 'anompb', &
               Ext_mode%anompb(:,:,taup1), domain=Dom%domain2d)       
      id_restart(3) = register_restart_field(Bar_restart, file_name, 'conv_rho_ud_t', &
               Ext_mode%conv_rho_ud_t(:,:,taup1), domain=Dom%domain2d)  
      id_restart(4) = register_restart_field(Bar_restart, file_name, 'eta_t_bar', &
               Ext_mode%eta_t_bar(:,:,taup1), domain=Dom%domain2d)  
      id_restart(5) = register_restart_field(Bar_restart, file_name, 'anompb_bar', &
               Ext_mode%anompb_bar(:,:,taup1), domain=Dom%domain2d)  
      id_restart(6) = register_restart_field(Bar_restart, file_name, 'eta_u', &
               Ext_mode%eta_u(:,:,taup1), domain=Dom%domain2d)  
      id_restart(7) = register_restart_field(Bar_restart, file_name, 'pbot_u', &
               Ext_mode%pbot_u(:,:,taup1), domain=Dom%domain2d)  
      id_restart(8) = register_restart_field(Bar_restart, file_name, 'patm_t', &
               Ext_mode%patm_t(:,:,taup1), domain=Dom%domain2d)  
      id_restart(9) = register_restart_field(Bar_restart, file_name, 'udrho', &
               Ext_mode%udrho(:,:,1, taup1), domain=Dom%domain2d) 
      id_restart(10)= register_restart_field(Bar_restart, file_name, 'vdrho', &
               Ext_mode%udrho(:,:,2, taup1), domain=Dom%domain2d) 
      id_restart(11)= register_restart_field(Bar_restart, file_name, 'eta_nonbouss',              &
               Ext_mode%eta_nonbouss(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(12)= register_restart_field(Bar_restart, file_name, 'eta_nonsteric',              &
               Ext_mode%eta_nonsteric(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(13)= register_restart_field(Bar_restart, file_name, 'eta_steric',              &
               Ext_mode%eta_steric(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(14)= register_restart_field(Bar_restart, file_name, 'eta_dynamic',              &
               Ext_mode%eta_dynamic(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(15)= register_restart_field(Bar_restart, file_name, 'eta_water',               &
               Ext_mode%eta_water(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(16)= register_restart_field(Bar_restart, file_name, 'eta_source',              &
               Ext_mode%eta_source(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(17)= register_restart_field(Bar_restart, file_name, 'eta_surf_temp',              &
               Ext_mode%eta_surf_temp(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(18)= register_restart_field(Bar_restart, file_name, 'eta_surf_salt',              &
               Ext_mode%eta_surf_salt(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(19)= register_restart_field(Bar_restart, file_name, 'eta_surf_water',              &
               Ext_mode%eta_surf_water(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
      id_restart(20)= register_restart_field(Bar_restart, file_name, 'eta_bott_temp',              &
               Ext_mode%eta_bott_temp(:,:,taup1), domain=Dom%domain2d, mandatory=.false.) 
  endif
  id_restart(21)= register_restart_field(Bar_restart, file_name, 'forcing_u_bt',Ext_mode%forcing_bt(:,:,1), &
       domain=Dom%domain2d) 
  id_restart(22)= register_restart_field(Bar_restart, file_name, 'forcing_v_bt',Ext_mode%forcing_bt(:,:,2), &
       domain=Dom%domain2d) 


  file_name = 'INPUT/ocean_barotropic.res.nc'
  if(.NOT. file_exist(trim(file_name)) ) return

  call restore_state(Bar_restart)

  if(tendency==THREE_LEVEL) then 

      write (stdoutunit,'(a)') &
      '  Reading THREE_LEVEL restart from INPUT/ocean_barotropic.res.nc'
      write (stdoutunit,'(a)') &
      '  If not an initial condition, then expect two time records for each restart field.'

      call mpp_update_domains(Ext_mode%eta_t(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%anompb(:,:,:), Dom%domain2d)
      Ext_mode%pbot_t(:,:,tau)   = Ext_mode%anompb(:,:,tau)   + rho0grav*Grd%ht(:,:) 
      Ext_mode%pbot_t(:,:,taup1) = Ext_mode%anompb(:,:,taup1) + rho0grav*Grd%ht(:,:) 
      call mpp_update_domains(Ext_mode%conv_rho_ud_t(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_t_bar(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%anompb_bar(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_u(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%pbot_u(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%patm_t(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%udrho(:,:,1,tau), Ext_mode%udrho(:,:,2,tau), &
                              Dom%domain2d,gridtype=BGRID_NE)
      call mpp_update_domains(Ext_mode%udrho(:,:,1,taup1),Ext_mode%udrho(:,:,2,taup1), &
                              Dom%domain2d,gridtype=BGRID_NE)
      call mpp_update_domains(Ext_mode%eta_nonbouss(:,:,:),   Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_nonsteric(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_steric(:,:,:),     Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_dynamic(:,:,:),    Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_water(:,:,:),      Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_source(:,:,:),     Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_surf_temp(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_surf_salt(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_surf_water(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_bott_temp(:,:,:),  Dom%domain2d)

      call mpp_update_domains(Ext_mode%forcing_bt(:,:,1), Ext_mode%forcing_bt(:,:,2),&
                              Dom%domain2d,gridtype=BGRID_NE)

  elseif(tendency==TWO_LEVEL) then 

      write (stdoutunit,'(/a)') &
      '  Reading TWO_LEVEL restart from INPUT/ocean_barotropic.res.nc'
      write (stdoutunit,'(a)')  &
      '  Expecting only one time record for each restart field.'

      call field_size(file_name,'eta_t', siz)
      if (siz(4) > 1) then
        write(stdoutunit,'(/a)') &
        '==>ERROR: Attempt to read ocean_barotropic.res.nc from 3-level time scheme (2 time records)'
        write(stdoutunit,'(a)')  &
        ' when running mom4 with 2-level timestepping (only need 1 time record in restart).'
        write(stdoutunit,'(a)')  &
        ' Reduce restart file to only a single time record in order to avoid confusion.'
        call mpp_error(FATAL, &
        'Reading 3-time level ocean_barotropic.res.nc (w/ 2 time records) while using 2-level (needs only 1 record)')
      endif

      call mpp_update_domains(Ext_mode%eta_t(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%anompb(:,:,:), Dom%domain2d)
      Ext_mode%pbot_t(:,:,taup1) = Ext_mode%anompb(:,:,taup1) + rho0grav*Grd%ht(:,:) 
      call mpp_update_domains(Ext_mode%conv_rho_ud_t(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_t_bar(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%anompb_bar(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_u(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%pbot_u(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%patm_t(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%udrho(:,:,1,taup1), Ext_mode%udrho(:,:,2,taup1), &
                              Dom%domain2d, gridtype=BGRID_NE)
      call mpp_update_domains(Ext_mode%eta_nonbouss(:,:,:),   Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_nonsteric(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_steric(:,:,:),     Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_dynamic(:,:,:),    Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_water(:,:,:),      Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_source(:,:,:),     Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_surf_temp(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_surf_salt(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_surf_water(:,:,:), Dom%domain2d)
      call mpp_update_domains(Ext_mode%eta_bott_temp(:,:,:),  Dom%domain2d)
      call mpp_update_domains(Ext_mode%forcing_bt(:,:,1), Ext_mode%forcing_bt(:,:,2),&
                              Dom%domain2d,gridtype=BGRID_NE)

  endif

  if(have_obc) then
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,tau),'M','s')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,tau),'M','i')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,tau),'Z','i')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,tau),'Z','s')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,taup1),'M','s')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,taup1),'M','i')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,1,taup1),'Z','i')
     call ocean_obc_update_boundary(Ext_mode%udrho(:,:,2,taup1),'Z','s')
     call ocean_obc_update_boundary(Ext_mode%eta_t(:,:,:), 'T')    
!!!     call ocean_obc_update_boundary(Ext_mode%eta_u(:,:,:), 'C')    
     call ocean_obc_update_boundary(Ext_mode%eta_t_bar(:,:,:), 'T')    
  endif

end subroutine read_barotropic
! </SUBROUTINE> NAME="read_barotropic"


!#######################################################################
! <SUBROUTINE NAME="ocean_barotropic_restart">
! <DESCRIPTION>
!
!  Write out restart files registered through register_restart_file
!  Call to reset_field_pointer only needed for fields with a time index. 
!
! </DESCRIPTION>
subroutine ocean_barotropic_restart(Time, Ext_mode, time_stamp)
  type(ocean_time_type),      intent(in)           :: Time
  type(ocean_external_mode_type), intent(inout)    :: Ext_mode
  character(len=*),           intent(in), optional :: time_stamp
   integer :: tau, taup1

  tau    = Time%tau
  taup1  = Time%taup1

  if(tendency==THREE_LEVEL) then 
     call reset_field_pointer(Bar_restart, id_restart(1),  Ext_mode%eta_t(:,:,tau), Ext_mode%eta_t(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(2),  Ext_mode%anompb(:,:,tau), Ext_mode%anompb(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(3),  Ext_mode%conv_rho_ud_t(:,:,tau), Ext_mode%conv_rho_ud_t(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(4),  Ext_mode%eta_t_bar(:,:,tau), Ext_mode%eta_t_bar(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(5),  Ext_mode%anompb_bar(:,:,tau), Ext_mode%anompb_bar(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(6),  Ext_mode%eta_u(:,:,tau), Ext_mode%eta_u(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(7),  Ext_mode%pbot_u(:,:,tau), Ext_mode%pbot_u(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(8),  Ext_mode%patm_t(:,:,tau), Ext_mode%patm_t(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(9),  Ext_mode%udrho(:,:,1,tau), Ext_mode%udrho(:,:,1,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(10), Ext_mode%udrho(:,:,2,tau), Ext_mode%udrho(:,:,2,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(11), Ext_mode%eta_nonbouss(:,:,tau), Ext_mode%eta_nonbouss(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(12), Ext_mode%eta_nonsteric(:,:,tau), Ext_mode%eta_nonsteric(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(13), Ext_mode%eta_steric(:,:,tau), Ext_mode%eta_steric(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(14), Ext_mode%eta_dynamic(:,:,tau), Ext_mode%eta_dynamic(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(15), Ext_mode%eta_water(:,:,tau), Ext_mode%eta_water(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(16), Ext_mode%eta_source(:,:,tau), Ext_mode%eta_source(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(17), Ext_mode%eta_surf_temp(:,:,tau), Ext_mode%eta_surf_temp(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(18), Ext_mode%eta_surf_salt(:,:,tau), Ext_mode%eta_surf_salt(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(19), Ext_mode%eta_surf_water(:,:,tau), Ext_mode%eta_surf_water(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(20), Ext_mode%eta_bott_temp(:,:,tau), Ext_mode%eta_bott_temp(:,:,taup1) )
  elseif(tendency==TWO_LEVEL) then
     call reset_field_pointer(Bar_restart, id_restart(1),  Ext_mode%eta_t(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(2),  Ext_mode%anompb(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(3),  Ext_mode%conv_rho_ud_t(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(4),  Ext_mode%eta_t_bar(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(5),  Ext_mode%anompb_bar(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(6),  Ext_mode%eta_u(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(7),  Ext_mode%pbot_u(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(8),  Ext_mode%patm_t(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(9),  Ext_mode%udrho(:,:,1,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(10), Ext_mode%udrho(:,:,2,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(11), Ext_mode%eta_nonbouss(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(12), Ext_mode%eta_nonsteric(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(13), Ext_mode%eta_steric(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(14), Ext_mode%eta_dynamic(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(15), Ext_mode%eta_water(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(16), Ext_mode%eta_source(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(17), Ext_mode%eta_surf_temp(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(18), Ext_mode%eta_surf_salt(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(19), Ext_mode%eta_surf_water(:,:,taup1) )
     call reset_field_pointer(Bar_restart, id_restart(20), Ext_mode%eta_bott_temp(:,:,taup1) )
  end if

  call save_restart(Bar_restart, time_stamp)

end subroutine ocean_barotropic_restart
! </SUBROUTINE> NAME="ocean_barotropic_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_barotropic_end">
!
! <DESCRIPTION>
!  Write out external mode fields to restart file. 
! </DESCRIPTION>
subroutine ocean_barotropic_end(Time, Ext_mode)
  
  type(ocean_time_type),          intent(in)           :: Time
  type(ocean_external_mode_type), intent(inout)        :: Ext_mode

  character*128 file_name
  integer :: tau, taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_barotropic_mod (ocean_barotropic_end): module must be initialized')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') &
    '==>Warning from ocean_barotropic_mod (ocean_barotropic_end): NO restart written.'
    call mpp_error(WARNING, &
    '==>Warning from ocean_barotropic_mod (ocean_barotropic_end): NO restart written.')
    return
  endif 

  call ocean_barotropic_restart(Time, Ext_mode)

  tau   = Time%tau
  taup1 = Time%taup1

  if(tendency==THREE_LEVEL) then 
      write(stdoutunit,*) ' '
      write(stdoutunit,*) &
      ' From ocean_barotropic_mod: ending external mode chksums at tau'
      call write_timestamp(Time%model_time)
      call barotropic_chksum(Ext_mode, tau)
  endif
  write(stdoutunit,*) ' '
  write(stdoutunit,*) &
  ' From ocean_barotropic_mod: ending external mode chksums at taup1'
  call write_timestamp(Time%model_time)
  call barotropic_chksum(Ext_mode, taup1)

  module_is_initialized = .FALSE.

  nullify(Grd)
  nullify(Dom)

end subroutine ocean_barotropic_end
! </SUBROUTINE> NAME="ocean_barotropic_end"


!#######################################################################
! <SUBROUTINE NAME="maximum_convrhoud">
!
! <DESCRIPTION>
! Compute maximum convergence(rho_ud,rho_vd).
! </DESCRIPTION>
!
subroutine maximum_convrhoud(Time, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  real    :: convrhoudt0, convrhoudu0, convrhoudt, convrhoudu, fudge
  integer :: i, j
  integer :: iconvrhoudt, jconvrhoudt, iconvrhoudu, jconvrhoudu
  integer :: tau

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_barotropic_mod (maximum_convrhoud): module needs initialization ')
  endif 

  tau = Time%tau

  ! find max convergence(rho_ud,rho_vd) on T-cells

  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when convrhoud is independent of processor
  convrhoudt=0.0; iconvrhoudt=isc; jconvrhoudt=jsc

  do j=jsc,jec
    do i=isc,iec
      
      if (abs(Ext_mode%conv_rho_ud_t(i,j,tau)) > abs(convrhoudt) .and. Grd%tmask(i,j,1) == 1.0) then
        convrhoudt  = Ext_mode%conv_rho_ud_t(i,j,tau)
        iconvrhoudt = i
        jconvrhoudt = j
      endif
    enddo
  enddo

  ! find max convergence(ud,vd) on U-cells
  tmp = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%conv_rho_ud_t(:,:,tau))

  convrhoudu=0.0; iconvrhoudu=isc; jconvrhoudu=jsc
  do j=jsc,jec
    do i=isc,iec
      if (abs(tmp(i,j)) > abs(convrhoudu) .and. Grd%umask(i,j,1) == 1.0) then 
        convrhoudu  = tmp(i,j)
        iconvrhoudu = i
        jconvrhoudu = j
      endif
    enddo
  enddo
  
  write (stdoutunit,'(//60x,a/)') &
  ' Convergence of depth integrated horz velocity summary:'

  convrhoudt  = convrhoudt*fudge
  convrhoudu  = convrhoudu*fudge
  convrhoudt0 = convrhoudt
  convrhoudu0 = convrhoudu
  convrhoudt  = abs(convrhoudt)
  convrhoudu  = abs(convrhoudu)
  call mpp_max(convrhoudt)
  call mpp_max(convrhoudu)

  if (abs(convrhoudt0) == convrhoudt .and. abs(convrhoudt0) /= 0.0) then
    convrhoudt = convrhoudt0
    write (unit,9112) convrhoudt/fudge, iconvrhoudt, jconvrhoudt, &
           Grd%xt(iconvrhoudt,jconvrhoudt), Grd%yt(iconvrhoudt,jconvrhoudt)
  endif
  
  if (abs(convrhoudu0) == convrhoudu .and. abs(convrhoudu0) /= 0.0) then
    convrhoudu = convrhoudu0
    write (unit,9113) convrhoudu/fudge, iconvrhoudu, jconvrhoudu, &
           Grd%xu(iconvrhoudu,jconvrhoudu), Grd%yu(iconvrhoudu,jconvrhoudu)
  endif


9112  format(/' Maximum at T-cell convrhoud_t (',es10.3,' (kg/m^3)*m/s) at (i,j) = ','(',i4,',',i4,'),',&
          ' (lon,lat) = (',f7.2,',',f7.2,')')
9113  format(' Maximum at U-cell convrhoud_u (',es10.3,' (kg/m^3)*m/s) at (i,j) = ','(',i4,',',i4,'),',&
          ' (lon,lat) = (',f7.2,',',f7.2,')'/)


end subroutine maximum_convrhoud
! </SUBROUTINE> NAME="maximum_convrhoud"


!#######################################################################
! <SUBROUTINE NAME="barotropic_chksum">
!
! <DESCRIPTION>
!  Compute checksum for external mode fields.
! </DESCRIPTION>
!
subroutine barotropic_chksum(Ext_mode, index)

  type(ocean_external_mode_type), intent(in) :: Ext_mode
  integer, intent(in) :: index

  integer :: stdoutunit 
  stdoutunit=stdout() 

  write(stdoutunit,*) 'chksum for eta_t        = ',mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_u        = ',mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for deta_dt      = ',mpp_chksum(Ext_mode%deta_dt(isc:iec,jsc:jec))
  write(stdoutunit,*) 'chksum for eta_t_bar    = ',mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,index))

  write(stdoutunit,*) 'chksum for pbot_t       = ',mpp_chksum(Ext_mode%pbot_t(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for pbot_u       = ',mpp_chksum(Ext_mode%pbot_u(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for dpbot_dt     = ',mpp_chksum(Ext_mode%dpbot_dt(isc:iec,jsc:jec))
  write(stdoutunit,*) 'chksum for anompb       = ',mpp_chksum(Ext_mode%anompb(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for anompb_bar   = ',mpp_chksum(Ext_mode%anompb_bar(isc:iec,jsc:jec,index))

  write(stdoutunit,*) 'chksum for patm_t       = ',mpp_chksum(Ext_mode%patm_t(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for dpatm_dt     = ',mpp_chksum(Ext_mode%dpatm_dt(isc:iec,jsc:jec))

  write(stdoutunit,*) 'chksum for ps           = ',mpp_chksum(Ext_mode%ps(isc:iec,jsc:jec))
  write(stdoutunit,*) 'chksum for grad_ps_1    = ',mpp_chksum(Ext_mode%grad_ps(isc:iec,jsc:jec,1))
  write(stdoutunit,*) 'chksum for grad_ps_2    = ',mpp_chksum(Ext_mode%grad_ps(isc:iec,jsc:jec,2))
  write(stdoutunit,*) 'chksum for grad_anompb_1= ',mpp_chksum(Ext_mode%grad_anompb(isc:iec,jsc:jec,1))
  write(stdoutunit,*) 'chksum for grad_anompb_2= ',mpp_chksum(Ext_mode%grad_anompb(isc:iec,jsc:jec,2))

  write(stdoutunit,*) 'chksum for udrho        = ',mpp_chksum(Ext_mode%udrho(isc:iec,jsc:jec,1,index))
  write(stdoutunit,*) 'chksum for vdrho        = ',mpp_chksum(Ext_mode%udrho(isc:iec,jsc:jec,2,index))
  write(stdoutunit,*) 'chksum for conv_rho_ud_t= ',mpp_chksum(Ext_mode%conv_rho_ud_t(isc:iec,jsc:jec,index))

  write(stdoutunit,*) 'chksum for source       = ',mpp_chksum(Ext_mode%source(isc:iec,jsc:jec))
  write(stdoutunit,*) 'chksum for eta smoother = ',mpp_chksum(Ext_mode%eta_smooth(isc:iec,jsc:jec))
  write(stdoutunit,*) 'chksum for pbot smoother= ',mpp_chksum(Ext_mode%pbot_smooth(isc:iec,jsc:jec))

  write(stdoutunit,*) 'chksum for eta_nonbouss   = ',mpp_chksum(Ext_mode%eta_nonbouss(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_nonsteric  = ',mpp_chksum(Ext_mode%eta_nonsteric(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_steric     = ',mpp_chksum(Ext_mode%eta_steric(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_dynamic    = ',mpp_chksum(Ext_mode%eta_dynamic(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_water      = ',mpp_chksum(Ext_mode%eta_water(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_source     = ',mpp_chksum(Ext_mode%eta_source(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_surf_temp  = ',mpp_chksum(Ext_mode%eta_surf_temp(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_surf_salt  = ',mpp_chksum(Ext_mode%eta_surf_salt(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_surf_water = ',mpp_chksum(Ext_mode%eta_surf_water(isc:iec,jsc:jec,index))
  write(stdoutunit,*) 'chksum for eta_bott_temp  = ',mpp_chksum(Ext_mode%eta_bott_temp(isc:iec,jsc:jec,index))

  write(stdoutunit,*) 'chksum for forcing_u_bt = ',mpp_chksum(Ext_mode%forcing_bt(isc:iec,jsc:jec,1))
  write(stdoutunit,*) 'chksum for forcing_v_bt = ',mpp_chksum(Ext_mode%forcing_bt(isc:iec,jsc:jec,2))

  write(stdoutunit,*) ' '

end subroutine barotropic_chksum 
! </SUBROUTINE> NAME="barotropic_chksum"


!#######################################################################
! <SUBROUTINE NAME="psi_compute">
!
! <DESCRIPTION>
!  Compute quasi-barotropic streamfunctions for diagnostic purposes.
!  When no fresh water and steady state, these two streamfunctions 
!  will be equal, and they will be equal to the rigid lid barotropic
!  streamfunction in the Boussinesq case.  
!
!  Original algorithm: Stephen.Griffies@noaa.gov
!  Modifications for parallel efficiency: Giang.Nong@noaa.gov
!
!  13MAR2007: Change units to 10^9 kg/s, which is a "mass Sv"
!  This is the natural unit of transport for a mass-based 
!  vertical coordinate model.  
!
!  Updated Dec2009 to be compatible with tx_trans and ty_trans
!  calculation.  
!
! </DESCRIPTION>
!
subroutine psi_compute(Time, Adv_vel, psiu, psiv)

  type(ocean_time_type),      intent(in)    :: Time
  type(ocean_adv_vel_type),   intent(in)    :: Adv_vel
  real, dimension(isd:,jsd:), intent(inout) :: psiu
  real, dimension(isd:,jsd:), intent(inout) :: psiv

  real     :: psiu2(isd:ied), psiv2(jsd:jed) 
  integer  :: lenx,leny
  integer  :: i,j,k,ii
  integer  :: tau

  tau   = Time%tau
  psiu  = 0.0
  psiv  = 0.0
  psiu2 = 0.0
  psiv2 = 0.0 
  lenx  = iec-isc+1
  leny  = jec-jsc+1

  wrk1_v2d(:,:,:) = 0.0 
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           wrk1_v2d(i,j,1) = wrk1_v2d(i,j,1) + Adv_vel%uhrho_et(i,j,k)*Grd%dyte(i,j)
           wrk1_v2d(i,j,2) = wrk1_v2d(i,j,2) + Adv_vel%vhrho_nt(i,j,k)*Grd%dxtn(i,j)
        enddo
     enddo
  enddo
  
  do j=jsc,jec
     do i=isc,iec
        psiu(i,j) = -wrk1_v2d(i,j,1)*Grd%tmask(i,j,1)*transport_convert
        psiv(i,j) =  wrk1_v2d(i,j,2)*Grd%tmask(i,j,1)*transport_convert
     enddo
  enddo

  if(Dom%joff == 0 .and. jsc==1) then 
    psiu(:,jsc)=0.0  
  endif
  do i=isc,iec
     do j=jsc+1,jec
        psiu(i,j) = psiu(i,j) + psiu(i,j-1)
     enddo
  enddo
  psiu2(:) = psiu(:,jec)

  ! send psiu2 to other PEs in the same column with higher rank 
  if(myid_y<size_y) then
     do ii = myid_y+1,size_y
        call mpp_send(psiu2(isc:iec),lenx,pelist_y(ii))
     enddo
  endif
  call mpp_sync_self

  ! receive psiu2 from other PEs in the same column with lower rank
  if(myid_y>1) then
     do ii = 1,myid_y-1
        call mpp_recv(psiu2(isc:iec),lenx,pelist_y(ii))
         do i=isc,iec
            do j=jsc,jec
               psiu(i,j) = psiu(i,j) + psiu2(i)
            enddo
         enddo  
     enddo
  endif

  ! compute psiv
  if(Dom%ioff == 0 .and. isc==1 ) then 
      psiv(1,:) = psiu(1,:)
  endif  
  do j=jsc,jec
     do i=isc+1,iec
        psiv(i,j) = psiv(i,j) + psiv(i-1,j)
     enddo
  enddo
  psiv2(:) = psiv(iec,:)

  ! send psiv2 to other PEs in the same row with higher rank
  if(myid_x<size_x) then
     do ii = myid_x+1,size_x
        call mpp_send(psiv2(jsc:jec),leny,pelist_x(ii))
     enddo
  endif
  call mpp_sync_self

  ! receive psiv2 from other PEs in the same row with lower rank 
  if(myid_x>1) then
     do ii = 1,myid_x-1
        call mpp_recv(psiv2(jsc:jec),leny,pelist_x(ii))         
            do j=jsc,jec
               do i=isc,iec 
               psiv(i,j) = psiv(i,j) + psiv2(j)
            enddo
         enddo  
     enddo
  endif


end subroutine psi_compute
! </SUBROUTINE> NAME="psi_compute"


!#######################################################################
! <SUBROUTINE NAME="eta_terms_diagnose">
!
! <DESCRIPTION>
!  Diagnose various contributions to the sea level. 
!
!  WARNING: The steric diagnostics from this subroutine are confusing
!  when evaluated in a Boussinesq model.  The reason is that volume
!  conserving Boussinesq models have spurious mass sources, which 
!  corrupt the bottom pressure signal.  One needs to apply corrections
!  to make sense of the Boussinesq models for purposes of studying 
!  mass budgets, including the local contribution to steric effects. 
!
!
!  A/ contributions from dynamics, eustatic, and steric:
!
!  eta_nonbouss = (eta_dynamic + eta_water + eta_source) + eta_steric 
!                 = eta_nonsteric + eta_steric  
!
!  The time stepping diagnosed in this subroutine can lead to individual 
!  components to etanb that are quite huge.  The resulting sum, however,
!  should be quite close to eta_t.  In particular, locations of fresh 
!  water input can create eta_water that grows unboundedely, with 
!  a compenstating eta_dynamic that is negative, accounting for the 
!  dynamical adjustment accuring in the presence of water introduced to 
!  the ocean.  
! 
!  For PRESSURE_BASED vertical coordinates, eta_smooth_tend has
!  already been computed in subroutine eta_smooth_diagnosed.
!  We do not add this contribution to eta_nonbouss, since this 
!  piece is not part of the tendencies affecting bottom pressure.
!  It is only added for cosmetic reasons.  It is for this reason
!  that eta_smooth is NOT included in the restart file.  
!
!  For calculation of the steric contribution, a single time step
!  scheme is assumed, which is the recommended time stepping in 
!  MOM4p1.  
!   
!  For DEPTH_BASED models, the smoothing of eta is included in 
!  Ext_mode%source, so eta_smooth_tend is zero for depth-based models.
!
!  For PRESSURE_BASED vertical coordinates, eta_nonbouss as computed 
!  in this routine is very close to the prognostic eta_t.  Differences
!  arise from any possible smoothing applied to the diagnosed eta_t.  
!
!  B/ contributions from boundary fluxes, with residual due to 
!  SGS fluxes and nonlinear equation of state. 
!
! </DESCRIPTION>
!
subroutine eta_terms_diagnose(Time, Dens, Thickness, Theta, Salinity, Adv_vel, Ext_mode, pme, river)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_density_type),       intent(in)    :: Dens
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_prog_tracer_type),   intent(in)    :: Theta
  type(ocean_prog_tracer_type),   intent(in)    :: Salinity 
  type(ocean_adv_vel_type),       intent(in)    :: Adv_vel
  type(ocean_external_mode_type), intent(inout) :: Ext_mode 
  real, dimension(isd:,jsd:),     intent(in)    :: pme
  real, dimension(isd:,jsd:),     intent(in)    :: river

  integer :: i, j, k, kbot
  integer :: taum1, tau, taup1

  real    :: eta_global
  real    :: eta_steric_global 
  real    :: eta_nonsteric_global 
  real    :: eta_nonbouss_global 
  real    :: eta_dynamic_global 
  real    :: eta_water_global 
  real    :: eta_source_global 
  real    :: eta_surf_temp_global
  real    :: eta_surf_salt_global
  real    :: eta_surf_water_global
  real    :: eta_bott_temp_global
  real    :: rho_inv
  real    :: rhobarz_taup1, rhobarz_inv
  real    :: global_tmp
  real    :: volume_tmp
  real    :: factor 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  eta_dynamic_tend(:,:) = 0.0
  eta_water_tend(:,:)   = 0.0
  eta_steric_tend(:,:)  = 0.0
  eta_source_tend(:,:)  = 0.0
  rhodzt(:,:)           = 0.0
  rhodzt_inv(:,:)       = 0.0 
  rhodzt_taup1(:,:)     = 0.0


  ! compute boundary forcing from water, temp, and salt fluxes.  These fluxes, when
  ! globally integrated, contribute to time changes in global ocean volume.   
  ! wrk1_2d = surf water forcing in units of meter, with wrk1_2d > 0 increasing ocean volume 
  ! wrk2_2d = surf temp  forcing in units of meter, with wrk2_2d > 0 increasing ocean volume 
  !                 note: alpha = -rho_inv*Dens%drhodT > 0
  !                 note: stf > 0 means heat enters ocean 
  ! wrk3_2d = surf salt  forcing in units of meter, with wrk3_2d > 0 increasing ocean volume 
  !                 note: beta  =  rho_inv*Dens%drhodS > 0 (units 1/psu)
  ! wrk4_2d = bott temp  forcing in units of meter, with wrk4_2d > 0 increasing ocean volume 
  !                 note: alpha = -rho_inv*Dens%drhodT > 0
  !                 note: btf < 0 means heat enters ocean 
  wrk1_2d = 0.0
  wrk2_2d = 0.0
  wrk3_2d = 0.0
  wrk4_2d = 0.0
  do j=jsc,jec
     do i=isc,iec
        k            =  1
        rho_inv      =  Grd%tmask(i,j,k)/(epsln+Dens%rho(i,j,k,tau))
        wrk1_2d(i,j) =  dtime*(pme(i,j)+river(i,j))*rho_inv
        wrk2_2d(i,j) = -dtime*Dens%drhodT(i,j,k)*rho_inv*rho_inv*Theta%stf(i,j)
        wrk3_2d(i,j) = -dtime*Dens%drhodS(i,j,k)*rho_inv*rho_inv*Salinity%stf(i,j)

        kbot = Grd%kmt(i,j)
        if(kbot > 0) then 
           k            = kbot
           rho_inv      = Grd%tmask(i,j,k)/(epsln+Dens%rho(i,j,k,tau))
           wrk4_2d(i,j) = dtime*Dens%drhodT(i,j,kbot)*rho_inv*rho_inv*Theta%btf(i,j)
        endif 

     enddo
  enddo
  
  ! compute vertically integrated density
  if(vert_coordinate_class==PRESSURE_BASED) then 
      k=1
      do j=jsc,jec
         do i=isc,iec
            rhodzt(i,j)       = Grd%tmask(i,j,k)*grav_r*(Ext_mode%pbot_t(i,j,tau)  -Ext_mode%patm_t(i,j,tau))
            rhodzt_taup1(i,j) = Grd%tmask(i,j,k)*grav_r*(Ext_mode%pbot_t(i,j,taup1)-Ext_mode%patm_t(i,j,taup1))
         enddo
      enddo
  elseif(vert_coordinate_class==DEPTH_BASED) then
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               rhodzt(i,j)       = rhodzt(i,j)       + Grd%tmask(i,j,k)*rho0r*Thickness%rho_dzt(i,j,k,tau)   &
                                                       *Dens%rho(i,j,k,tau)
               rhodzt_taup1(i,j) = rhodzt_taup1(i,j) + Grd%tmask(i,j,k)*rho0r*Thickness%rho_dzt(i,j,k,taup1) &
                                                       *Dens%rho(i,j,k,taup1)
            enddo
         enddo
      enddo
  endif

  ! compute inverse 
  do j=jsc,jec
     do i=isc,iec
       rhodzt_inv(i,j) = Grd%tmask(i,j,1)/(epsln+rhodzt(i,j))
     enddo
  enddo


  wrk1_v2d(:,:,:) = 0.0

  ! diagnose tendency terms contributing to sea level changes 
  if(vert_coordinate_class==PRESSURE_BASED) then 
      do j=jsc,jec
         do i=isc,iec
            rhobarz_inv             =  rhodzt_inv(i,j)  *(Grd%ht(i,j) + Ext_mode%eta_t(i,j,tau))           
            eta_dynamic_tend(i,j)   =  Grd%tmask(i,j,1)*dtime*rhobarz_inv*Ext_mode%conv_rho_ud_t(i,j,tau)
            eta_water_tend(i,j)     =  Grd%tmask(i,j,1)*dtime*rhobarz_inv*(pme(i,j) + river(i,j))
            eta_source_tend(i,j)    =  Grd%tmask(i,j,1)*dtime*rhobarz_inv*Ext_mode%source(i,j)

            factor                  =  (1.0+Ext_mode%eta_t(i,j,tau)*Grd%htr(i,j)) &
                                      /(1.0+Ext_mode%eta_t(i,j,taup1)*Grd%htr(i,j))

            wrk1_v2d(i,j,1)         = factor 
            wrk1_v2d(i,j,2)         = rhodzt_taup1(i,j)*rhodzt_inv(i,j)

            eta_steric_tend(i,j)    =  Grd%tmask(i,j,1)*(Grd%ht(i,j) + Ext_mode%eta_t(i,j,tau)) &
                                      *(1.0-rhodzt_taup1(i,j)*rhodzt_inv(i,j)*factor)   

            eta_nonsteric_tend(i,j) =  eta_dynamic_tend(i,j) + eta_water_tend(i,j) &
                                      +eta_source_tend(i,j)
          enddo
      enddo
  elseif(vert_coordinate_class==DEPTH_BASED) then  
      do j=jsc,jec
         do i=isc,iec
            eta_dynamic_tend(i,j)   =  Grd%tmask(i,j,1)*dtime*rho0r*Ext_mode%conv_rho_ud_t(i,j,tau)
            eta_water_tend(i,j)     =  Grd%tmask(i,j,1)*dtime*rho0r*(pme(i,j) + river(i,j))
            eta_source_tend(i,j)    =  Grd%tmask(i,j,1)*dtime*rho0r*Ext_mode%source(i,j)

            factor                  =  (1.0+Ext_mode%eta_t(i,j,tau)*Grd%htr(i,j)) &
                                      /(1.0+Ext_mode%eta_t(i,j,taup1)*Grd%htr(i,j))
            eta_steric_tend(i,j)    =  Grd%tmask(i,j,1)*(Grd%ht(i,j) + Ext_mode%eta_t(i,j,tau)) &
                                      *(1.0-rhodzt_taup1(i,j)*rhodzt_inv(i,j)*factor)   

            eta_nonsteric_tend(i,j) =  eta_dynamic_tend(i,j) + eta_water_tend(i,j) &
                                      +eta_source_tend(i,j) 
         enddo
      enddo
  endif

  ! update eta_nonbouss and its contributions 
  do j=jsc,jec
     do i=isc,iec

        Ext_mode%eta_surf_water(i,j,taup1) = Ext_mode%eta_surf_water(i,j,taum1) &
            + Grd%dat(i,j)*area_total_r*wrk1_2d(i,j)

        Ext_mode%eta_surf_temp(i,j,taup1) = Ext_mode%eta_surf_temp(i,j,taum1) &
            + Grd%dat(i,j)*area_total_r*wrk2_2d(i,j)

        Ext_mode%eta_surf_salt(i,j,taup1) = Ext_mode%eta_surf_salt(i,j,taum1) &
            + Grd%dat(i,j)*area_total_r*wrk3_2d(i,j)

        Ext_mode%eta_bott_temp(i,j,taup1) = Ext_mode%eta_bott_temp(i,j,taum1) &
            + Grd%dat(i,j)*area_total_r*wrk4_2d(i,j)

        Ext_mode%eta_dynamic(i,j,taup1)  = Ext_mode%eta_dynamic(i,j,taum1) + eta_dynamic_tend(i,j)

        Ext_mode%eta_water(i,j,taup1)    = Ext_mode%eta_water(i,j,taum1)   + eta_water_tend(i,j)

        Ext_mode%eta_source(i,j,taup1)   = Ext_mode%eta_source(i,j,taum1)  + eta_source_tend(i,j)

        Ext_mode%eta_nonsteric(i,j,taup1)= Ext_mode%eta_nonsteric(i,j,taum1)  &
                                           + eta_dynamic_tend(i,j)            &
                                           + eta_water_tend(i,j)              &
                                           + eta_source_tend(i,j)            

        Ext_mode%eta_steric(i,j,taup1)   = Ext_mode%eta_steric(i,j,taum1) &
                                           + eta_steric_tend(i,j) 

        Ext_mode%eta_nonbouss(i,j,taup1) = Ext_mode%eta_nonbouss(i,j,taum1)              &
                                           + eta_dynamic_tend(i,j) + eta_water_tend(i,j) &
                                           + eta_source_tend(i,j)  + eta_steric_tend(i,j)
 
     enddo
  enddo
  call mpp_update_domains(Ext_mode%eta_nonbouss(:,:,taup1),    Dom%domain2d, complete=.false.)
  call mpp_update_domains(Ext_mode%eta_nonsteric(:,:,taup1),   Dom%domain2d, complete=.false.)
  call mpp_update_domains(Ext_mode%eta_steric(:,:,taup1),      Dom%domain2d, complete=.false.)
  call mpp_update_domains(Ext_mode%eta_dynamic(:,:,taup1),     Dom%domain2d, complete=.false.)
  call mpp_update_domains(Ext_mode%eta_water(:,:,taup1),       Dom%domain2d, complete=.false.)
  call mpp_update_domains(Ext_mode%eta_source(:,:,taup1),      Dom%domain2d, complete=.false. )
  call mpp_update_domains(Ext_mode%eta_surf_temp(:,:,taup1),   Dom%domain2d, complete=.false. )
  call mpp_update_domains(Ext_mode%eta_surf_salt(:,:,taup1),   Dom%domain2d, complete=.false. )
  call mpp_update_domains(Ext_mode%eta_surf_water(:,:,taup1),  Dom%domain2d, complete=.false. )
  call mpp_update_domains(Ext_mode%eta_bott_temp(:,:,taup1),   Dom%domain2d, complete=.true. )

  if(id_eta_nonbouss > 0)  used = send_data (id_eta_nonbouss, Ext_mode%eta_nonbouss(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                              &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_nonsteric > 0)  used = send_data (id_eta_nonsteric, Ext_mode%eta_nonsteric(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                                 &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_steric > 0)  used = send_data (id_eta_steric, Ext_mode%eta_steric(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                        &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_water > 0)  used = send_data (id_eta_water, Ext_mode%eta_water(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                     &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_dynamic > 0)  used = send_data (id_eta_dynamic, Ext_mode%eta_dynamic(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                           &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_source > 0)  used = send_data (id_eta_source, Ext_mode%eta_source(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                        &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if(id_eta_surf_temp > 0)  used = send_data (id_eta_surf_temp, Ext_mode%eta_surf_temp(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                                 &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_surf_salt > 0)  used = send_data (id_eta_surf_salt, Ext_mode%eta_surf_salt(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                                 &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec) 
  if(id_eta_surf_water > 0)  used = send_data (id_eta_surf_water, Ext_mode%eta_surf_water(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                                    &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_bott_temp > 0)  used = send_data (id_eta_bott_temp, Ext_mode%eta_bott_temp(:,:,tau),&
       Time%model_time, rmask=Grd%tmask(:,:,1),                                                 &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


  ! area averaged values eta and its pieces 
  if(id_eta_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = (Ext_mode%eta_t(i,j,tau)+grav_rho0r*Ext_mode%patm_t(i,j,tau))*Grd%dat(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_global, eta_global, Time%model_time)
  endif 
  if(id_eta_nonbouss_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_nonbouss(i,j,tau)*Grd%dat(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_nonbouss_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_nonbouss_global, eta_nonbouss_global, Time%model_time)
  endif 
  if(id_eta_nonsteric_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_nonsteric(i,j,tau)*Grd%dat(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_nonsteric_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_nonsteric_global, eta_nonsteric_global, Time%model_time)
  endif 
  if(id_eta_steric_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_steric(i,j,tau)*Grd%dat(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_steric_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_steric_global, eta_steric_global, Time%model_time)
  endif 
  if(id_eta_dynamic_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_dynamic(i,j,tau)*Grd%dat(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_dynamic_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_dynamic_global, eta_dynamic_global, Time%model_time)
  endif 
  if(id_eta_water_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_water(i,j,tau)*Grd%dat(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_water_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_water_global, eta_water_global, Time%model_time)
  endif 
  if(id_eta_source_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_source(i,j,tau)*Grd%dat(i,j)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_source_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_source_global, eta_source_global, Time%model_time)
  endif 


  ! area averaged values for boundary pieces of eta 
  if(id_eta_surf_water_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_surf_water(i,j,tau)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_surf_water_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)
      used = send_data (id_eta_surf_water_global, eta_surf_water_global, Time%model_time)
  endif 
  if(id_eta_surf_temp_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_surf_temp(i,j,tau)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_surf_temp_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)
      used = send_data (id_eta_surf_temp_global, eta_surf_temp_global, Time%model_time)
  endif 
  if(id_eta_surf_salt_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_surf_salt(i,j,tau)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_surf_salt_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)
      used = send_data (id_eta_surf_salt_global, eta_surf_salt_global, Time%model_time)
  endif 
  if(id_eta_bott_temp_global>0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Ext_mode%eta_bott_temp(i,j,tau)*Grd%tmask(i,j,1)
         enddo
      enddo
      eta_bott_temp_global = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)
      used = send_data (id_eta_bott_temp_global, eta_bott_temp_global, Time%model_time)
  endif 

  ! maps of tendencies 
 
  if(id_eta_tend > 0) used = send_data (id_eta_tend,                                         &
       Grd%tmask(:,:,1)*(Ext_mode%eta_nonbouss(:,:,taup1)-Ext_mode%eta_nonbouss(:,:,taum1)), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                              & 
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_dynamic_tend > 0) used = send_data (id_eta_dynamic_tend, eta_dynamic_tend(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                             &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_water_tend > 0)   used = send_data (id_eta_water_tend, eta_water_tend(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                         &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_nonsteric_tend > 0)  used = send_data (id_eta_nonsteric_tend, eta_nonsteric_tend(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                                    &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_source_tend > 0)  used = send_data (id_eta_source_tend, eta_source_tend(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                           &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_smooth_tend > 0)  used = send_data (id_eta_smooth_tend, eta_smooth_tend(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                           &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_steric_tend > 0)  used = send_data (id_eta_steric_tend, eta_steric_tend(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                           &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_steric_tend_A > 0)  used = send_data (id_eta_steric_tend_A, wrk1_v2d(:,:,1), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                          &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_steric_tend_B > 0)  used = send_data (id_eta_steric_tend_B, wrk1_v2d(:,:,2), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                          &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)


  ! maps of boundary forcing 
  if(id_eta_water_forcing > 0)  used = send_data (id_eta_water_forcing, wrk1_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                       &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_stemp_forcing > 0)  used = send_data (id_eta_stemp_forcing, wrk2_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                       &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_salt_forcing > 0)  used = send_data (id_eta_salt_forcing, wrk3_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                     &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if(id_eta_btemp_forcing > 0)  used = send_data (id_eta_btemp_forcing, wrk4_2d(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),                                       &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)



  ! area averaged tendencies 
  if(id_eta_tend_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j) &
                       *(Ext_mode%eta_nonbouss(i,j,taup1)-Ext_mode%eta_nonbouss(i,j,taum1))
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_tend_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_dynamic_tend_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*eta_dynamic_tend(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_dynamic_tend_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_water_tend_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*eta_water_tend(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_water_tend_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_steric_tend_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*eta_steric_tend(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_steric_tend_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_nonsteric_tend_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*eta_nonsteric_tend(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_nonsteric_tend_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_source_tend_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*eta_source_tend(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_source_tend_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_smooth_tend_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*eta_smooth_tend(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_smooth_tend_global, global_tmp, Time%model_time)
  endif 

  ! area averaged forcing 
  if(id_eta_water_forcing_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*wrk1_2d(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_water_forcing_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_stemp_forcing_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*wrk2_2d(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_stemp_forcing_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_salt_forcing_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*wrk3_2d(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_salt_forcing_global, global_tmp, Time%model_time)
  endif 
  if(id_eta_btemp_forcing_global > 0) then 
      do j=jsd,jed
         do i=isd,ied
            tmp(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*wrk4_2d(i,j)
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,tmp(:,:),NON_BITWISE_EXACT_SUM)*area_total_r
      used = send_data (id_eta_btemp_forcing_global, global_tmp, Time%model_time)
  endif 


  ! vertically averaged in-situ density 
  if(id_rhobarz > 0)  then 
      wrk1_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            if(Grd%tmask(i,j,1) > 0.0) then 
                wrk1_2d(i,j) = rhodzt(i,j)/(Grd%ht(i,j) + Ext_mode%eta_t(i,j,tau))
            endif
         enddo
      enddo
      used = send_data (id_rhobarz, wrk1_2d(:,:),  &
           Time%model_time, rmask=Grd%tmask(:,:,1),&
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! global average in-situ density 
  if(id_rhoavg > 0) then 
      wrk1_2d(:,:) = 0.0
      wrk2_2d(:,:) = 0.0
      do j=jsd,jed
         do i=isd,ied
            wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*rhodzt(i,j)
            wrk2_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*(Grd%ht(i,j) + Ext_mode%eta_t(i,j,tau))
         enddo
      enddo
      global_tmp = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:),NON_BITWISE_EXACT_SUM)
      volume_tmp = mpp_global_sum(Dom%domain2d,wrk2_2d(:,:),NON_BITWISE_EXACT_SUM)
      global_tmp = global_tmp/(epsln+volume_tmp)
      used = send_data (id_rhoavg, global_tmp, Time%model_time)
  endif 

  ! vertically integrated grav*w/sound_speed2.
  ! for vert_coor=pressure, this term gives a contribution to 
  ! the steric sea level evolution. For other coordinates,
  ! we are only approximating the contribution.  
  ! It should generally be very small. 
  if(id_eta_dpress_dt > 0) then 
      do j=jsc,jec
         do i=isc,iec
            wrk1_2d(i,j) = 0.0 
            do k=1,nk
               if(Grd%tmask(i,j,k) > 0.0) then 
                   wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%dzt(i,j,k) &
                    *(Adv_vel%wrho_bt(i,j,k)/Dens%rho(i,j,k,tau))*Dens%drhodP(i,j,k)
               endif
            enddo
            wrk1_2d(i,j) = wrk1_2d(i,j)*grav
         enddo
      enddo
      used = send_data (id_eta_dpress_dt, wrk1_2d(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


end subroutine eta_terms_diagnose
! </SUBROUTINE> NAME="eta_terms_diagnose"


!#######################################################################
! <SUBROUTINE NAME="eta_truncate">
!
! <DESCRIPTION>
!
!  Truncate eta_t to keep 
!
!  dzt(1) + eta_t >= frac_crit_cell_height*dzt(1)
!
!  and 
!
!  eta_t <= eta_max
!
! May be needed when run GEOPOTENTIAL models.
!
! </DESCRIPTION>
!
subroutine eta_truncate(Time, Ext_mode)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: i, j
  integer :: taup1
  real    :: cell_thickness

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_barotropic_mod (eta_truncate): module must be initialized')
  endif 

  taup1 = Time%taup1

  do j=jsc,jec
     do i=isc,iec
        cell_thickness = Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,taup1)+Grd%dzt(1)
        if (cell_thickness < Grd%dzt(1)*frac_crit_cell_height) then
            if (verbose_truncate) then 
                write(*,*) 'WARNING from ocean_barotropic_mod: eta_t(',i+Dom%ioff,',',j+Dom%joff,') = ', &
                Ext_mode%eta_t(i,j,taup1), '(m) has been truncated to', -Grd%dzt(1)*(1.0-frac_crit_cell_height)
            endif
            Ext_mode%eta_t(i,j,taup1) = -Grd%dzt(1)*(1.0-frac_crit_cell_height)
        elseif (Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,taup1) > eta_max) then
            if (verbose_truncate) then
                write(*,*) 'WARNING from ocean_barotropic_mod: eta_t(',i+Dom%ioff,',',j+Dom%joff,') = ', &
                Ext_mode%eta_t(i,j,taup1), '(m) has been truncated to', eta_max
            endif
            Ext_mode%eta_t(i,j,taup1) = eta_max
        endif
     enddo
  enddo


end subroutine eta_truncate
! </SUBROUTINE> NAME="eta_truncate"


!#######################################################################
! <SUBROUTINE NAME="eta_check">
!
! <DESCRIPTION>
! Perform diagnostic check on top cell thickness.  Useful when 
! when use GEOPOTENTIAL vertical coordinate.   
!
! </DESCRIPTION>
subroutine eta_check(Time, Ext_mode)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode

  integer :: i_thin_t_cell, j_thin_t_cell, i_thick_t_cell, j_thick_t_cell
  integer :: i_thin_u_cell, j_thin_u_cell, i_thick_u_cell, j_thick_u_cell
  integer :: tau, taup1
  integer :: i, j
  real :: thin_t_cell, thin_u_cell, thick_t_cell, thick_u_cell, cell_height, crit_cell_height
  real :: thin_t0, thin_u0, thick_t0, thick_u0, fudge

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_barotropic_mod (eta_check): module must be initialized')
  endif 

  tau   = Time%tau
  taup1 = Time%taup1

  thin_t_cell=Grd%dzt(1); thick_t_cell=Grd%dzt(1); thin_u_cell=Grd%dzt(1); thick_u_cell=Grd%dzt(1)
  
  i_thin_t_cell=isc;  j_thin_t_cell=jsc; i_thick_t_cell=isc; j_thick_t_cell=jsc
  i_thin_u_cell=isc;  j_thin_u_cell=jsc; i_thick_u_cell=isc; j_thick_u_cell=jsc
  crit_cell_height = frac_crit_cell_height*Grd%dzt(1)

  do j=jsc,jec
    do i=isc,iec
      cell_height = Ext_mode%eta_t(i,j,taup1)+Grd%dzt(1)    
      if (cell_height < thin_t_cell .and. Grd%tmask(i,j,1) == 1) then
        i_thin_t_cell = i
        j_thin_t_cell = j
        thin_t_cell   = cell_height
      endif 
    
      if (cell_height > thick_t_cell) then
        i_thick_t_cell = i
        j_thick_t_cell = j
        thick_t_cell   = cell_height
      endif 

      cell_height = Ext_mode%eta_u(i,j,taup1)+Grd%dzt(1)    
      if (cell_height < thin_u_cell .and. Grd%umask(i,j,1) == 1) then
        i_thin_u_cell = i
        j_thin_u_cell = j
        thin_u_cell   = cell_height
      endif 
    
      if (cell_height > thick_u_cell) then
        i_thick_u_cell = i
        j_thick_u_cell = j
        thick_u_cell   = cell_height
      endif         
    enddo
  enddo

  write (stdoutunit,'(//60x,a/)') ' Surface cell thickness summary:  '
  fudge = 1 + 1.e-12*mpp_pe()  ! to distinguish processors when tracer is independent of processor
  thin_t_cell  = thin_t_cell*fudge
  thin_u_cell  = thin_u_cell*fudge
  thick_t_cell = thick_t_cell*fudge
  thick_u_cell = thick_u_cell*fudge
  
  thin_t0  = thin_t_cell
  thin_u0  = thin_u_cell
  thick_t0 = thick_t_cell
  thick_u0 = thick_u_cell
  call mpp_min(thin_t_cell)
  call mpp_min(thin_u_cell)
  call mpp_max(thick_t_cell)
  call mpp_max(thick_u_cell)

  if (thin_t0 == thin_t_cell) then
    write (unit,7000) ' Minimum thickness surface T cell = ',thin_t_cell/fudge, &
                        i_thin_t_cell+Dom%ioff, j_thin_t_cell+Dom%joff, &
                        Grd%xt(i_thin_t_cell,j_thin_t_cell), Grd%yt(i_thin_t_cell,j_thin_t_cell)
  endif
  
  if (thick_t0 == thick_t_cell) then
    write (unit,7000) ' Maximum thickness surface T cell = ',thick_t_cell/fudge, &
                        i_thick_t_cell+Dom%ioff, j_thick_t_cell+Dom%joff, &
                        Grd%xt(i_thick_t_cell,j_thick_t_cell), Grd%yt(i_thick_t_cell,j_thick_t_cell)
  endif

  if (thin_u0 == thin_u_cell) then
    write (unit,7000) ' Minimum thickness surface U cell = ',thin_u_cell/fudge, &
                        i_thin_u_cell+Dom%ioff, j_thin_u_cell+Dom%joff, &
                        Grd%xu(i_thin_u_cell,j_thin_u_cell), Grd%yu(i_thin_u_cell,j_thin_u_cell)
  endif

  if (thick_u0 == thick_u_cell) then
    write (unit,7000) ' Maximum thickness surface U cell = ',thick_u_cell/fudge, &
                        i_thick_u_cell+Dom%ioff, j_thick_u_cell+Dom%joff, &
                        Grd%xu(i_thick_u_cell,j_thick_u_cell), Grd%yu(i_thick_u_cell,j_thick_u_cell)
  endif

  if (thin_t_cell < crit_cell_height) then
    write (stdoutunit,*) '==>Error: minimum surface cell thickness is less than ',crit_cell_height,' meters.'
    call mpp_error(FATAL, '==>Error: minimum surface cell thickness is below critical value.')
  endif 

  7000 format (1x,a,e14.7,' m at (i,j) = (',i4,',',i4,'), (lon,lat) = (',f7.2,',',f7.2,')')

end subroutine eta_check
! </SUBROUTINE> NAME="eta_check"



!#######################################################################
! <SUBROUTINE NAME="tidal_forcing_init">
!
! <DESCRIPTION>
! Initialize fields needed for lunar and solar tidal forcing.  
! </DESCRIPTION>
!
subroutine tidal_forcing_init(Time)

  type(ocean_time_type), intent(in) :: Time

  real :: xcenter, ycenter
  real :: xwidth, ywidth
  integer :: i,j

     allocate (eta_eq_tidal(isd:ied,jsd:jed))
     allocate (cos2lon(isd:ied,jsd:jed))
     allocate (coslon(isd:ied,jsd:jed))
     allocate (coslat2(isd:ied,jsd:jed))
     allocate (sin2lon(isd:ied,jsd:jed))
     allocate (sinlon(isd:ied,jsd:jed))
     allocate (sin2lat(isd:ied,jsd:jed))
     allocate (ideal_amp(isd:ied,jsd:jed))

     eta_eq_tidal = 0.0
     cos2lon      = 0.0 
     coslon       = 0.0 
     coslat2      = 0.0 
     sin2lon      = 0.0 
     sinlon       = 0.0 
     sin2lat      = 0.0
     ideal_amp    = 0.0

     if(tidal_forcing_m2) then
       call mpp_error(NOTE, &
        '==>Note from tidal_forcing_init: adding M2 lunar tidal forcing to ext-mode')
     endif 
     if(tidal_forcing_8) then
       call mpp_error(NOTE, &
       '==>Note from tidal_forcing_init: adding tidal forcing to ext-mode from 8-lunar/solar constituents')
     endif 
     if(tidal_forcing_m2 .and. tidal_forcing_8) then 
       call mpp_error(NOTE, &
       '==>Note from tidal_forcing_init: both tidal_forcing_m2 and tidal_forcing_8 = .true. Will use tidal_forcing_8.')
     endif 

     if(tidal_forcing_m2 .or. tidal_forcing_8 .or. tidal_forcing_ideal) then 
       tidal_forcing=.true.
     else 
       tidal_forcing=.false.
     endif 

     if(tidal_forcing) then 
       alphat = 0.948
     else 
       alphat = 1.0
     endif 

     tidal_omega_K1 = 0.72921e-4*3600.*24.     ! K1-tide
     tidal_omega_O1 = 0.67598e-4*3600.*24.     ! O1-tide
     tidal_omega_P1 = 0.72523e-4*3600.*24.     ! P1-tide
     tidal_omega_Q1 = 0.64959e-4*3600.*24.     ! Q1-tide

     tidal_omega_M2 = 1.40519e-4*3600.*24.     ! M2-tide
     tidal_omega_S2 = 1.45444e-4*3600.*24.     ! S2-tide
     tidal_omega_N2 = 1.37880e-4*3600.*24.     ! N2-tide
     tidal_omega_K2 = 1.45842e-4*3600.*24.     ! K2-tide

     Love_M2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_S2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_N2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_K2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_K1 =   1.0 + 0.256 - 0.520    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_O1 =   1.0 + 0.298 - 0.603    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_P1 =   1.0 + 0.287 - 0.581    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_Q1 =   1.0 + 0.298 - 0.603    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 

     amp_M2 =   0.242334 ! amplitude of M2 equilibrium tide (m)
     amp_S2 =   0.112743 ! amplitude of S2 equilibrium tide (m)
     amp_N2 =   0.046397 ! amplitude of N2 equilibrium tide (m)
     amp_K2 =   0.030684 ! amplitude of K2 equilibrium tide (m)
     amp_K1 =   0.141565 ! amplitude of K1 equilibrium tide (m)
     amp_O1 =   0.100661 ! amplitude of O1 equilibrium tide (m)
     amp_P1 =   0.046848 ! amplitude of P1 equilibrium tide (m)
     amp_Q1 =   0.019273 ! amplitude of Q1 equilibrium tide (m)

     rradian = 1./radian
     cos2lon(:,:) = cos(2.0*Grd%xt(:,:)*rradian)
     coslon (:,:) = cos(    Grd%xt(:,:)*rradian)
     coslat2(:,:) = cos(    Grd%yt(:,:)*rradian)**2
     sin2lon(:,:) = sin(2.0*Grd%xt(:,:)*rradian)
     sinlon (:,:) = sin(    Grd%xt(:,:)*rradian)
     sin2lat(:,:) = sin(2.0*Grd%yt(:,:)*rradian)


     ! Gaussian bump profile, set to zero near global boundaries  
     if(tidal_forcing_ideal) then 
         xcenter = Grd%grid_x_t(int(Grd%ni/2)) 
         xwidth  = Grd%grid_x_t(int(Grd%ni/4))
         ycenter = Grd%grid_y_t(int(Grd%nj/2)) 
         ywidth  = Grd%grid_y_t(int(Grd%nj/4))
         ideal_amp(:,:) = 0.0
         do j=jsd,jed
            do i=isd,ied
               if(i+Dom%ioff < Grd%ni-2 .and. i+Dom%ioff > 2 .and. &
                  j+Dom%joff < Grd%nj-2 .and. j+Dom%joff > 2) then 
               ideal_amp(i,j) = Love_M2*amp_M2 &
                *exp( -((Grd%xt(i,j)-xcenter)/xwidth)**2 - ((Grd%yt(i,j)-ycenter)/ywidth)**2 )
               endif 
            enddo
         enddo
     endif

     id_eta_eq_tidal  = register_diag_field ('ocean_model', 'eta_eq_tidal', Grd%tracer_axes(1:2), &
                        Time%model_time, 'equilibrium tidal potential', 'meter',                  &
                        missing_value=missing_value, range=(/-1e3,1e3 /)) 
 
end subroutine tidal_forcing_init
! </SUBROUTINE> NAME="tidal_forcing_init"


!#######################################################################
! <SUBROUTINE NAME="geoid_forcing_init">
!
! <DESCRIPTION>
! Initialize fields needed for modifying the geoid, relative to the 
! standard geoid. 
! </DESCRIPTION>
!
subroutine geoid_forcing_init(Time)

  type(ocean_time_type), intent(in) :: Time

  allocate (eta_geoid(isd:ied,jsd:jed))
  eta_geoid(:,:) = 0.0 
  
  if(geoid_forcing) then
      call mpp_error(NOTE, &
           '==>Note from geoid_forcing_init: modifying the geoid through a time independent tide-like forcing')

      if( file_exist('INPUT/eta_geoid.nc') ) then
          call read_data('INPUT/eta_geoid.nc', 'eta_geoid',eta_geoid(:,:), Dom%domain2d,timelevel=1)
          eta_geoid(:,:) = Grd%tmask(:,:,1)*eta_geoid(:,:)
          call mpp_update_domains (eta_geoid(:,:), Dom%domain2d)
      else 
          call mpp_error(FATAL,&
               '==>Error in ocean_barotropic_init: could not find INPUT/eta_geoid.nc, yet geoid_forcing=.true.')
      endif

  endif

  id_eta_geoid = register_static_field ('ocean_model', 'eta_geoid', Grd%tracer_axes(1:2), &
                   'perturbation of the geoid', 'metre',                                  &
                   missing_value=missing_value, range=(/-1e3,1e3/))
  if (id_eta_geoid > 0) then
     used = send_data(id_eta_geoid, eta_geoid(isc:iec,jsc:jec), Time%model_time)
  endif 
 

end subroutine geoid_forcing_init
! </SUBROUTINE> NAME="geoid_forcing_init"


!#######################################################################
! <SUBROUTINE NAME="get_tidal_forcing">
!
! <DESCRIPTION>
! Compute equilibrium tidal forcing.
! </DESCRIPTION>
!
subroutine get_tidal_forcing(Time, dayr)

  type(ocean_time_type), intent(in) :: Time 
  real,                  intent(in) :: dayr

  real    :: cosomegat_M2,sinomegat_M2
  real    :: cosomegat_S2,sinomegat_S2
  real    :: cosomegat_N2,sinomegat_N2
  real    :: cosomegat_K2,sinomegat_K2
  real    :: cosomegat_K1,sinomegat_K1
  real    :: cosomegat_O1,sinomegat_O1
  real    :: cosomegat_P1,sinomegat_P1
  real    :: cosomegat_Q1,sinomegat_Q1

  sinomegat_K1=sin(tidal_omega_K1*dayr)
  sinomegat_O1=sin(tidal_omega_O1*dayr)
  sinomegat_P1=sin(tidal_omega_P1*dayr)
  sinomegat_Q1=sin(tidal_omega_Q1*dayr)

  cosomegat_K1=cos(tidal_omega_K1*dayr)
  cosomegat_O1=cos(tidal_omega_O1*dayr)
  cosomegat_P1=cos(tidal_omega_P1*dayr)
  cosomegat_Q1=cos(tidal_omega_Q1*dayr)

  sinomegat_M2=sin(tidal_omega_M2*dayr)
  sinomegat_s2=sin(tidal_omega_S2*dayr)
  sinomegat_n2=sin(tidal_omega_N2*dayr)
  sinomegat_K2=sin(tidal_omega_K2*dayr)

  cosomegat_M2=cos(tidal_omega_M2*dayr)
  cosomegat_s2=cos(tidal_omega_S2*dayr)
  cosomegat_n2=cos(tidal_omega_N2*dayr)
  cosomegat_K2=cos(tidal_omega_K2*dayr)

  eta_eq_tidal(:,:) = 0.0

  ! M2 constituent only
  if(tidal_forcing_m2) then 
    eta_eq_tidal(:,:) = Love_M2*amp_M2*(coslat2(:,:))*(cosomegat_M2*cos2lon(:,:)-sinomegat_M2*sin2lon(:,:))
  endif 

  ! 8 principle constituents
  if(tidal_forcing_8) then 
     eta_eq_tidal(:,:) =  Love_M2*amp_M2*(coslat2(:,:))*(cosomegat_M2*cos2lon(:,:) - sinomegat_M2*sin2lon(:,:))+&
                          Love_S2*amp_S2*(coslat2(:,:))*(cosomegat_S2*cos2lon(:,:) - sinomegat_S2*sin2lon(:,:))+&
                          Love_N2*amp_N2*(coslat2(:,:))*(cosomegat_N2*cos2lon(:,:) - sinomegat_N2*sin2lon(:,:))+&
                          Love_K2*amp_K2*(coslat2(:,:))*(cosomegat_K2*cos2lon(:,:) - sinomegat_K2*sin2lon(:,:))+&
                          Love_K1*amp_K1*(sin2lat(:,:))*(cosomegat_K1*coslon (:,:) - sinomegat_K1*sinlon (:,:))+&
                          Love_O1*amp_O1*(sin2lat(:,:))*(cosomegat_O1*coslon (:,:) - sinomegat_O1*sinlon (:,:))+&
                          Love_P1*amp_P1*(sin2lat(:,:))*(cosomegat_P1*coslon (:,:) - sinomegat_P1*sinlon (:,:))+&
                          Love_Q1*amp_Q1*(sin2lat(:,:))*(cosomegat_Q1*coslon (:,:) - sinomegat_Q1*sinlon (:,:))
  endif 

  ! ideal tidal forcing
  if(tidal_forcing_ideal) then 
     eta_eq_tidal(:,:)  =  ideal_amp(:,:)*cosomegat_M2
  endif 


  if (id_eta_eq_tidal > 0) then 
    used = send_data (id_eta_eq_tidal, eta_eq_tidal(:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,1),       &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec) 
  endif 


end subroutine get_tidal_forcing
! </SUBROUTINE> NAME="get_tidal_forcing"


!#######################################################################
! <SUBROUTINE NAME="ideal_initialize_eta">
!
! <DESCRIPTION>
! Idealized initial condition for eta. 
! </DESCRIPTION>
!
subroutine ideal_initialize_eta

  integer :: i,j

  do j=jsd,jed
     do i=isd,ied
        eta_t_init(i,j) =  Grd%tmask(i,j,1)*ideal_initial_eta_amplitude     &
                          *sin(2.0*pi*Grd%xt(i,j)/ideal_initial_eta_xwidth) &
                          *sin(2.0*pi*Grd%yt(i,j)/ideal_initial_eta_ywidth) 
     enddo
  enddo



end subroutine ideal_initialize_eta
! </SUBROUTINE> NAME="ideal_initialize_eta"


!#######################################################################
! <FUNCTION NAME="REMAP_BT_TO_BU_LOCAL">
!
! <DESCRIPTION>
! Local version of the operator, needed here for initialization
! when read in eta information from an initialization file.
! Since barotropic is initialized prior to operators, we need
! to have this operator here locally.  
! </DESCRIPTION>
!
function REMAP_BT_TO_BU_LOCAL(a) 

    real, dimension(isd:,jsd:), intent(in) :: a
    real, dimension(isd:ied,jsd:jed)       :: REMAP_BT_TO_BU_LOCAL
    integer :: i, j

    do j=jsc-halo,jec+halo-1
       do i=isc-halo,iec+halo-1
          REMAP_BT_TO_BU_LOCAL(i,j) = (a(i,j)*Grd%dte(i,j)*Grd%dus(i,j) + a(i+1,j)*Grd%dtw(i+1,j)*Grd%dus(i,j)&
               + a(i,j+1)*Grd%dte(i,j+1)*Grd%dun(i,j) + a(i+1,j+1)*Grd%dtw(i+1,j+1)*Grd%dun(i,j))*Grd%daur(i,j)
       enddo
    enddo
    REMAP_BT_TO_BU_LOCAL(iec+halo,:) = 0.0
    REMAP_BT_TO_BU_LOCAL(:,jec+halo) = 0.0

end function REMAP_BT_TO_BU_LOCAL
! </FUNCTION> NAME="REMAP_BT_TO_BU_LOCAL"



end module ocean_barotropic_mod
