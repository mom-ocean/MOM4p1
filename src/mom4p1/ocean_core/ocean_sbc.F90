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
module ocean_sbc_mod
!  
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison
!</CONTACT>
!
! <REVIEWER EMAIL="Tony.Rosati@noaa.gov"> A. Rosati 
! </REVIEWER>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
! </REVIEWER>
!
! <REVIEWER EMAIL="V.Balaji@noaa.gov">
! V. Balaji
! </REVIEWER>
!
!<OVERVIEW>
! Set up the surface boundary conditions for mom4. 
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module sets up the surface boundary conditions for the model. 
! Also fill Ocean_sfc derived-type used to pass information to other 
! component models.
!
! The surface temperature should be the surface insitu temperature,
! which is the same as the surface potential temperature.  When the 
! model prognostic temperature variable is conservative temperature, 
! then the surface potential temperature is carried in T_diag(index_diag_temp).
! The resulting heat flux is potential enthalpy, which is the correct 
! field to be forcing the T_prog(index_temp) field when the prognostic
! temperature field is the conservative temperature.   
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! "Potential enthalpy: A conservative oceanic variable for evaluating
!  heat content and heat fluxes"
!  Trevor J McDougall, Journal of Physical Oceanography, 
!  vol 33, pages 945--963.
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_sbc_nml">
!
!  <DATA NAME="use_waterflux" TYPE="logical">
!  Set to true when wish to use real fresh water flux as opposed to virtual 
!  salt fluxes.   
!  </DATA> 
!  <DATA NAME="waterflux_tavg" TYPE="logical">
!  Set to true when aiming to suppress the leap-frog computational mode
!  by setting pme and river equal to a time averaged value over the 
!  present and previous time step.  This method requires an extra
!  field in the restart file.  This method is not needed when using
!  the TWO_LEVEL time tendency.  It remains for those who wish to 
!  use the leap-frog THREE_LEVEL time stepping scheme.  
!  Note that it does not lead to simple checks of conservation across
!  model components, since there is a time averaging performed for 
!  the water flux added to the ocean model.  It is generally NOT 
!  recommended.  Default waterflux_tavg=.false. 
!  </DATA> 
!
!  <DATA NAME="temp_restore_tscale" UNITS="day" TYPE="real">
!  Time scale in days for restoring temperature within the top model 
!  grid cell. 
!  </DATA> 
!  <DATA NAME="salt_restore_tscale" UNITS="day" TYPE="real">
!  Time scale in days for restoring salinity within the top model 
!  grid cell. 
!  </DATA> 
!  <DATA NAME="salt_restore_as_salt_flux" TYPE="logical">
!  When running a use_waterflux=.true. model, we may choose to add the 
!  salinity from a restoring condition as a salt flux or convert to 
!  a fresh water flux. The addition of salt does not alter the sea 
!  level nor does it alter the concentration of other tracers, whereas
!  converting to an implied water flux will alter sea level and other
!  concentrations.  So we generally recommend the default   
!  salt_restore_as_salt_flux=.true. 
!  </DATA> 
!  <DATA NAME="max_delta_salinity_restore" UNITS="ppt" TYPE="real">
!  When computing the restoring flux for salinity, we can define
!  a maximum absolute value for the difference between salinity(k=1)
!  and the restoring salinity from a dataset.  This approach is useful
!  especially in NAtl western boundary, where poor Gulf Stream separation
!  can lead to large salinity biases.  If restore too much the salinity
!  field, we can spuriously transport large amounts of fresh water to the 
!  subpoloar gyre, thus impacting the overturning circulation too much.  
!  If max_delta_salinity_restore < 0.0, then will NOT provide a max to the 
!  delta salinity; will instead compute an unbounded restoring flux.  
!  Default max_delta_salinity_restore=-0.50.
!  </DATA> 
!
!  <DATA NAME="read_restore_mask" TYPE="logical">
!  For reading in a mask that selects regions of the domain 
!  that are restored (mask=1) or not restored (mask=0).
!  Default  read_restore_mask=.false., whereby restore_mask
!  is set to tmask(k=1). 
!  </DATA> 
!  <DATA NAME="restore_mask_gfdl" TYPE="logical">
!  For modifying the restore mask based on reading in 
!  the GFDL regional mask. Default restore_mask_gfdl=.false.
!  </DATA> 
!  <DATA NAME="salinity_ref" UNITS="psu" TYPE="real">
!  Reference salinity used for converting fresh water flux
!  to salt flux. 
!  </DATA> 
!  <DATA NAME="salt_restore_under_ice" TYPE="logical">
!  Logical indicating whether to restore salinity under sea ice or not.
!  When .false. then will not restore salinity  in regions where we 
!  use a "frazil" condition as a proxy for where sea-ice is present.
!  Do not use sea ice extent from a sea ice model since we generally do 
!  not pass information regarding ice extent between the sea ice model 
!  and the ocean model.     
!  </DATA> 
!  <DATA NAME="zero_net_salt_restore" TYPE="logical">
!  Logical indicating whether to remove the area mean of the salinity 
!  restore flux so there is a net zero input of salt to the ocean
!  associated with restoring.    
!  </DATA> 
!  <DATA NAME="zero_net_salt_correction" TYPE="logical">
!  Logical indicating whether to remove the area mean of the salinity 
!  correction flux so there is a net zero input of salt to the ocean
!  associated with salt correction. 
!  </DATA> 
!  <DATA NAME="zero_net_water_restore" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  restore flux so there is a net zero input of water to the ocean
!  associated with restoring.    
!  </DATA> 
!  <DATA NAME="zero_net_water_correction" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  correction flux so there is a net zero input of water to the ocean
!  associated with water correction.
!  </DATA> 
!  <DATA NAME="zero_net_water_coupler" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  passed through the coupler so there is a net zero input of 
!  fresh water to the ocean associated with p-e+r. Do so by removing 
!  area mean from pme--keep river values unchanged. Note that a choice
!  must be made whether to remove the area mean from rivers or pme.  
!  We choose pme since it is more evenly distributed than rivers.   
!  Also note that we DO NOT include the ice melt in this normalization.
!  </DATA> 
!  <DATA NAME="zero_net_water_couple_restore" TYPE="logical">
!  This logical keeps the total water forcing on the ocean+ice system
!  to a global mean of zero at each time step.  We DO NOT include
!  the ice melt in this normalization.  
!  Setting zero_net_water_couple_restore to true may be appropriate when 
!  running an ice-ocean model using a bulk formulae to compute
!  evaporation (e.g., CORE) and when only providing a weak (or zero)
!  salinity restoring.  It is not appropriate when running a coupled
!  ocean-atmosphere model, where the moisture budget should be 
!  conserved without an artificial removal of the global mean.  
!  </DATA> 
!
!  <DATA NAME="land_model_heat_fluxes" TYPE="logical">
!  For the case where land model passes through the coupler the heat flux 
!  associated with the liquid runoff and calving land ice fields.
!  This heat flux is computed relative to 0C, and takes the form 
!  heat flux = mass flux of water * temp of water * heat capacity, 
!  where the water can be either liquid or solid.  For many coupled models,
!  the water temperature is assumed to be that of the SST.  But 
!  more complete land models now carry the heat of its water relative to 0C,
!  in which case the ocean model does not need to assume anything about the 
!  heat content of the land water. 
!  Default land_model_heat_fluxes=.false.      
!  </DATA> 
!
!  <DATA NAME="debug_water_fluxes" TYPE="logical">
!  Logical for debugging water fluxes. Must be true for any of the 
!  options zero_water_fluxes, zero_calving_fluxes, zero_pme_fluxes
!  or zero_runoff_fluxes to be enabled.  
!  Default debug_water_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_water_fluxes" TYPE="logical">
!  Logical for debugging to zero the pme, river, and pme_taum1 into 
!  ocean, over-riding any input from Ice_ocean_boundary. 
!  Default zero_water_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_calving_fluxes" TYPE="logical">
!  Logical for debugging to zero the calving flux passed into the ocean.
!  Default zero_calving_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_pme_fluxes" TYPE="logical">
!  Logical for debugging to zero the pme flux passed into the ocean.
!  Default zero_pme_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_runoff_fluxes" TYPE="logical">
!  Logical for debugging to zero the runoff flux passed into the ocean.
!  Default zero_runoff_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="zero_river_fluxes" TYPE="logical">
!  Logical for debugging to zero the river (calving+runoff) flux passed into the ocean.
!  Default zero_river_fluxes=.false.      
!  </DATA> 
!  <DATA NAME="convert_river_to_pme" TYPE="logical">
!  Logical for debugging.  Here we add the river water input (calving+runoff)
!  to pme, then set river=calving=runoff=0.0.
!  Default convert_river_to_pme=.false.      
!  </DATA> 
!
!  <DATA NAME="zero_heat_fluxes" TYPE="logical">
!  Logical for debugging to set all heat fluxes into the ocean to zero, 
!  over-riding any input from Ice_ocean_boundary.  Default is .false.      
!  </DATA> 
!  <DATA NAME="zero_surface_stress" TYPE="logical">
!  Logical for debugging to zero all surface stress applied to the ocean,
!  over-riding any input from Ice_ocean_boundary.  Default is .false.      
!  </DATA> 
!
!  <DATA NAME="rotate_winds" TYPE="logical">
!  Set to true when need to rotate the winds onto the ocean model grid.
!  This is needed for cases where the winds are on a spherical grid and 
!  the ocean model uses tripolar=.true.  If generate the wind data on 
!  the ocean model grid, then do not need to rotate, since the rotation 
!  has already been done.  
!  </DATA> 
!
!  <DATA NAME="max_ice_thickness" UNITS="m" TYPE="real">
!  When coupling mom4 to an ice model, the sea ice thickness may need
!  to be restricted to prevent vanishing top-level in mom4. Set 
!  max_ice_thickness (meters) < dzt(k=1) to restrict. This truncation 
!  avoids the numerical problem but we loose mass conservation in the coupled
!  sea ice and ocean system. We also alter the pressure felt on the ocean 
!  as applied by the sea ice. Different vertical coordinates are needed 
!  to do the problem more realistically.   
!
!  Note that the problem of vanishing top layer is removed when use
!  either ZSTAR or PSTAR as vertical coordinate.  
!  </DATA> 
!
!  <DATA NAME="ice_salt_concentration" UNITS="kg salt / kg ice" TYPE="real">
!  The salt concentration of sea ice.  This is taken as a bulk value, and should 
!  be the same as that used by the ice model. Default is ice_salt_concentration=0.005,
!  as that is the value used in the GFDL coupled climate model. 
!  </DATA> 
!
!  <DATA NAME="runoff_salinity" UNITS="g salt / kg runoff water (ppt)" TYPE="real">
!  The salinity of river runoff water. Default is runoff_salinity=0.0.
!  </DATA> 
!  <DATA NAME="runoff_temp_min" UNITS="DegC" TYPE="real">
!  The minimum temperature that river runoff into the ocean is assigned. 
!  Default runoff_temp_min=0.0.
!  </DATA> 
!
!  <DATA NAME="runoffspread" TYPE="logical">
!  Set to true if wish to use the spread_river_horz algorithm to spread 
!  the river runoff flux horizontally over an area into the ocean wider than 
!  set by the coupler.  This option requires the setup of a table for 
!  determining the points over which we spread. 
!  Default runoffspread=.false.
!  </DATA> 
!  <DATA NAME="calvingspread" TYPE="logical">
!  Set to true if wish to use the spread_river_horz algorithm to spread 
!  the calving flux horizontally over an area into the ocean wider than 
!  set by the coupler.  This option requires the setup of a table for 
!  determining the points over which we spread. 
!  Default calvingspread=.false.
!  </DATA> 
!
!  <DATA NAME="avg_sfc_velocity" TYPE="logical">
!  If set to true, the u and v fields passed up to the sea ice
!  are averaged over a coupling interval. TRUE by default.
!  </DATA> 
!  <DATA NAME="avg_sfc_temp_salt_eta" TYPE="logical">
!  If set to true, the t, s and sea_level fields passed up to the sea ice
!  are averaged over a coupling interval. TRUE by default.
!  </DATA> 
!
!  <DATA NAME="use_full_patm_for_sea_level" TYPE="logical">
! The option use_full_patm_for_sea_level allows for the passing 
! of the sea level including the full weight of sea ice back to
! the ice model.  This approach maintains the max weight on the liquid
! ocean according to the nml variable max_ice_thickness.  But it does 
! allow the sea ice to know when there is actually more sea ice than that
! set by max_ice_thickness.  This option then provides for a negative
! feedback on the runaway growth of sea ice, since the full pressure acting to 
! make the ice flow will be correctly felt.  This is a new option, and is not
! fully tested, So the default is use_full_patm_for_sea_level=.false
!  </DATA> 
!
!  <DATA NAME="do_flux_correction" TYPE="logical">
!  For applying surface flux correction to to a tracer or wind stress field. 
!  This code is used at GFDL for idealized perturbation experiments, such 
!  as when one wishes to artificially enhance the wind stress to test 
!  model sensitivity.  It is also appropriate for coupled models that 
!  may require a modification to the fluxes arrising from a coupled model,
!  via reading in information from a pre-defined
!  data file, 
!  Default do_flux_correction=.false.
!  </DATA> 
!  <DATA NAME="temp_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for temperature.  
!  Default temp_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="salt_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for salinity.
!  Default salt_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="tau_x_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for tau_x.
!  Default tau_x_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="tau_y_correction_scale" UNITS="dimensionless" TYPE="real">
!  A scale multiplying the flux correction for tau_y.
!  Default tau_y_correction_scale=0.0.
!  </DATA> 
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency.
!  The default value is do_bitwise_exact_sum=.true. in order to ensure answers
!  do not change when alter processors.  But if wish to enhance the efficiency
!  of coupled ocean-ice models that use one of the global normalization options
!  zero_net_salt_restore        =.true.
!  zero_net_salt_correction     =.true.
!  zero_net_water_restore       =.true.
!  zero_net_water_correction    =.true.
!  zero_net_water_coupler       =.true.
!  zero_net_water_couple_restore=.true.
!  then one may wish to consider setting do_bitwise_exact_sum=.false.
!  </DATA>
!
!</NAMELIST>
!

#include <fms_platform.h>

use constants_mod,            only: epsln, grav, hlv, hlf, kelvin
use diag_manager_mod,         only: register_diag_field, register_static_field, send_data
use fms_mod,                  only: open_namelist_file, check_nml_error, file_exist
use fms_mod,                  only: close_file, read_data, write_version_number
use fms_io_mod,               only: register_restart_field, save_restart, restore_state, restart_file_type
use mpp_domains_mod,          only: mpp_update_domains, BGRID_NE, mpp_define_domains, mpp_get_compute_domain
use mpp_domains_mod,          only: mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_domains_mod,          only: mpp_define_io_domain
use mpp_mod,                  only: mpp_error, FATAL, stdout, stdlog
use time_interp_external_mod, only: time_interp_external, init_external_field
use time_manager_mod,         only: time_type

use ocean_domains_mod,        only: get_local_indices
use ocean_parameters_mod,     only: missing_value
use ocean_parameters_mod,     only: rho_cp, cp_ocean, cp_liquid_runoff, cp_solid_runoff, rho0, rho0r 
use ocean_parameters_mod,     only: CONSERVATIVE_TEMP, POTENTIAL_TEMP 
use ocean_riverspread_mod,    only: spread_river_horz
use ocean_tpm_mod,            only: ocean_tpm_sum_sfc, ocean_tpm_avg_sfc, ocean_tpm_sbc
use ocean_tpm_mod,            only: ocean_tpm_zero_sfc, ocean_tpm_sfc_end
use ocean_types_mod,          only: ocean_grid_type, ocean_domain_type, ocean_public_type
use ocean_types_mod,          only: ocean_time_type, ocean_thickness_type 
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,          only: ocean_external_mode_type, ocean_velocity_type 
use ocean_types_mod,          only: ice_ocean_boundary_type, ocean_density_type
use ocean_types_mod,          only: ocean_public_type
use ocean_workspace_mod,      only: wrk1_2d, wrk2_2d, wrk3_2d
  
implicit none

private

! for restoring input field
integer, allocatable, dimension(:) :: id_restore

! for flux corrections input fields  
integer :: id_tau_x_correction =-1
integer :: id_tau_y_correction =-1
integer, allocatable, dimension(:) :: id_correction

integer :: index_temp         =-1
integer :: index_salt         =-1
integer :: index_diag_temp    =-1
integer :: index_frazil       =-1
integer :: prog_temp_variable =-1

integer :: memuse
integer :: num_prog_tracers
integer :: num_diag_tracers
integer :: global_sum_flag

! ids for diagnostic manager 
logical :: used
integer, allocatable, dimension(:) :: id_stf_coupler
integer, allocatable, dimension(:) :: id_stf_restore
integer, allocatable, dimension(:) :: id_stf_correct
integer, allocatable, dimension(:) :: id_stf_total
integer, allocatable, dimension(:) :: id_stf_runoff
integer, allocatable, dimension(:) :: id_stf_calving
integer, allocatable, dimension(:) :: id_stf_pme
integer, allocatable, dimension(:) :: id_stf_prec
integer, allocatable, dimension(:) :: id_stf_evap
integer, allocatable, dimension(:) :: id_trunoff
integer, allocatable, dimension(:) :: id_tcalving 
integer, allocatable, dimension(:) :: id_triver 

integer, allocatable, dimension(:) :: id_total_ocean_stf_coupler
integer, allocatable, dimension(:) :: id_total_ocean_stf_runoff
integer, allocatable, dimension(:) :: id_total_ocean_stf_calving 
integer, allocatable, dimension(:) :: id_total_ocean_stf_pme
integer, allocatable, dimension(:) :: id_total_ocean_stf_prec
integer, allocatable, dimension(:) :: id_total_ocean_stf_evap
integer, allocatable, dimension(:) :: id_total_ocean_stf_restore
integer, allocatable, dimension(:) :: id_total_ocean_stf_correct
integer, allocatable, dimension(:) :: id_total_ocean_stf_sum

integer :: id_tau_x_flux_correction=-1
integer :: id_tau_y_flux_correction=-1

integer :: id_salt_flux_ice      =-1
integer :: id_total_salt_flux_ice=-1
integer :: id_temp_runoff_eff    =-1
integer :: id_temp_calving_eff   =-1

integer :: id_tau_x          =-1
integer :: id_tau_y          =-1
integer :: id_tau_curl       =-1
integer :: id_ekman_we       =-1
integer :: id_ekman_heat     =-1
integer :: id_swflx          =-1
integer :: id_swflx_vis      =-1

integer :: id_lw_heat            =-1
integer :: id_sens_heat          =-1
integer :: id_fprec_melt_heat    =-1
integer :: id_calving_melt_heat  =-1
integer :: id_evap_heat          =-1

integer :: id_fprec          =-1
integer :: id_lprec          =-1
integer :: id_river          =-1
integer :: id_calving        =-1
integer :: id_runoff         =-1
integer :: id_melt           =-1
integer :: id_evap           =-1
integer :: id_pme_sbc        =-1
integer :: id_pme_river      =-1
integer :: id_pme_restore    =-1
integer :: id_pme_correct    =-1
integer :: id_pme_net        =-1
integer :: id_ice_mask       =-1
integer :: id_open_ocean_mask=-1
integer :: id_restore_mask   =-1

! ids for scalar fields 

integer :: id_total_ocean_swflx             =-1
integer :: id_total_ocean_swflx_vis         =-1
integer :: id_total_ocean_evap_heat         =-1
integer :: id_total_ocean_lw_heat           =-1
integer :: id_total_ocean_sens_heat         =-1
integer :: id_total_ocean_river_heat        =-1
integer :: id_total_ocean_pme_heat          =-1
integer :: id_total_ocean_fprec_melt_heat   =-1
integer :: id_total_ocean_calving_melt_heat =-1

integer :: id_total_ocean_river      =-1
integer :: id_total_ocean_evap       =-1
integer :: id_total_ocean_melt       =-1
integer :: id_total_ocean_pme_sbc    =-1
integer :: id_total_ocean_pme_restore=-1
integer :: id_total_ocean_pme_correct=-1
integer :: id_total_ocean_pme_net    =-1
integer :: id_total_ocean_pme_river  =-1

integer :: id_total_ocean_fprec   =-1
integer :: id_total_ocean_lprec   =-1
integer :: id_total_ocean_calving =-1
integer :: id_total_ocean_runoff  =-1


#include <ocean_memory.h>

#ifdef  MOM4_STATIC_ARRAYS
real, dimension(isd:ied,jsd:jed) :: data
real, dimension(isd:ied,jsd:jed) :: pme_taum1    ! mass flux (kg/(m^2 sec)) of precip-evap from coupler at taum1 time step 
real, dimension(isd:ied,jsd:jed) :: river_taum1  ! mass flux of river water (liquid+solid) from coupler at taum1 time step 
real, dimension(isd:ied,jsd:jed) :: pme_river    ! mass flux of water into ocean from pme+river-melt
real, dimension(isd:ied,jsd:jed) :: restore_mask ! mask for setting regions that are restored 
real, dimension(isd:ied,jsd:jed) :: runoff       ! mass flux of liquid river runoff 
real, dimension(isd:ied,jsd:jed) :: calving      ! mass flux of calving land ice into ocean 
#else
real, allocatable, dimension(:,:) :: data
real, allocatable, dimension(:,:) :: pme_taum1    ! mass flux (kg/(m^2 sec)) of precip-evap from coupler at taum1 time step 
real, allocatable, dimension(:,:) :: river_taum1  ! mass flux of river water (liquid+solid) from coupler at taum1 time step
real, allocatable, dimension(:,:) :: pme_river    ! mass flux of water into ocean from pme+river-melt
real, allocatable, dimension(:,:) :: restore_mask ! mask for setting regions that are restored 
real, allocatable, dimension(:,:) :: runoff       ! mass flux of liquid river runoff  
real, allocatable, dimension(:,:) :: calving      ! mass flux of calving land ice into ocean 
#endif


! ice-ocean-boundary fields are allocated using absolute
! indices (regardless of whether ocean allocations are static)
integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
integer :: i_shift, j_shift                      ! shift isc_bnd to isc and jsc_bnd to jsc

real :: grav_rho0_r
real :: cp_liquid_runoff_r
real :: cp_solid_runoff_r
real :: cp_ocean_r
real :: ice_salt_concentration_r
real :: dtime 

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()
type(restart_file_type), save    :: Sbc_restart
type(restart_file_type), save    :: Sfc_restart

public :: ocean_sbc_init
public :: sum_ocean_sfc
public :: avg_ocean_sfc
public :: zero_ocean_sfc
public :: flux_adjust
public :: get_ocean_sbc
public :: initialize_ocean_sfc
public :: ocean_sfc_end
public :: ocean_sfc_restart

character(len=128) :: version=&
     '$Id: ocean_sbc.F90,v 16.0.2.11.8.1.2.3.10.5.10.1.4.2 2009/12/01 21:41:00 smg Exp $'
character (len=128) :: tagname = &
     '$Name: mom4p1_pubrel_dec2009_nnz $'

logical :: module_is_initialized        =.false.
logical :: use_waterflux                =.true.
logical :: waterflux_tavg               =.false.
logical :: rotate_winds                 =.false.
logical :: runoffspread                 =.false.
logical :: calvingspread                =.false.
logical :: salt_restore_under_ice       =.true.
logical :: salt_restore_as_salt_flux    =.true.
logical :: zero_net_salt_restore        =.false.
logical :: zero_net_salt_correction     =.false.
logical :: zero_net_water_restore       =.false.
logical :: zero_net_water_correction    =.false.
logical :: zero_net_water_coupler       =.false.
logical :: zero_net_water_couple_restore=.false.
logical :: debug_water_fluxes           =.false.
logical :: zero_water_fluxes            =.false. 
logical :: zero_pme_fluxes              =.false. 
logical :: zero_calving_fluxes          =.false. 
logical :: zero_runoff_fluxes           =.false. 
logical :: zero_river_fluxes            =.false. 
logical :: convert_river_to_pme         =.false.
logical :: zero_heat_fluxes             =.false. 
logical :: zero_surface_stress          =.false.
logical :: read_restore_mask            =.false. 
logical :: restore_mask_gfdl            =.false.
logical :: land_model_heat_fluxes       =.false. 
logical :: do_flux_correction           =.false.

real    :: ice_salt_concentration     = 0.005  ! kg/kg
real    :: runoff_salinity            = 0.0    ! psu
real    :: runoff_temp_min            = 0.0    ! degC
real    :: temp_restore_tscale        = -30.
real    :: salt_restore_tscale        = -30.
real    :: max_ice_thickness          = 5.0 
real    :: salinity_ref               = 35.0
real    :: max_delta_salinity_restore = -0.5
real    :: temp_damp_factor
real    :: salt_damp_factor
real    :: temp_correction_scale      = 0.0
real    :: salt_correction_scale      = 0.0
real    :: tau_x_correction_scale     = 0.0
real    :: tau_y_correction_scale     = 0.0

logical :: avg_sfc_velocity           =.true.
logical :: avg_sfc_temp_salt_eta      =.true.
logical :: use_full_patm_for_sea_level=.false. 
logical :: do_bitwise_exact_sum       = .true.

namelist /ocean_sbc_nml/ temp_restore_tscale, salt_restore_tscale, salt_restore_under_ice, salt_restore_as_salt_flux,        &
         rotate_winds, use_waterflux, waterflux_tavg, max_ice_thickness, runoffspread, calvingspread,                        &
         salinity_ref, zero_net_salt_restore, zero_net_water_restore, zero_net_water_coupler, zero_net_water_couple_restore, &
         zero_net_salt_correction, zero_net_water_correction,                                                                &
         debug_water_fluxes, zero_water_fluxes, zero_calving_fluxes, zero_pme_fluxes, zero_runoff_fluxes, zero_river_fluxes, &
         convert_river_to_pme, zero_heat_fluxes, zero_surface_stress, avg_sfc_velocity, avg_sfc_temp_salt_eta,               &
         ice_salt_concentration, runoff_salinity, runoff_temp_min, read_restore_mask, restore_mask_gfdl,                     &
         land_model_heat_fluxes, use_full_patm_for_sea_level, max_delta_salinity_restore, do_flux_correction,                &
         temp_correction_scale, salt_correction_scale, tau_x_correction_scale, tau_y_correction_scale, do_bitwise_exact_sum

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_sbc_init">
!
! <DESCRIPTION>
! Initialize the ocean sbc module 
! </DESCRIPTION>
!
subroutine ocean_sbc_init(Grid, Domain, Time, T_prog, T_diag, Velocity, &
                          Ocean_sfc, time_tendency, dtime_t)

  type(ocean_grid_type),          intent(in),    target :: Grid
  type(ocean_domain_type),        intent(in),    target :: Domain
  type(ocean_time_type),          intent(in)            :: Time
  type(ocean_prog_tracer_type),   intent(inout), target :: T_prog(:)
  type(ocean_diag_tracer_type),   intent(inout), target :: T_diag(:)
  type(ocean_velocity_type),      intent(in),    target :: Velocity
  type(ocean_public_type),        intent(inout)         :: Ocean_sfc
  character(len=32),              intent(in)            :: time_tendency 
  real,                           intent(in)            :: dtime_t

  integer            :: ioun, ierr, io_status
  integer            :: i,j,n
  integer            :: taup1, id_field
  real               :: secday
  character(len=128) :: name, filename

  integer :: stdoutunit,stdlogunit 
  stdoutunit=stdout();stdlogunit=stdlog() 

  if ( module_is_initialized ) then
    call mpp_error(FATAL, '==>Error ocean_sbc_init: module has been initialized')
  endif 
  module_is_initialized = .TRUE.

  dtime = dtime_t 

  call write_version_number(version, tagname)

  ioun = open_namelist_file()
  read  (ioun, ocean_sbc_nml,iostat=io_status)
  write (stdoutunit,'(/)')
  write (stdoutunit, ocean_sbc_nml)  
  write (stdlogunit, ocean_sbc_nml)
  ierr = check_nml_error(io_status,'ocean_sbc_nml')
  call close_file (ioun)

  if(do_bitwise_exact_sum) then
     global_sum_flag = BITWISE_EXACT_SUM
  else
     global_sum_flag = NON_BITWISE_EXACT_SUM
  endif

#ifndef MOM4_STATIC_ARRAYS  
  call get_local_indices(Domain,isd, ied, jsd, jed, isc, iec, jsc, jec)
  allocate(data(isd:ied,jsd:jed))
  allocate(pme_taum1(isd:ied,jsd:jed))
  allocate(river_taum1(isd:ied,jsd:jed))
  allocate(pme_river(isd:ied,jsd:jed))
  allocate(restore_mask(isd:ied,jsd:jed))
  allocate(runoff(isd:ied,jsd:jed))
  allocate(calving(isd:ied,jsd:jed))
#endif
  data              = 0.0
  pme_taum1         = 0.0
  river_taum1       = 0.0
  pme_river         = 0.0
  restore_mask(:,:) = Grid%tmask(:,:,1)
  runoff            = 0.0
  calving           = 0.0

  Dom => Domain
  Grd => Grid

  grav_rho0_r        = rho0r/grav
  secday             = 1.0/(60.0*1440.0)
  taup1              = Time%taup1
  cp_liquid_runoff_r = 1.0/cp_liquid_runoff
  cp_solid_runoff_r  = 1.0/cp_solid_runoff
  cp_ocean_r         = 1.0/cp_ocean 

  if(ice_salt_concentration > 0.0) then
     ice_salt_concentration_r = 1.0/(epsln+ice_salt_concentration) 
  else
     ice_salt_concentration_r = 0.0
  endif 

  call mpp_define_domains((/1,Grd%ni,1,Grd%nj/),Dom%layout,Ocean_sfc%Domain,maskmap=Dom%maskmap, name='sbc', &
                          x_cyclic_offset = Domain%x_cyclic_offset, y_cyclic_offset = Domain%y_cyclic_offset )  
  call mpp_define_io_domain(Ocean_sfc%Domain, Dom%io_layout)
  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)

  i_shift = isc - isc_bnd
  j_shift = jsc - jsc_bnd

  allocate ( Ocean_sfc%t_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%s_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%u_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%v_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%sea_lev(isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%area   (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%frazil (isc_bnd:iec_bnd,jsc_bnd:jec_bnd))

  Ocean_sfc%t_surf  = 0.0  ! time averaged sst (Kelvin) passed to atmosphere/ice model
  Ocean_sfc%s_surf  = 0.0  ! time averaged sss (psu) passed to atmosphere/ice models
  Ocean_sfc%u_surf  = 0.0  ! time averaged u-current (m/sec) passed to atmosphere/ice models
  Ocean_sfc%v_surf  = 0.0  ! time averaged v-current (m/sec)  passed to atmosphere/ice models 
  Ocean_sfc%sea_lev = 0.0  ! time averaged thickness of top model grid cell (m) plus patm/(grav*rho0) 
  Ocean_sfc%frazil  = 0.0  ! time accumulated frazil (J/m^2) passed to ice model

  Ocean_sfc%area    = Grid%dat(isc:iec, jsc:jec) * Grid%tmask(isc:iec, jsc:jec, 1) !grid cell area

  ! set restore_mask=0.0 in those regions where restoring is NOT applied
  ! and restore_mask=1.0 in regions where restoring is applied.  Default is 
  ! to restore everywhere (restore_mask(:,:) = Grd%tmask(:,:,1)).  
  if(read_restore_mask) then 
      call read_data('INPUT/restore_mask','restore_mask',data,Domain%domain2d)
      do j=jsc,jec
         do i=isc,iec
            restore_mask(i,j) = data(i,j)
         enddo
      enddo

      ! the following is specific to the mask used at GFDL
      if(restore_mask_gfdl) then 
          do j=jsc,jec
             do i=isc,iec
                if(restore_mask(i,j)==0.0 .or. restore_mask(i,j)>=6.0) then 
                    restore_mask(i,j)=0.0
                else 
                    restore_mask(i,j)=1.0
                endif
             enddo
          enddo
      endif

      call mpp_update_domains(restore_mask(:,:), Dom%domain2d)
  endif

  num_prog_tracers = size(T_prog)
  num_diag_tracers = size(T_diag)

  ! for file ids 
  allocate( id_restore    (num_prog_tracers) )
  allocate( id_correction (num_prog_tracers) )
  id_restore   (:)  = -1
  id_correction(:)  = -1

  allocate( id_stf_coupler           (num_prog_tracers) )
  allocate( id_stf_restore           (num_prog_tracers) )
  allocate( id_stf_correct           (num_prog_tracers) )
  allocate( id_stf_total             (num_prog_tracers) )
  allocate( id_stf_runoff            (num_prog_tracers) )
  allocate( id_stf_calving           (num_prog_tracers) )
  allocate( id_stf_pme               (num_prog_tracers) )
  allocate( id_stf_prec              (num_prog_tracers) )
  allocate( id_stf_evap              (num_prog_tracers) )
  allocate( id_trunoff               (num_prog_tracers) )
  allocate( id_tcalving              (num_prog_tracers) )
  allocate( id_triver                (num_prog_tracers) )

  allocate( id_total_ocean_stf_coupler(num_prog_tracers) )
  allocate( id_total_ocean_stf_runoff (num_prog_tracers) )
  allocate( id_total_ocean_stf_calving(num_prog_tracers) )
  allocate( id_total_ocean_stf_pme    (num_prog_tracers) )
  allocate( id_total_ocean_stf_prec   (num_prog_tracers) )
  allocate( id_total_ocean_stf_evap   (num_prog_tracers) )
  allocate( id_total_ocean_stf_restore(num_prog_tracers) )
  allocate( id_total_ocean_stf_correct(num_prog_tracers) )
  allocate( id_total_ocean_stf_sum    (num_prog_tracers) )

  id_stf_coupler           (:)  = -1
  id_stf_restore           (:)  = -1
  id_stf_correct           (:)  = -1
  id_stf_total             (:)  = -1
  id_stf_runoff            (:)  = -1
  id_stf_calving           (:)  = -1
  id_stf_pme               (:)  = -1
  id_stf_prec              (:)  = -1
  id_stf_evap              (:)  = -1
  id_trunoff               (:)  = -1
  id_tcalving              (:)  = -1
  id_triver                (:)  = -1

  id_total_ocean_stf_coupler(:) = -1
  id_total_ocean_stf_runoff (:) = -1
  id_total_ocean_stf_calving(:) = -1
  id_total_ocean_stf_pme    (:) = -1
  id_total_ocean_stf_prec   (:) = -1
  id_total_ocean_stf_evap   (:) = -1
  id_total_ocean_stf_restore(:) = -1
  id_total_ocean_stf_correct(:) = -1
  id_total_ocean_stf_sum    (:) = -1

  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  
  do n=1,num_diag_tracers
     if (T_diag(n)%name == 'frazil')   index_frazil    = n
     if (T_diag(n)%name == 'con_temp') index_diag_temp = n
     if (T_diag(n)%name == 'pot_temp') index_diag_temp = n
  enddo

  if (index_temp == -1 .or. index_salt == -1) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_sbc_mod (ocean_sbc_init): temp and/or salt not identified')
  endif 
  if (index_diag_temp == -1) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_sbc_mod (ocean_sbc_init): diagnostic temp not identified')
  endif 
  
  if(T_prog(index_temp)%longname=='Conservative temperature') prog_temp_variable = CONSERVATIVE_TEMP
  if(T_prog(index_temp)%longname=='Potential temperature')    prog_temp_variable = POTENTIAL_TEMP


  ! get file indices for wind stress flux corrections 
  name = 'INPUT/tau_x_correction.nc'
  if (file_exist(trim(name)) .and. do_flux_correction) then
      id_tau_x_correction = init_external_field(name, "tau_x", domain=Dom%domain2d)
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: applying adjust surface i-directed wind stress from tau_x to sfc_momentum fluxes'
      if (id_tau_x_correction == -1) then 
         call mpp_error(FATAL,'==>Error in ocean_sbc_mod: failed to find tau_x field in INPUT/tau_x_correction.nc') 
      endif 
  endif
  name = 'INPUT/tau_y_correction.nc'
  if (file_exist(trim(name)) .and. do_flux_correction ) then
      id_tau_y_correction = init_external_field(name, "tau_y", domain=Dom%domain2d)
      write(stdoutunit,*) &
      '==>Note from ocean_sbc_mod: applying adjust surface j-directed wind stress from tau_y to sfc_momentum fluxes'
      if (id_tau_y_correction == -1) then 
         call mpp_error(FATAL,'==>Error in ocean_sbc_mod: failed to find tau_y field in INPUT/tau_y_correction.nc') 
      endif 
  endif


  do n = 1, num_prog_tracers
     ! init_external_field(file_name,field_name,domain)
#ifndef MOM4_STATIC_ARRAYS         
        allocate(T_prog(n)%stf(isd:ied,jsd:jed))
        allocate(T_prog(n)%tpme(isd:ied,jsd:jed))
        allocate(T_prog(n)%triver(isd:ied,jsd:jed))
        allocate(T_prog(n)%trunoff(isd:ied,jsd:jed))
        allocate(T_prog(n)%tcalving(isd:ied,jsd:jed))
        allocate(T_prog(n)%runoff_tracer_flux(isd:ied,jsd:jed))
        allocate(T_prog(n)%calving_tracer_flux(isd:ied,jsd:jed))
        allocate(T_prog(n)%riverdiffuse(isd:ied,jsd:jed))
#endif

        ! get file indices for restoring fields on temp and salinity  
        name = 'INPUT/'//trim(T_prog(n)%name)//'_sfc_restore.nc'
        if (file_exist(trim(name))) then
            id_restore(n) = init_external_field(name, T_prog(n)%name, domain=Dom%domain2d)
            write(stdoutunit,*) &
            '==>Note from ocean_sbc_mod: applying surface restoring to '//trim(T_prog(n)%name)
            if (id_restore(n) == -1) then
               call mpp_error(FATAL,'==>ocean_sbc_mod: failure to find sfc_restore field in INPUT/_sfc_restore.nc') 
            endif
        elseif(trim(T_prog(n)%name)=='temp' .and. temp_restore_tscale > 0.0) then 
            call mpp_error(FATAL, &
            '==>ocean_sbc_mod: temp_restore_tscale > 0.0 but cannot find INPUT/temp_sfc_restore.nc') 
        elseif(trim(T_prog(n)%name)=='salt' .and. salt_restore_tscale > 0.0) then 
            call mpp_error(FATAL, &
            '==>ocean_sbc_mod: salt_restore_tscale > 0.0 but cannot find INPUT/salt_sfc_restore.nc') 
        endif

        ! get file indices for temp and pme flux correction 
        name = 'INPUT/'//trim(T_prog(n)%name)//'_sfc_correction.nc'
        if (file_exist(trim(name)) .and. do_flux_correction) then
            if  (n == index_temp) then
                id_correction(n) = init_external_field(name, "sfc_hflux", domain=Dom%domain2d)
                write(stdoutunit,*) '==>Note from ocean_sbc_mod: applying surface heat flux correction to '//trim(T_prog(n)%name)
                if (id_correction(n) == -1) then 
                  call mpp_error(FATAL,&
                  '==>Error in ocean_sbc_mod: failure to find temp_sfc_correction field in INPUT/temp_sfc_correction.nc') 
                endif 
            endif
            if  (n == index_salt) then
                id_correction(n) = init_external_field(name, "pme", domain=Dom%domain2d)
                write(stdoutunit,*) '==>Note from ocean_sbc_mod: applying surface pme flux correction to '//trim(T_prog(n)%name)
                if (id_correction(n) == -1) then
                   call mpp_error(FATAL, &
                   '==>Error in ocean_sbc_mod: failure to find salt_sfc_correction field in INPUT/salt_sfc_correction.nc') 
                endif 
            endif
        endif

        T_prog(n)%stf                 = 0.0
        T_prog(n)%tpme                = 0.0
        T_prog(n)%triver              = 0.0
        T_prog(n)%trunoff             = 0.0
        T_prog(n)%tcalving            = 0.0
        T_prog(n)%runoff_tracer_flux  = 0.0
        T_prog(n)%calving_tracer_flux = 0.0
        T_prog(n)%riverdiffuse        = 0.0
        if (n==index_salt) then
           T_prog(n)%trunoff(:,:) = runoff_salinity*Grd%tmask(:,:,1)
        endif 
        write(stdoutunit,*) &
        '==>Note from ocean_sbc_mod: if inputting river water, enable rivermix_mod to get river tracers into ocean.'

        if (n == index_temp) then
            id_stf_coupler(n) = register_diag_field('ocean_model','sfc_hflux_coupler', &
                 Grd%tracer_axes(1:2),                                                 &
                 Time%model_time, 'surface heat flux coming through coupler',          &
                 'Watts/m^2' ,                                                         &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_stf_restore(n) = register_diag_field('ocean_model','sfc_hflux_restore', &
                 Grd%tracer_axes(1:2),                                                 &
                 Time%model_time, 'surface heat flux from restoring',                  &
                 'Watts/m^2' ,                                                         &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            id_stf_correct(n) = register_diag_field('ocean_model','sfc_hflux_correct', &
                 Grd%tracer_axes(1:2),                                                 &
                 Time%model_time, 'surface heat flux from flux correction',            &
                 'Watts/m^2' ,                                                         &
                 missing_value=missing_value,range=(/-1.e4,1.e4/),                     &
                 standard_name='heat_flux_correction')            
            id_stf_total(n) = register_diag_field('ocean_model','sfc_hflux_total', &
                 Grd%tracer_axes(1:2),                                             &
                 Time%model_time, 'total surface heat flux',                       &
                 'Watts/m^2' ,                                                     &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            id_stf_runoff(n) = register_diag_field('ocean_model','sfc_hflux_from_runoff',&
                 Grd%tracer_axes(1:2),                                                   &
                 Time%model_time, 'heat flux (relative to 0C) from liquid river runoff', &    
                 'Watts/m^2' ,                                                           &
                 missing_value=missing_value,range=(/-1.e4,1.e4/),                       &
                 standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')            
            id_stf_calving(n) = register_diag_field('ocean_model','sfc_hflux_from_calving',       &
                 Grd%tracer_axes(1:2),                                                            &
                 Time%model_time, 'heat flux (relative to 0C) from solid land ice entering ocean',&    
                 'Watts/m^2' ,                                                             &
                 missing_value=missing_value,range=(/-1.e4,1.e4/),                         &
                 standard_name='temperature_flux_due_to_icebergs_expressed_as_heat_flux_into_sea_water')            
            id_stf_pme(n) = register_diag_field('ocean_model','sfc_hflux_pme',                                  &
                 Grd%tracer_axes(1:2),                                                                          &
                 Time%model_time, 'heat flux (relative to 0C) from pme transfer of water across ocean surface', &
                 'Watts/m^2' , missing_value=missing_value,range=(/-1.e4,1.e4/))            
            id_stf_prec(n) = register_diag_field('ocean_model','sfc_hflux_from_water_prec',      &
                 Grd%tracer_axes(1:2),                                                           &
                 Time%model_time, 'heat flux from precip transfer of water across ocean surface',&
                 'Watts/m^2' , missing_value=missing_value,range=(/-1.e4,1.e4/),                 &
                 standard_name='temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water')            
            id_stf_evap(n) = register_diag_field('ocean_model','sfc_hflux_from_water_evap',    &
                 Grd%tracer_axes(1:2),                                                         &
                 Time%model_time, 'heat flux from evap transfer of water across ocean surface',&
                 'Watts/m^2' , missing_value=missing_value,range=(/-1.e4,1.e4/),               &
                 standard_name='temperature_flux_due_to_evaporation_expressed_as_heat_flux_into_sea_water')            
            id_trunoff(n) = register_diag_field('ocean_model','temp_runoff',              &
                 Grd%tracer_axes(1:2),                                                    &
                 Time%model_time, 'temperature of liquid river runoff entering the ocean',&
                 'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            
            id_tcalving(n) = register_diag_field('ocean_model','temp_calving', &
                 Grd%tracer_axes(1:2),                                         &
                 Time%model_time, 'temperature of land ice calving into ocean',&
                 'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            
            id_triver(n) = register_diag_field('ocean_model','temp_river',                      &
                 Grd%tracer_axes(1:2),                                                          &
                 Time%model_time, 'temperature of river water (=runoff+calving) entering ocean',&
                 'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            

            id_total_ocean_stf_coupler(n) = register_diag_field('ocean_model','total_ocean_hflux_coupler', &
                 Time%model_time, 'total surface heat flux passed through coupler', 'Watts/1e15' ,         &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_restore(n) = register_diag_field('ocean_model','total_ocean_hflux_restore', &
                 Time%model_time, 'total surface heat flux adjustment from restoring', 'Watts/1e15',       &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_correct(n) = register_diag_field('ocean_model','total_ocean_hflux_correct', &
                 Time%model_time, 'total surface heat flux adjustment from correction', 'Watts/1e15',      &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_sum(n) = register_diag_field('ocean_model','total_ocean_hflux_sum', &
                 Time%model_time, 'total surface heat flux', 'Watts/1e15' ,                        &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_runoff(n) = register_diag_field('ocean_model','total_ocean_runoff_heat',&
                 Time%model_time, 'total heat flux from liquid river runoff', 'Watts/1e15',            &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_calving(n) = register_diag_field('ocean_model','total_ocean_calving_heat',&
                 Time%model_time, 'total heat flux from calving land ice', 'Watts/1e15',                 &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_pme(n) = register_diag_field('ocean_model','total_ocean_hflux_pme', &
                 Time%model_time, 'total heat flux from pme transferring water across surface',    &
                 'Watts/1e15', missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_prec(n) = register_diag_field('ocean_model','total_ocean_hflux_prec', &
                 Time%model_time, 'total heat flux from precip transferring water across surface',   &
                 'Watts/1e15', missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_total_ocean_stf_evap(n) = register_diag_field('ocean_model','total_ocean_hflux_evap', &
                 Time%model_time, 'total heat flux from evap transferring water across surface',     &
                 'Watts/1e15', missing_value=missing_value,range=(/-1.e4,1.e4/))

        elseif(n == index_salt) then 

            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
            id_stf_coupler(n) = register_diag_field('ocean_model',trim(name), &
                 Grd%tracer_axes(1:2),                                        &
                 Time%model_time, trim(name)//': flux from the coupler',      &
                 'kg/(m^2*sec)' ,                                             &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_restore'
            id_stf_restore(n) = register_diag_field('ocean_model',         &
                 trim(name),                                               &
                 Grd%tracer_axes(1:2),                                     &
                 Time%model_time, trim(name)//': flux from restoring term',&
                 'kg/(m^2*sec)' ,                                          &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_correct'
            id_stf_correct(n) = register_diag_field('ocean_model',                &
                 trim(name),                                                      &
                 Grd%tracer_axes(1:2),                                            &
                 Time%model_time, trim(name)//': flux correction from data file', &
                 'kg/(m^2*sec)' ,                                                 &
                 missing_value=missing_value,range=(/-1.e4,1.e4/),                &
                 standard_name='virtual_salt_flux_correction')            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_total'
            id_stf_total(n) = register_diag_field('ocean_model',&
                 trim(name),                                    &
                 Grd%tracer_axes(1:2),                          &
                 Time%model_time, trim(name),                   &
                 'kg/(m^2*sec)' ,                               &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
            id_stf_runoff(n) = register_diag_field('ocean_model',&
                 trim(name),                                     &
                 Grd%tracer_axes(1:2),                           &
                 Time%model_time, trim(name),                    &
                 'kg/(m^2*sec)' ,                                &
                 missing_value=missing_value,range=(/-1.e4,1.e4/), &
                 standard_name='salt_flux_into_sea_water_from_rivers')            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_calving'
            id_stf_calving(n) = register_diag_field('ocean_model',&
                 trim(name),                                      &
                 Grd%tracer_axes(1:2),                            &
                 Time%model_time, trim(name),                     &
                 'kg/(m^2*sec)' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_ice'
            id_salt_flux_ice = register_diag_field('ocean_model',&
                 trim(name),                                     &
                 Grd%tracer_axes(1:2),                           &
                 Time%model_time, trim(name),                    &
                 'kg/(m^2*sec)' ,                                &
                 missing_value=missing_value,range=(/-1.e4,1.e4/), &
                 standard_name='downward_sea_ice_basal_salt_flux')            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_pme'
            id_stf_pme(n) = register_diag_field('ocean_model',&
                 trim(name),                                  &
                 Grd%tracer_axes(1:2),                        &
                 Time%model_time, trim(name),                 &
                 'kg/(m^2*sec)' ,                             &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_prec'
            id_stf_prec(n) = register_diag_field('ocean_model',&
                 trim(name),                                   & 
                 Grd%tracer_axes(1:2),                         &
                 Time%model_time, trim(name),                  &
                 'kg/(m^2*sec)' ,                              &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_evap'
            id_stf_evap(n) = register_diag_field('ocean_model',&
                 trim(name),                                   &
                 Grd%tracer_axes(1:2),                         &
                 Time%model_time, trim(name),                  &
                 'kg/(m^2*sec)' ,                              & 
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = trim(T_prog(n)%name)//'_runoff'
            id_trunoff(n) = register_diag_field('ocean_model',&
                 trim(name),                                  &
                 Grd%tracer_axes(1:2),                        &
                 Time%model_time, trim(name),                 &
                 trim(T_prog(n)%units),                       &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = trim(T_prog(n)%name)//'_calving'
            id_tcalving(n) = register_diag_field('ocean_model',&
                 trim(name),                                   &
                 Grd%tracer_axes(1:2),                         &
                 Time%model_time, trim(name),                  &
                 trim(T_prog(n)%units),                        &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = trim(T_prog(n)%name)//'_river'
            id_triver(n) = register_diag_field('ocean_model',&
                 trim(name),                                 &
                 Grd%tracer_axes(1:2),                       &
                 Time%model_time, trim(name),                &
                 trim(T_prog(n)%units),                      &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            

            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_ice'
            id_total_salt_flux_ice =  register_diag_field('ocean_model',    &
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
            id_total_ocean_stf_coupler(n) =  register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',   &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_restore'
            id_total_ocean_stf_restore(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_correct'
            id_total_ocean_stf_correct(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_sum'
            id_total_ocean_stf_sum(n) = register_diag_field('ocean_model',  &
                 trim(name),Time%model_time, trim(name), 'kg/sec (*1e-15)', &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
            id_total_ocean_stf_runoff(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_calving'
            id_total_ocean_stf_calving(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_pme'
            id_total_ocean_stf_pme(n) = register_diag_field('ocean_model',  &
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_prec'
            id_total_ocean_stf_prec(n) = register_diag_field('ocean_model', &
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_evap'
            id_total_ocean_stf_evap(n) = register_diag_field('ocean_model', &
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            

        else

            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
            id_stf_coupler(n) = register_diag_field('ocean_model',trim(name), &
                 Grd%tracer_axes(1:2),                                        &
                 Time%model_time, trim(name)//': flux from the coupler',      &
                 'kg/(m^2*sec)' ,                                             &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_restore'
            id_stf_restore(n) = register_diag_field('ocean_model',        &
                 trim(name),                                              &
                 Grd%tracer_axes(1:2),                                    &
                 Time%model_time, trim(name)//': flux from the restoring',&
                 'kg/(m^2*sec)' ,                                         &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_correct'
            id_stf_correct(n) = register_diag_field('ocean_model',          &
                 trim(name),                                                &
                 Grd%tracer_axes(1:2),                                      &
                 Time%model_time, trim(name)//': flux from flux correction',&
                 'kg/(m^2*sec)' ,                                           &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_total'
            id_stf_total(n) = register_diag_field('ocean_model',&
                 trim(name),                                    &
                 Grd%tracer_axes(1:2),                          &
                 Time%model_time, trim(name),                   &
                 'kg/(m^2*sec)' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
            id_stf_runoff(n) = register_diag_field('ocean_model',&
                 trim(name),                                     &
                 Grd%tracer_axes(1:2),                           &
                 Time%model_time, trim(name),                    &
                 'kg/(m^2*sec)' ,                                &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_calving'
            id_stf_calving(n) = register_diag_field('ocean_model',&
                 trim(name),                                      &
                 Grd%tracer_axes(1:2),                            &
                 Time%model_time, trim(name),                     &
                 'kg/(m^2*sec)' ,                                 &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_pme'
            id_stf_pme(n) = register_diag_field('ocean_model',&
                 trim(name),                                  &
                 Grd%tracer_axes(1:2),                        &
                 Time%model_time, trim(name),                 &
                 'kg/(m^2*sec)' ,                             &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_prec'
            id_stf_prec(n) = register_diag_field('ocean_model',&
                 trim(name),                                   &
                 Grd%tracer_axes(1:2),                         &
                 Time%model_time, trim(name),                  &
                 'kg/(m^2*sec)' ,                              &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_evap'
            id_stf_evap(n) = register_diag_field('ocean_model',&
                 trim(name),                                   &
                 Grd%tracer_axes(1:2),                         &
                 Time%model_time, trim(name),                  &
                 'kg/(m^2*sec)' ,                              &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = trim(T_prog(n)%name)//'_runoff'
            id_trunoff(n) = register_diag_field('ocean_model',&
                 trim(name),                                  &
                 Grd%tracer_axes(1:2),                        &
                 Time%model_time, trim(name),                 &
                 trim(T_prog(n)%units),                       &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = trim(T_prog(n)%name)//'_calving'
            id_tcalving(n) = register_diag_field('ocean_model',&
                 trim(name),                                   &
                 Grd%tracer_axes(1:2),                         &
                 Time%model_time, trim(name),                  &
                 trim(T_prog(n)%units),                        &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = trim(T_prog(n)%name)//'_river'
            id_triver(n) = register_diag_field('ocean_model',&
                 trim(name),                                 &
                 Grd%tracer_axes(1:2),                       &
                 Time%model_time, trim(name),                &
                 trim(T_prog(n)%units),                      &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            

            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_coupler'
            id_total_ocean_stf_coupler(n) =  register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',   &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_restore'
            id_total_ocean_stf_restore(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_correct'
            id_total_ocean_stf_correct(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_sum'
            id_total_ocean_stf_sum(n) = register_diag_field('ocean_model', &
                 trim(name),Time%model_time, trim(name), 'kg/sec (*1e-15)',&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_runoff'
            id_total_ocean_stf_runoff(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_calving'
            id_total_ocean_stf_calving(n) = register_diag_field('ocean_model',&
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' , &
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_pme'
            id_total_ocean_stf_pme(n) = register_diag_field('ocean_model',   &
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_prec'
            id_total_ocean_stf_prec(n) = register_diag_field('ocean_model', &
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'total_ocean_sfc_'//trim(T_prog(n)%name)//'_flux_evap'
            id_total_ocean_stf_evap(n) = register_diag_field('ocean_model', &
                 trim(name), Time%model_time, trim(name), 'kg/sec (*1e-15)',&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            

        endif
  enddo

  ! dimensions (kg/m^3)*(m/s)
  if (temp_restore_tscale > 0.0) then
     temp_damp_factor = rho0*Grd%dzt(1)*secday/temp_restore_tscale 
  else
     temp_damp_factor = -1.0
     write(stdoutunit,*) &
     '==>Note from ocean_sbc_mod: temp_restore_tscale < 0. no surface restoring for temp'
  endif

  ! dimensions (kg/m^3)*(m/s)
  if (salt_restore_tscale > 0.0) then
     salt_damp_factor = rho0*Grd%dzt(1)*secday/salt_restore_tscale  

     if(salt_restore_under_ice) then 
       write(stdoutunit,*) &
       '==>Note from ocean_sbc_mod: salt_restore_under_ice=.true. => sss restore even under ice.'
     else 
       write(stdoutunit,*) &
       '==>Note from ocean_sbc_mod: salt_restore_under_ice=.false. => no sss restore under ice.'
     endif 

     if(debug_water_fluxes) then 
         if(zero_water_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the moisture fluxes: pme=river=calving=0.0.'
         endif
         if(zero_calving_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the calving fluxes: calving=0.0.'
         endif
         if(zero_runoff_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the river runoff fluxes: runoff=0.0.'
         endif
         if(zero_pme_fluxes) then 
             write(stdoutunit,*) &
                  '==>Warning: Over-riding Ice_ocean_boundary to zero the pme fluxes: pme=0.0.'
         endif
     endif

     if(use_waterflux) then 
       if(zero_net_water_restore) then 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_restore=.true.=>zero net restoring water put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_restore=.false.=>nonzero net restoring water put in ocean.'
       endif 
       if(zero_net_water_correction) then 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_correction=.true.=>zero net correction water put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_correction=.false.=>nonzero net correction water put in ocean.'
       endif 
       if(zero_net_water_coupler) then 
          write(stdoutunit,*) &
           '==>Note from ocean_sbc_mod: zero_net_water_coupler=.true.=>zero water into ocean via coupler (sans the sea ice).'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_coupler=.false.=>nonzero net water into ocean via coupler.'
       endif 
       if(zero_net_water_couple_restore) then 
          write(stdoutunit,*) &
          '==>ocean_sbc_mod: zero_net_water_couple_restore=.true.=>zero water into ocean from restore + coupler (sans the sea ice).'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_water_couple_restore=.false.'
       endif 
     else 
       if(zero_net_salt_restore) then 
          write(stdoutunit,*) &
         '==>Note from ocean_sbc_mod: zero_net_salt_restore=.true.=>zero net restoring salt put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_salt_restore=.false.=>nonzero net restoring salt put in ocean.'
       endif 
       if(zero_net_salt_correction) then 
          write(stdoutunit,*) &
         '==>Note from ocean_sbc_mod: zero_net_salt_correction=.true.=>zero net correction of salt put in ocean.'
       else 
          write(stdoutunit,*) &
          '==>Note from ocean_sbc_mod: zero_net_salt_correction=.false.=>nonzero net correction salt put in ocean.'
       endif 
     endif 

  else

     salt_damp_factor = -1.0
     write(stdoutunit,*) &
     '==>Note from ocean_sbc_mod: salt_restore_tscale < 0. no surface restoring for salt.'

 endif

 if(use_full_patm_for_sea_level) then 
     write(stdoutunit,*) &
     '==>NOTE: Allowing for the ice model to feel the fully depressed sea level for purposes of its dynamics.'
 endif 

 if(land_model_heat_fluxes) then 
     write(stdoutunit,*) &
     '==>NOTE: Assuming the land model carries heat of liquid runoff and solid calving.'
     write(stdoutunit,*) &
     '   Be sure to make the appropriate changes ALSO in ocean_rivermix.F90 and ocean_vert_kpp.F90'
 endif 

 if(zero_heat_fluxes) then 
     write(stdoutunit,*) &
     '==>Warning: Over-riding Ice_ocean_boundary to zero the heat fluxes: stf(temp)=0.0.'
 endif

 if(zero_surface_stress) then 
     write(stdoutunit,*) &
     '==>Warning: Over-riding Ice_ocean_boundary to zero the surface stress: smf=0.0.'
 endif

 if(avg_sfc_velocity) then 
    write(stdoutunit,*) &
    '==>If coupling, then avg_sfc_velocity=.true. means will pass averaged ocean velocity to ice model.'
 else 
    write(stdoutunit,*) &
    '==>If coupling, then avg_sfc_velocity=.false. means will pass most recent ocean velocity to ice model.'
 endif 
 if(avg_sfc_temp_salt_eta) then 
    write(stdoutunit,*) &
    '==>If coupling, then avg_sfc_temp_salt_eta=.true. means will pass averaged sst, sss, eta to ice model.'
 else 
    write(stdoutunit,*) &
    '==>If coupling, then avg_sfc_velocity=.false. means will pass most recent sst, sss, eta to ice model.'
 endif 


 id_restore_mask = register_static_field('ocean_model','restore_mask', &
                   Grd%tracer_axes(1:2),'restoring mask', 'none' ,     &
                   missing_value=missing_value,range=(/-10.,10./))
 if (id_restore_mask > 0) then
     used = send_data(id_restore_mask, restore_mask(:,:), Time%model_time, &
            rmask=Grd%tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
 endif

 id_temp_runoff_eff = register_diag_field('ocean_model','temp_runoff_eff',        &
      Grd%tracer_axes(1:2),                                                       &
      Time%model_time, 'effective temp of liquid river runoff entering the ocean',&
      'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            

 id_temp_calving_eff = register_diag_field('ocean_model','temp_calving_eff',&
      Grd%tracer_axes(1:2),                                                 &
      Time%model_time, 'effective temp of land ice calving into ocean',     & 
      'degC' , missing_value=missing_value,range=(/-1.e4,1.e4/))            

 id_tau_x = register_diag_field('ocean_model','tau_x', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'i-directed wind stress', 'N/m^2',                  &
       missing_value=missing_value,range=(/-10.,10./),                      &
       standard_name='surface_downward_x_stress')

 id_tau_x_flux_correction = register_diag_field('ocean_model','tau_x_flux_correction',&
       Grd%vel_axes_uv(1:2),                                                          &
       Time%model_time, 'i-directed wind stress flux corretion', 'N/m^2',             &
       missing_value=missing_value,range=(/-10.,10./),                                &
       standard_name='surface_downward_x_stress_correction')

 id_tau_y = register_diag_field('ocean_model','tau_y', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'j-directed wind stress', 'N/m^2' ,                 &
       missing_value=missing_value,range=(/-10.,10./),                      &
       standard_name='surface_downward_y_stress') 

 id_tau_y_flux_correction = register_diag_field('ocean_model','tau_x_flux_correction',&
       Grd%vel_axes_uv(1:2),                                                          &
       Time%model_time, 'j-directed wind stress flux corretion', 'N/m^2',             &
       missing_value=missing_value,range=(/-10.,10./),                                &
       standard_name='surface_downward_y_stress_correction')

 id_tau_curl = register_diag_field('ocean_model','tau_curl', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'wind stress curl averaged to U-point', 'N/m^3',          &
       missing_value=missing_value,range=(/-10.,10./)) 

 id_ekman_we = register_diag_field('ocean_model','ekman_we', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'Ekman vertical velocity averaged to U-point', 'm/s',     &
       missing_value=missing_value,range=(/-100.,100./)) 

 id_ekman_heat = register_diag_field('ocean_model','ekman_heat',                            &
                 Grd%tracer_axes(1:2), Time%model_time, 'Ekman Component to heat transport',&
                 'Watts', missing_value=missing_value,range=(/-1.e4,1.e4/))           

 id_ice_mask = register_diag_field('ocean_model','ice_mask', Grd%tracer_axes(1:2),&
       Time%model_time, 'ice mask according to near-frazil condition', 'none' ,   &
       missing_value=missing_value,range=(/-10.,10./))

 id_open_ocean_mask = register_diag_field('ocean_model','open_ocean_mask', Grd%tracer_axes(1:2),&
       Time%model_time, 'open-ocean mask according to near-frazil condition', 'none' ,          &
       missing_value=missing_value,range=(/-10.,10./))

 id_river = register_diag_field('ocean_model','river', Grd%tracer_axes(1:2),    &
       Time%model_time, 'mass flux of river (runoff + calving) entering ocean', &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))

 id_runoff = register_diag_field('ocean_model','runoff', Grd%tracer_axes(1:2),&
       Time%model_time, 'mass flux of liquid river runoff entering ocean ',   &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),    &
       standard_name='water_flux_into_sea_water_from_rivers')

 id_calving = register_diag_field('ocean_model','ice_calving', Grd%tracer_axes(1:2),&
       Time%model_time, 'mass flux of land ice calving into ocean',                 &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),          &
       standard_name='water_flux_into_sea_water_from_icebergs')

 id_evap = register_diag_field('ocean_model','evap', Grd%tracer_axes(1:2),&
       Time%model_time, 'evaporative mass flux (>0 leaves ocean)',        &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),&
       standard_name='water_evaporation_flux' )

 id_melt = register_diag_field('ocean_model','melt', Grd%tracer_axes(1:2),   &
       Time%model_time, 'liquid water melted from sea ice (>0 enters ocean)',&
       '(kg/m^3)*(m/sec)',missing_value=missing_value,range=(/-1e6,1e6/),    &
       standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics')

 id_pme_river= register_diag_field('ocean_model','pme_river', Grd%tracer_axes(1:2),            &
       Time%model_time, 'mass flux of precip-evap+river via sbc (liquid, frozen, evaporation)',&
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/),                     &
       standard_name='water_flux_into_sea_water')

 id_pme_sbc = register_diag_field('ocean_model','pme_sbc', Grd%tracer_axes(1:2), &
       Time%model_time, 'precip-evap via sbc (liquid, frozen, evaporation)',     &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))

 id_pme_restore = register_diag_field('ocean_model','pme_restore', Grd%tracer_axes(1:2), &
       Time%model_time, 'precip-evap from restore (>0 enters ocean)',                    &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))   

 id_pme_correct = register_diag_field('ocean_model','pme_correct', Grd%tracer_axes(1:2), &
       Time%model_time, 'precip-evap from flux correction (>0 enters ocean)',            &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))   

 id_pme_net = register_diag_field('ocean_model','pme_net', Grd%tracer_axes(1:2), &
       Time%model_time, 'precip-evap into ocean (total w/ restore + normalize)', &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1e6,1e6/))   

 id_fprec = register_diag_field('ocean_model','fprec', Grd%tracer_axes(1:2),   &
       Time%model_time, 'snow falling onto ocean (>0 enters ocean)',           &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1.e10,1.e10/), &
       standard_name='snowfall_flux')   

 id_lprec = register_diag_field('ocean_model','lprec', Grd%tracer_axes(1:2),  &
       Time%model_time, 'liquid precip into ocean (>0 enters ocean)',         &
       '(kg/m^3)*(m/sec)', missing_value=missing_value,range=(/-1.e10,1.e10/),&
       standard_name='rainfall_flux')   

 id_swflx = register_diag_field('ocean_model','swflx', Grd%tracer_axes(1:2),   &
       Time%model_time, 'shortwave flux into ocean (>0 heats ocean)', 'W/m^2', &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                     &
       standard_name='surface_net_downward_shortwave_flux')   

 id_swflx_vis = register_diag_field('ocean_model','swflx_vis', Grd%tracer_axes(1:2),&
       Time%model_time, 'visible shortwave into ocean (>0 heats ocean)', 'W/m^2' ,  &
       missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_evap_heat = register_diag_field('ocean_model','evap_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'latent heat flux into ocean (<0 cools ocean)', 'W/m^2',    &
       missing_value=missing_value,range=(/-1e10,1e10/),                            &
       standard_name='surface_downward_latent_heat_flux')

 id_lw_heat = register_diag_field('ocean_model','lw_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'longwave flux into ocean (<0 cools ocean)', 'W/m^2' ,  &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                      &
       standard_name='surface_net_downward_longwave_flux' )   

 id_sens_heat = register_diag_field('ocean_model','sens_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'sensible heat into ocean (<0 cools ocean)', 'W/m^2' ,      &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                          &
       standard_name='surface_downward_sensible_heat_flux')   

 id_fprec_melt_heat = register_diag_field('ocean_model','fprec_melt_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'heat flux to melt frozen precip (<0 cools ocean)', 'W/m^2' , &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                            &
       standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics')   

 id_calving_melt_heat = register_diag_field('ocean_model','calving_melt_heat', Grd%tracer_axes(1:2),&
       Time%model_time, 'heat flux needed to melt calving ice (<0 cools ocean)', 'W/m^2',           &
       missing_value=missing_value,range=(/-1.e10,1.e10/),                                          &
       standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics')   

 id_total_ocean_river = register_diag_field('ocean_model','total_ocean_river',                  &
       Time%model_time, 'total liquid river water and calving ice entering ocean', 'kg/sec/1e15',&
       missing_value=missing_value,range=(/-1e6,1e6/))

 id_total_ocean_evap = register_diag_field('ocean_model','total_ocean_evap',   &
       Time%model_time, 'total evaporative ocean mass flux (>0 leaves ocean)', &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1e6,1e6/))

 id_total_ocean_melt = register_diag_field('ocean_model','total_ocean_melt',        &
       Time%model_time, 'total liquid water melted from sea ice (>0 enters ocean)', &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1e6,1e6/))

 id_total_ocean_pme_river = register_diag_field('ocean_model','total_ocean_pme_river',         &
       Time%model_time, 'total ocean precip-evap+river via sbc (liquid, frozen, evaporation)', &
       'kg/sec/1e15' , missing_value=missing_value,range=(/-1e6,1e6/))

 id_total_ocean_pme_sbc = register_diag_field('ocean_model','total_ocean_pme_sbc',       &
       Time%model_time, 'total ocean precip-evap via sbc (liquid, frozen, evaporation)', &
       'kg/sec/1e15' , missing_value=missing_value,range=(/-1e6,1e6/))

 id_total_ocean_pme_restore = register_diag_field('ocean_model','total_ocean_pme_restore', &
       Time%model_time, 'total precip-evap from restore (>0 enters ocean)',                &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1e6,1e6/))   

 id_total_ocean_pme_correct = register_diag_field('ocean_model','total_ocean_pme_correct', &
       Time%model_time, 'total precip-evap from flux correction (>0 enters ocean)',        &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1e6,1e6/))   

 id_total_ocean_pme_net = register_diag_field('ocean_model','total_ocean_pme_net',     &
       Time%model_time, 'total precip-evap into ocean (total w/ restore + normalize)', &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1e6,1e6/))   

 id_total_ocean_fprec = register_diag_field('ocean_model','total_ocean_fprec',  &
       Time%model_time, 'total snow falling onto ocean (>0 enters ocean)',      &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_lprec = register_diag_field('ocean_model','total_ocean_lprec',  &
       Time%model_time, 'total liquid precip into ocean (>0 enters ocean)',     &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_runoff = register_diag_field('ocean_model','total_ocean_runoff',&
       Time%model_time, 'total liquid river runoff (>0 water enters ocean)',    &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_calving = register_diag_field('ocean_model','total_ocean_calving', &
       Time%model_time, 'total water entering ocean from calving land ice',        &
       'kg/sec/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_swflx = register_diag_field('ocean_model','total_ocean_swflx', &
       Time%model_time, 'total shortwave flux into ocean (>0 heats ocean)',    &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_swflx_vis = register_diag_field('ocean_model','total_ocean_swflx_vis', &
       Time%model_time, 'total visible shortwave into ocean (>0 heats ocean)',         &
       'Watts/1e15' , missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_evap_heat = register_diag_field('ocean_model','total_ocean_evap_heat', &
       Time%model_time, 'total latent heat flux into ocean (<0 cools ocean)',          &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))

 id_total_ocean_lw_heat = register_diag_field('ocean_model','total_ocean_lw_heat',  &
       Time%model_time, 'total longwave flux into ocean (<0 cools ocean)',          &
       'Watts/1e15',missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_sens_heat = register_diag_field('ocean_model','total_ocean_sens_heat', &
       Time%model_time, 'total sensible heat into ocean (<0 cools ocean)',             &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_fprec_melt_heat = register_diag_field('ocean_model','total_ocean_fprec_melt_heat',&
       Time%model_time, 'total heat flux to melt frozen precip (<0 cools ocean)',                 &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_calving_melt_heat = register_diag_field('ocean_model','total_ocean_calving_melt_heat',&
       Time%model_time, 'total heat flux to melt frozen land ice (<0 cools ocean)',                   &
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   

 id_total_ocean_river_heat = register_diag_field('ocean_model','total_ocean_river_heat',       &
       Time%model_time, 'total heat flux into ocean from liquid+solid runoff (<0 cools ocean)',&
       'Watts/1e15', missing_value=missing_value,range=(/-1.e10,1.e10/))   


  if(waterflux_tavg) then 
    write(stdoutunit,*) &
    '==>Note: waterflux_tavg sets pme+river = avg of ice_ocean_boundary values over'
    write(stdoutunit,*) &
    '         tau and taum1 time steps. This may damp splitting between leap-frog modes.'
    write(stdoutunit,*) &
    '         However, it compromises the conservation of mass between component models.'

    if(time_tendency=='twolevel') then 
      write(stdoutunit,'(/a)') &
      '==>Warning in ocean_sbc_mod: waterflux_tavg=.true. unnecessary with time_tendency==twolevel.'
      write(stdoutunit,'(/a)') &
      '   Strongly recommend setting waterflux_tavg=.false. to conserve mass between component models.'
    endif 

    filename = 'ocean_waterflux.res.nc'
    id_field = register_restart_field(Sbc_restart, filename, 'pme_taum1',   pme_taum1,  Domain%domain2d)
    id_field = register_restart_field(Sbc_restart, filename, 'river_taum1', river_taum1,Domain%domain2d)
    if (file_exist('INPUT/ocean_waterflux.res.nc')) then
       call restore_state(Sbc_restart)
    endif

  endif 
  
  return

end subroutine ocean_sbc_init
! </SUBROUTINE> NAME="ocean_sbc_init"


 !#######################################################################
! <SUBROUTINE NAME="initialize_ocean_sfc">
!
! <DESCRIPTION>
! Initialize the ocean surface type, which passes information between ocean 
! and other component models. 
!
! Note that ocean model sst passed to the atmosphere must be the surface
! potential temperature (which is equated to surface in situ temperature).
! If the ocean prognostic temperature variable is conservative temperature,
! then the sst is carried in T_diag(index_diag_temp).  If the prognostic 
! temperature is potential temperature, then the sst is carried in 
! T_prog(index_temp). 
!
!  Ocean_sfc%t_surf  = time averaged sst (Kelvin) passed to atmosphere/ice model
!  Ocean_sfc%s_surf  = time averaged sss (psu) passed to atmosphere/ice models
!  Ocean_sfc%u_surf  = time averaged u-current (m/sec) passed to atmosphere/ice models
!  Ocean_sfc%v_surf  = time averaged v-current (m/sec)  passed to atmosphere/ice models 
!  Ocean_sfc%sea_lev = time averaged ocean free surface height (m) plus patm/(grav*rho0) 
!  Ocean_sfc%frazil  = time accumulated frazil (J/m^2) passed to ice model.  time averaging 
!                      not performed, since ice model needs the frazil accumulated over the 
!                      ocean time steps.  Note that Ocean_sfc%frazil is accumulated, whereas 
!                      T_diag%frazil (saved in diagnostic tracer restart file) is instantaneous. 
!
! </DESCRIPTION>
!
subroutine initialize_ocean_sfc(Time, Thickness, T_prog, T_diag, Velocity, Ocean_sfc)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_prog_tracer_type),  intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),  intent(in)    :: T_diag(:)
  type(ocean_velocity_type),     intent(in)    :: Velocity
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc

  integer                          :: taup1, id_field
  real, dimension(isd:ied,jsd:jed) :: sst
  character(len=128)               :: filename

  taup1 = Time%taup1
 
  Ocean_sfc%avg_kount = 0
  Ocean_sfc%t_surf    = kelvin
  sst                 = 0.0

  if(prog_temp_variable==CONSERVATIVE_TEMP) then
    sst(:,:) = T_diag(index_diag_temp)%field(:,:,1)
  else
    sst(:,:) = T_prog(index_temp)%field(:,:,1,taup1)
  endif 

  where (Grd%tmask(isc:iec,jsc:jec,1) == 1.0)
    Ocean_sfc%t_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = sst(isc:iec,jsc:jec) + kelvin
    Ocean_sfc%s_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = T_prog(index_salt)%field(isc:iec,jsc:jec,1,taup1)
    Ocean_sfc%u_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = Velocity%u(isc:iec,jsc:jec,1,1,taup1)
    Ocean_sfc%v_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = Velocity%u(isc:iec,jsc:jec,1,2,taup1)
    Ocean_sfc%sea_lev(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)= Thickness%sea_lev(isc:iec,jsc:jec)
    Ocean_sfc%frazil(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = 0.0
  end where 

  filename = 'ocean_sbc.res.nc'
  id_field = register_restart_field(Sfc_restart, filename, 't_surf', Ocean_sfc%t_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 's_surf', Ocean_sfc%s_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'u_surf', Ocean_sfc%u_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'v_surf', Ocean_sfc%v_surf,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'sea_lev',Ocean_sfc%sea_lev,Ocean_sfc%Domain)
  id_field = register_restart_field(Sfc_restart, filename, 'frazil', Ocean_sfc%frazil,Ocean_sfc%Domain)
  if (file_exist('INPUT/ocean_sbc.res.nc')) then
     call restore_state(Sfc_restart)
  endif

  return

end subroutine initialize_ocean_sfc
! </SUBROUTINE> NAME="initialize_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="sum_ocean_sfc">
!
! <DESCRIPTION>
! Accumulate the ocean_sfc derived type over the course of the 
! ocean component sub-cycling used when coupling to other models. 
!
! Note that ocean model sst passed to the atmosphere must be the surface
! potential temperature (which is equated to surface in situ temperature).
! If the ocean prognostic temperature variable is conservative temperature,
! then the sst is carried in T_diag(index_diag_temp).  If the prognostic 
! temperature is potential temperature, then the sst is carried in 
! T_prog(index_temp). 
!
! Note that this routine is called after eta_and_pbot_diagnose,
! so Thickness%eta is eta_t(taup1).  
!
! </DESCRIPTION>
!  
subroutine sum_ocean_sfc(Time, Thickness, T_prog, T_diag, Dens, Velocity, Ocean_sfc)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_prog_tracer_type),  intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),  intent(in)    :: T_diag(:)
  type(ocean_density_type),      intent(in)    :: Dens
  type(ocean_velocity_type),     intent(in)    :: Velocity
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc
  
  real, dimension(isd:ied,jsd:jed) :: sst
  integer                          :: taup1, i, j, ii, jj

  if (Ocean_sfc%avg_kount == 0) call zero_ocean_sfc(Ocean_sfc)

  taup1 = Time%taup1

  Ocean_sfc%avg_kount = Ocean_sfc%avg_kount + 1

  if(prog_temp_variable==CONSERVATIVE_TEMP) then
    sst(:,:) = T_diag(index_diag_temp)%field(:,:,1)
  else
    sst(:,:) = T_prog(index_temp)%field(:,:,1,taup1)
  endif 

  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j) + sst(ii,jj) 
        Ocean_sfc%s_surf(i,j)  = Ocean_sfc%s_surf(i,j) + T_prog(index_salt)%field(ii,jj,1,taup1)
        Ocean_sfc%u_surf(i,j)  = Ocean_sfc%u_surf(i,j) + Velocity%u(ii,jj,1,1,taup1)
        Ocean_sfc%v_surf(i,j)  = Ocean_sfc%v_surf(i,j) + Velocity%u(ii,jj,1,2,taup1)
        Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j)+ Thickness%sea_lev(ii,jj)
     enddo
  enddo

  if(index_frazil > 0) then 
     do j = jsc_bnd,jec_bnd
        do i = isc_bnd,iec_bnd
           ii = i + i_shift
           jj = j + j_shift
           Ocean_sfc%frazil(i,j) = Ocean_sfc%frazil(i,j) + T_diag(index_frazil)%field(ii,jj,1)
        enddo
     enddo
  endif 

  call ocean_tpm_sum_sfc(Dom, T_prog(:), Dens, Ocean_sfc, Time, Grd, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)
      
end subroutine sum_ocean_sfc
! </SUBROUTINE> NAME="sum_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="zero_ocean_sfc">
!
! <DESCRIPTION>
! Zero the elements of the Ocean_sfc derived type.  
! </DESCRIPTION>
! 
subroutine zero_ocean_sfc(Ocean_sfc)

  type(ocean_public_type), intent(inout), target :: Ocean_sfc

  integer :: i, j

  Ocean_sfc%avg_kount = 0

  do j = jsc_bnd, jec_bnd
     do i = isc_bnd, iec_bnd
        Ocean_sfc%t_surf(i,j) = 0.0
        Ocean_sfc%s_surf(i,j) = 0.0
        Ocean_sfc%u_surf(i,j) = 0.0
        Ocean_sfc%v_surf(i,j) = 0.0
        Ocean_sfc%sea_lev(i,j)= 0.0
        Ocean_sfc%frazil(i,j) = 0.0
     enddo
  enddo

  call ocean_tpm_zero_sfc(Ocean_sfc)

end subroutine zero_ocean_sfc
! </SUBROUTINE> NAME="zero_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="avg_ocean_sfc">
!
! <DESCRIPTION>
! Compute average of ocean surface quantities.  This is for coupling, 
! where pass time averaged information from ocean to other component
! models. Note that Ocean_sfc%frazil is NOT time averaged.  Rather, it 
! is accumulated from T_diag(index_frazil)%field in subroutine sum_ocean_sfc.
! Doing so is necessary for heat conservation between ocean and sea 
! ice systems.  Since it is not time averaged, frazil is not part of 
! this averaging subroutine.  
!
! Note that ocean model SST passed to the atmosphere is the surface
! potential temperature (which is equal to surface in situ temperature).
! If the ocean prognostic temperature variable is conservative temperature,
! then the sst is carried in T_diag(index_diag_temp).  If the prognostic 
! temperature is potential temperature, then the sst is carried in 
! T_prog(index_temp). 
!
! Note that if one removes the averaging, then we take only the 
! latest values of the surface fields.  This approach has been 
! found useful to stabilize the "concurrent" coupling approach.  
!
! Note that this routine is called after eta_and_pbot_diagnose,
! so Thickness%eta is eta_t(taup1).  
!
! </DESCRIPTION>
!
subroutine avg_ocean_sfc(Time, Thickness, T_prog, T_diag, Velocity, Ocean_sfc)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_thickness_type),    intent(in)    :: Thickness
  type(ocean_prog_tracer_type),  intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type),  intent(in)    :: T_diag(:)
  type(ocean_velocity_type),     intent(in)    :: Velocity
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc
  
  real, dimension(isd:ied,jsd:jed) :: sst
  real                             :: divid
  integer                          :: taup1, i, j, ii, jj

  taup1 = Time%taup1

  if ( Ocean_sfc%avg_kount == 0) then 
    call mpp_error (FATAL,&
    '==>Error from ocean_sbc_mod (avg_ocean_sfc): no ocean surface quantities have been time averaged')
  endif 

  divid = 1./float(Ocean_sfc%avg_kount)

  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        if(Grd%tmask(ii,jj,1) == 1.0) then
           Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j)*divid + kelvin  !C --> K
           Ocean_sfc%s_surf(i,j)  = Ocean_sfc%s_surf(i,j)*divid
           Ocean_sfc%u_surf(i,j)  = Ocean_sfc%u_surf(i,j)*divid
           Ocean_sfc%v_surf(i,j)  = Ocean_sfc%v_surf(i,j)*divid
           Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j)*divid 
        endif
     enddo
  enddo


  !replace time-averaged t,s,sealev with latest value
  if(.NOT. avg_sfc_temp_salt_eta) then 

      if(prog_temp_variable==CONSERVATIVE_TEMP) then
          sst(:,:) = T_diag(index_diag_temp)%field(:,:,1)
      else
          sst(:,:) = T_prog(index_temp)%field(:,:,1,taup1)
      endif

      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            if(Grd%tmask(ii,jj,1) == 1.0) then
                Ocean_sfc%t_surf(i,j) = sst(ii,jj) + kelvin
                Ocean_sfc%s_surf(i,j) = T_prog(index_salt)%field(ii,jj,1,taup1)
                Ocean_sfc%sea_lev(i,j)= Thickness%sea_lev(ii,jj)
            endif
         enddo
      enddo

  end if


  !replace time-averaged u,v with latest value
  if(.NOT. avg_sfc_velocity) then 
     do j = jsc_bnd,jec_bnd
        do i = isc_bnd,iec_bnd
           ii = i + i_shift
           jj = j + j_shift
           if(Grd%tmask(ii,jj,1) == 1.0) then
              Ocean_sfc%u_surf(i,j) = Velocity%u(ii,jj,1,1,taup1)
              Ocean_sfc%v_surf(i,j) = Velocity%u(ii,jj,1,2,taup1)
           endif
        enddo
     enddo
  end if

  call ocean_tpm_avg_sfc(Dom, Ocean_sfc, Grd, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)

  !set count to zero (surface quantities will be zeroed out before next sum)
  Ocean_sfc%avg_kount = 0   

  
end subroutine avg_ocean_sfc
! </SUBROUTINE> NAME="avg_ocean_sfc"

!#######################################################################
! <SUBROUTINE NAME="ocean_sfc_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_sfc_restart(time_stamp)
  character(len=*),           intent(in), optional :: time_stamp

  call save_restart(Sfc_restart, time_stamp)
  if(waterflux_tavg) call save_restart(sbc_restart, time_stamp)

end subroutine ocean_sfc_restart
! </SUBROUTINE> NAME="ocean_sbc_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_sfc_end">
!
! <DESCRIPTION>
! Save information from Ocean_sfc to restarts. Note that it is 
! important in general to distinguish the time accumulated quantity 
! Ocean_sfc%frazil, saved here, from the instantaneous quantity 
! T_diag%frazil, which is saved in the diagnostic tracer restart file.  
! </DESCRIPTION>
!
subroutine ocean_sfc_end(Ocean_sfc)
  type(ocean_public_type), intent(in), target :: Ocean_sfc

    call ocean_sfc_restart

    call ocean_tpm_sfc_end

end subroutine ocean_sfc_end
! </SUBROUTINE> NAME="ocean_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="get_ocean_sbc">
!
! <DESCRIPTION>
! Subroutine to get the surface fluxes passed into the ocean from 
! other component models. 
! </DESCRIPTION>
!
subroutine get_ocean_sbc(Time, Ice_ocean_boundary, Thickness, Ext_mode, T_prog, Velocity, &
                         pme, melt, river, runoff, calving, upme, uriver, swflx, swflx_vis, patm)

  type(ocean_time_type),          intent(in)    :: Time 
  type(ice_ocean_boundary_type),  intent(in)    :: Ice_ocean_boundary
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  type(ocean_prog_tracer_type),   intent(inout) :: T_prog(:)
  type(ocean_velocity_type),      intent(inout) :: Velocity
  real, dimension(isd:,jsd:),     intent(inout) :: pme
  real, dimension(isd:,jsd:),     intent(inout) :: melt
  real, dimension(isd:,jsd:),     intent(inout) :: river
  real, dimension(isd:,jsd:),     intent(inout) :: runoff
  real, dimension(isd:,jsd:),     intent(inout) :: calving 
  real, dimension(isd:,jsd:),     intent(inout) :: swflx
  real, dimension(isd:,jsd:),     intent(inout) :: swflx_vis
  real, dimension(isd:,jsd:),     intent(inout) :: patm
  real, dimension(isd:,jsd:,:),   intent(inout) :: upme
  real, dimension(isd:,jsd:,:),   intent(inout) :: uriver

  real, dimension(isd:ied,jsd:jed) :: tmp_flux
  real, dimension(isd:ied,jsd:jed) :: tmp_patm

  real, dimension(isd:ied,jsd:jed) :: liquid_precip
  real, dimension(isd:ied,jsd:jed) :: evaporation


  real    :: tmp_x, tmp_y
  real    :: totz, tbrz, maskt, factor 
  real    :: pme_river_total, var
  real    :: total_stuff 
  real    :: umask_norm, tracer_input
  real    :: tmp_runoff, tmp_calving 
  integer :: tau, taup1, n, i, j, k, ii, jj

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1

  melt         = 0.0
  pme_river    = 0.0 
  tmp_flux     = 0.0 
  liquid_precip= 0.0
  evaporation  = 0.0

  ! i and j momentum flux (Newton/m^2) crossing ocean surface
  if(.not. zero_surface_stress) then 
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            Velocity%smf(ii,jj,1) = Ice_ocean_boundary%u_flux(i,j)*Grd%umask(ii,jj,1)
            Velocity%smf(ii,jj,2) = Ice_ocean_boundary%v_flux(i,j)*Grd%umask(ii,jj,1)
         enddo
      enddo
  endif

  if (rotate_winds) then
     do j=jsc,jec
        do i=isc,iec
           tmp_x =  Grd%cos_rot(i,j)*Velocity%smf(i,j,1) + Grd%sin_rot(i,j)*Velocity%smf(i,j,2)
           tmp_y = -Grd%sin_rot(i,j)*Velocity%smf(i,j,1) + Grd%cos_rot(i,j)*Velocity%smf(i,j,2)
           Velocity%smf(i,j,1) = tmp_x
           Velocity%smf(i,j,2) = tmp_y
        enddo
     enddo
  endif

  call mpp_update_domains(Velocity%smf(:,:,1),Velocity%smf(:,:,2),Dom%domain2d,gridtype=BGRID_NE)


  ! Set temperature for water in evaporation and precipitation equal to the 
  ! ocean surface value. Default for other tracers is to have zero concentration 
  ! in evaporation and precipitation.
  do j=jsc,jec
     do i=isc,iec
       T_prog(index_temp)%tpme(i,j) = T_prog(index_temp)%field(i,j,1,tau)
     enddo
  enddo


  ! start of long if-block for use_waterflux true or false. 
  if (use_waterflux) then

      ! water flux in (kg/m^3)*(m/s) from liquid, frozen, and q_flux (evaporation)  
      if(waterflux_tavg) then 
         do j = jsc_bnd,jec_bnd
            do i = isc_bnd,iec_bnd
               ii = i + i_shift
               jj = j + j_shift
               pme(ii,jj) = 0.5*pme_taum1(ii,jj)
               pme_taum1(ii,jj) = (Ice_ocean_boundary%lprec(i,j) + Ice_ocean_boundary%fprec(i,j) &
                                  -Ice_ocean_boundary%q_flux(i,j))*Grd%tmask(ii,jj,1) 
               pme(ii,jj) = pme(ii,jj) + 0.5*pme_taum1(ii,jj)
               liquid_precip(ii,jj) =  Ice_ocean_boundary%lprec(i,j)*Grd%tmask(ii,jj,1) 
               evaporation(ii,jj)   = -Ice_ocean_boundary%q_flux(i,j)*Grd%tmask(ii,jj,1) 
            enddo
         enddo
      else 
         do j = jsc_bnd,jec_bnd
            do i = isc_bnd,iec_bnd
               ii = i + i_shift
               jj = j + j_shift
                  pme(ii,jj) = (Ice_ocean_boundary%lprec(i,j) + Ice_ocean_boundary%fprec(i,j) &
                               -Ice_ocean_boundary%q_flux(i,j))*Grd%tmask(ii,jj,1) 
                  liquid_precip(ii,jj) =  Ice_ocean_boundary%lprec(i,j)*Grd%tmask(ii,jj,1) 
                  evaporation(ii,jj)   = -Ice_ocean_boundary%q_flux(i,j)*Grd%tmask(ii,jj,1) 
            enddo
         enddo
      endif 

      ! water flux in (kg/m^3)*(m/s) from liquid runoff and calving land ice  
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            runoff(ii,jj)  = Ice_ocean_boundary%runoff(i,j)*Grd%tmask(ii,jj,1)
            calving(ii,jj) = Ice_ocean_boundary%calving(i,j)*Grd%tmask(ii,jj,1)
         enddo
      enddo
      if(runoffspread)  call spread_river_horz(runoff)
      if(calvingspread) call spread_river_horz(calving)

      ! river=runoff+calving is the water mass flux 
      ! entering ocean, other than through pme. 
      ! calving is immediately melted, so it is appropriate
      ! to add it to runoff for purposes of getting a mass
      ! flux of liquid, from land, into the ocean.   
      if(waterflux_tavg) then 
          do j=jsc,jec
             do i=isc,iec
                river(i,j)       = 0.5*river_taum1(i,j)
                river_taum1(i,j) = runoff(i,j) + calving(i,j)
                river(i,j)       = river(i,j) + 0.5*river_taum1(i,j)
             enddo
          enddo
      else 
          do j=jsc,jec
             do i=isc,iec
                river(i,j) = runoff(i,j) + calving(i,j)
             enddo
          enddo
      endif


      ! Set the temperature flux associated with the water 
      ! entering ocean from land. This flux equals to the mass 
      ! flux of water times the temperature of water, whether 
      ! this be liquid or solid water.  
      if(land_model_heat_fluxes) then  

          ! for cases where the land model computes the heat flux (W/m^2)
          ! (relative to 0degC) associated with the liquid runoff and 
          ! solid calving ice.  In this case, we have heat_land/cp_ocean
          ! as the temperature flux going into evolving the ocean temp. 
          ! Notably, we should NOT be dividing by cp_liquid_runoff nor 
          ! cp_solid_runoff to get this temperature flux.  
          do j = jsc_bnd,jec_bnd
             do i = isc_bnd,iec_bnd
                ii = i + i_shift
                jj = j + j_shift
                T_prog(index_temp)%runoff_tracer_flux(ii,jj) = cp_ocean_r &
                     *Ice_ocean_boundary%runoff_hflx(i,j)*Grd%tmask(ii,jj,1)
                T_prog(index_temp)%calving_tracer_flux(ii,jj) = cp_ocean_r &
                     *Ice_ocean_boundary%calving_hflx(i,j)*Grd%tmask(ii,jj,1)
             enddo
          enddo

          ! For diagnostic purposes, compute the effective temperatures 
          ! trunoff and tcalving.  These are NOT the temperatures that the
          ! land model may have, as those temperatures are obtained by using
          ! cp_liquid_runoff and cp_solid_runoff. Instead, this is an effective
          ! temperature resulting from dividing the land heat flux by the ocean
          ! heat capacity. We do so for purposes of diagnostics transparency 
          ! in ocean_tracer_diag...
          !
          ! also note the code to avoid non-representable numbers in 
          ! trunoff/tcalving diagnostics when runoff/calving are small.
          do j=jsc,jec
             do i=isc,iec
                T_prog(index_temp)%trunoff(i,j)  = 0.0
                T_prog(index_temp)%tcalving(i,j) = 0.0

                if(runoff(i,j) > epsln) then
                    T_prog(index_temp)%trunoff(i,j) = &
                    min(100.0, max(-273.15,T_prog(index_temp)%runoff_tracer_flux(i,j)/runoff(i,j)))
                endif
                if(calving(i,j) > epsln) then
                    T_prog(index_temp)%tcalving(i,j)= &
                    min(100.0, max(-273.15,T_prog(index_temp)%calving_tracer_flux(i,j)/calving(i,j)))
                endif
             enddo
          enddo

      else 

          ! for cases where the land model does not carry heat flux associated with 
          ! the liquid runoff and solid calving ice. We assign a temperature to the 
          ! liquid and solid runoff.  
          ! Set temperature for liquid runoff water to ocean surface value, but no 
          ! less than runoff_temp_min. For other tracers, by default keep them   
          ! set to their initial concentration.
          ! For calving temperature, set this equal to the SST.  
          do j=jsc,jec
             do i=isc,iec
                T_prog(index_temp)%trunoff(i,j)  = max(runoff_temp_min,T_prog(index_temp)%field(i,j,1,tau))
                T_prog(index_temp)%tcalving(i,j) = T_prog(index_temp)%field(i,j,1,tau)
                T_prog(index_temp)%runoff_tracer_flux(i,j)  = Grd%tmask(i,j,1)*T_prog(index_temp)%trunoff(i,j)*runoff(i,j)
                T_prog(index_temp)%calving_tracer_flux(i,j) = Grd%tmask(i,j,1)*T_prog(index_temp)%tcalving(i,j)*calving(i,j)
             enddo
          enddo

      endif

      ! compute temperature and salinity of "river" according to 
      ! temperature and salinity of calving and runoff. 
      ! allow for possibility of runoff and/or calving to be negative,
      ! which may exist for certain land models that "suck" water from 
      ! the oceans to address limitations in their formulation. 
      ! We do not include such negative points when computing triver. 
      ! Also compute the salinity flux associated with liquid and solid runoff. 
      do j=jsc,jec
         do i=isc,iec
            tmp_runoff = max(0.0,runoff(i,j))
            tmp_calving= max(0.0,calving(i,j))
            T_prog(index_temp)%triver(i,j) = Grd%tmask(i,j,1)    &
                 *(tmp_runoff *T_prog(index_temp)%trunoff(i,j)   &
                  +tmp_calving*T_prog(index_temp)%tcalving(i,j)) &
                 /(epsln + tmp_runoff + tmp_calving)

            T_prog(index_salt)%triver(i,j) = Grd%tmask(i,j,1)    &
                 *(tmp_runoff *T_prog(index_salt)%trunoff(i,j)   &
                  +tmp_calving*T_prog(index_salt)%tcalving(i,j)) &
                 /(epsln + tmp_runoff + tmp_calving)
         enddo
      enddo


      ! for debugging
      if(debug_water_fluxes) then 

          write(stdoutunit,'(a)')

          wrk1_2d(:,:) = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*cp_ocean &
                               *T_prog(index_temp)%runoff_tracer_flux(i,j)*dtime
             enddo
          enddo
          tracer_input = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), global_sum_flag)

          ! truncate out the lower order bits if using non_bitwise_reproducible sums
          if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) tracer_input = real(tracer_input, kind=FLOAT_KIND)
          write (stdoutunit,'(a,es24.17,a)') ' Heat input via liquid runoff in ocean_sbc          = ',&
                                             tracer_input,' Joule'
          wrk1_2d(:,:) = 0.0
          do j=jsc,jec
             do i=isc,iec
                wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*cp_ocean &
                               *T_prog(index_temp)%calving_tracer_flux(i,j)*dtime
             enddo
          enddo
          tracer_input = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), global_sum_flag)

          ! truncate out the lower order bits if using non_bitwise_reproducible sums
          if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) tracer_input = real(tracer_input, kind=FLOAT_KIND)
          write (stdoutunit,'(a,es24.17,a)') ' Heat input via calving land ice in ocean_sbc       = ',&
                                           tracer_input,' Joule'

          write(stdoutunit,'(a)')

          if(zero_water_fluxes) then 
              river     = 0.0
              pme       = 0.0
              pme_taum1 = 0.0
              runoff    = 0.0
              calving   = 0.0
          endif
          if(zero_river_fluxes) then 
              river     = 0.0
          endif
          if(zero_pme_fluxes) then 
              pme       = 0.0
              pme_taum1 = 0.0
          endif
          if(zero_runoff_fluxes) then 
              runoff    = 0.0
          endif
          if(zero_calving_fluxes) then 
              calving   = 0.0
          endif
          if(convert_river_to_pme) then 
              do j=jsc,jec
                 do i=isc,iec 
                    pme(i,j) = pme(i,j) + river(i,j)
                 enddo
              enddo
              river   = 0.0
              calving = 0.0
              runoff  = 0.0
          endif

      endif

      ! set riverdiffuse to determine where to enhance diff_cbt inside ocean_rivermix_mod
      do n=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%riverdiffuse(i,j) = river(i,j) 
            enddo
         enddo
      enddo

      ! Allow for nonzero salt flux from, say, ice melt.  
      ! flux from ice-melt is in units of (kg/m^3)*(m/sec).
      ! to convert to rho*dz*salinity units for the model, 
      ! need to multiply by 1000.0.
      ! 
      ! The melt field is deduced from knowing that 
      ! salt transferred between the ice and ocean
      ! occurs in the presence of liquid water transport
      ! between the ice and ocean. When salt is added to 
      ! the ocean from the ice melt, fresh liquid water is 
      ! also added.
      ! 
      ! The minus sign arises since the ice model produces a 
      ! salt_flux > 0 when there is salt added to the ice model,
      ! and hence taken away from the ocean model.  
      melt = 0.0
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift      
            T_prog(index_salt)%stf(ii,jj) = -Ice_ocean_boundary%salt_flux(i,j)*1000.0
            melt(ii,jj)                   = -Ice_ocean_boundary%salt_flux(i,j)*ice_salt_concentration_r
         enddo
      enddo

      ! produce a zero area average of pme + river. 
      ! note that we remove the ice melt water  
      ! since the coupler has already included the 
      ! ice melt within the pme field. 
      pme_river = 0.0 
      if(zero_net_water_coupler) then 
         do j=jsc,jec
            do i=isc,iec
               pme_river(i,j) = pme(i,j) + river(i,j) - melt(i,j)
            enddo
         enddo
         pme_river_total = mpp_global_sum(Dom%domain2d,pme_river(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1),&
                                          global_sum_flag)/Grd%tcellsurf

         ! truncate out the lower order bits if using non_bitwise_reproducible sums
         if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_river_total = real(pme_river_total, kind=FLOAT_KIND)
         do j=jsc,jec
            do i=isc,iec
               pme(i,j)  = pme(i,j) - pme_river_total*Grd%tmask(i,j,1)
            enddo
         enddo
      endif 
      
      ! when have a nonzero salinity of runoff, then there is a 
      ! nonzero salt flux (kg/(m^2*sec)) into the ocean with the runoff. 
      do j=jsc,jec
         do i=isc,iec
            T_prog(index_salt)%runoff_tracer_flux(i,j) = &
            T_prog(index_salt)%trunoff(i,j)*runoff(i,j)*T_prog(index_salt)%conversion 
         enddo
      enddo


  else  ! now code for use_waterflux=.false.  

     
      ! water flux in (kg/m^3)*(m/s) from rivers and calving land glaciers  
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            runoff(ii,jj)  = Ice_ocean_boundary%runoff(i,j)*Grd%tmask(ii,jj,1)
            calving(ii,jj) = Ice_ocean_boundary%calving(i,j)*Grd%tmask(ii,jj,1)
         enddo
      enddo
      if(runoffspread)  call spread_river_horz(runoff)
      if(calvingspread) call spread_river_horz(calving)     
      do j=jsc,jec
         do i=isc,iec
            river(i,j) = runoff(i,j) + calving(i,j)
         enddo
      enddo

      ! convert salt flux from ice to psu with multiplication by 1000.0. 
      ! convert freshwater mass fluxes (kg/m^3)*(m/sec) into virtual salt
      ! fluxes with multiplication by salinity_ref. 
         do j = jsc_bnd, jec_bnd
            do i = isc_bnd, iec_bnd
               ii = i + i_shift
               jj = j + j_shift
               T_prog(index_salt)%stf(ii,jj) = -Ice_ocean_boundary%salt_flux(i,j)*1000.0     -&
                                               (Ice_ocean_boundary%lprec(i,j)                +&
                                                Ice_ocean_boundary%fprec(i,j) + river(ii,jj) -&
                                                Ice_ocean_boundary%q_flux(i,j))*salinity_ref*Grd%tmask(ii,jj,1) 
         enddo
      enddo
      ! set the riverdiffuse "mask" to determine where to enhance diff_cbt
      do n=1,num_prog_tracers
         do j = jsc,jec
            do i = isc, iec
               T_prog(n)%riverdiffuse(i,j) = river(i,j) 
            enddo
         enddo
      enddo

      ! set pme and river to zero since use_waterflux=.false. 
      pme     = 0.0
      river   = 0.0
      runoff  = 0.0
      calving = 0.0


  endif    ! end of long if-block for if(use_waterflux) then 


  call mpp_update_domains(pme(:,:)  , Dom%domain2d)
  call mpp_update_domains(river(:,:), Dom%domain2d)

  
  if(.not. zero_heat_fluxes) then 
    do j = jsc_bnd, jec_bnd
       do i = isc_bnd, iec_bnd
          ii = i + i_shift
          jj = j + j_shift  
          T_prog(index_temp)%stf(ii,jj) = (Ice_ocean_boundary%sw_flux_vis_dir(i,j) + &
                                           Ice_ocean_boundary%sw_flux_vis_dif(i,j) + &
                                           Ice_ocean_boundary%sw_flux_nir_dir(i,j) + &
                                           Ice_ocean_boundary%sw_flux_nir_dif(i,j) + &
                                           Ice_ocean_boundary%lw_flux(i,j)         - &
                                           (Ice_ocean_boundary%fprec(i,j)          + &
                                            calving(ii,jj))*hlf                    - &
                                           Ice_ocean_boundary%t_flux(i,j)          - &
                                           Ice_ocean_boundary%q_flux(i,j)*hlv        &
                                          )/cp_ocean*Grd%tmask(ii,jj,1) 
       enddo
    enddo
  endif 


  ! set velocity of pme and river water to that of upper ocean cell.   
  ! generalizations may be suitable with refined component models.   
  do n = 1, size(upme,3)
     do j = jsc, jec
        do i = isc, iec
           upme(i,j,n)   = Velocity%u(i,j,1,n,tau)   ! velocity of precip-evap water
           uriver(i,j,n) = Velocity%u(i,j,1,n,tau)   ! velocity of river water 
        enddo
     enddo
  enddo

  ! apply atmosphere and ice pressure to ocean only when there is
  ! mass transfer allowed between ocean and atmos/ice via fresh water. 
  ! if do not use fresh water, then cannot allow ocean to feel 
  ! weight of atmosphere and ice.  
  !
  ! NOTE: the option use_full_patm_for_sea_level allows for the passing 
  ! of the sea level including the full weight of sea ice back to
  ! the ice model.  This approach maintains the max weight on the liquid
  ! ocean according to the nml variable max_ice_thickness.  But it does 
  ! allow the sea ice to know when there is actually more sea ice than that
  ! set by max_ice_thickness.  This option then provides for a negative
  ! feedback on the runaway growth of sea ice, since the full pressure acting to 
  ! make the ice flow will be correctly felt.  This is a new option, and is not
  ! fully tested, so the default is use_full_patm_for_sea_level=.false.
  ! 
  tmp_patm(:,:) = 0.0
  if(use_waterflux) then 
      var = grav*rho0*max_ice_thickness
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift
            patm(ii,jj)     = min(Ice_ocean_boundary%p(i,j),var)
            tmp_patm(ii,jj) = Ice_ocean_boundary%p(i,j)
         enddo
      enddo
      call mpp_update_domains(patm(:,:), Dom%domain2d)
      call mpp_update_domains(tmp_patm(:,:), Dom%domain2d)
      if(use_full_patm_for_sea_level) then 
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%patm_t(i,j,taup1)     = patm(i,j) 
                Ext_mode%patm_for_sea_lev(i,j) = tmp_patm(i,j)
             enddo
          enddo
      else
          do j=jsd,jed
             do i=isd,ied
                Ext_mode%patm_t(i,j,taup1)     = patm(i,j) 
                Ext_mode%patm_for_sea_lev(i,j) = patm(i,j) 
             enddo
          enddo
      endif
  endif

  ! shortwave flux (W/m^2) into ocean
  if(.not. zero_heat_fluxes) then 
      swflx(isc:iec,jsc:jec) = Grd%tmask(isc:iec,jsc:jec,1)                  &
       *(Ice_ocean_boundary%sw_flux_nir_dir(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       + Ice_ocean_boundary%sw_flux_nir_dif(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       + Ice_ocean_boundary%sw_flux_vis_dir(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       + Ice_ocean_boundary%sw_flux_vis_dif(isc_bnd:iec_bnd,jsc_bnd:jec_bnd))

      swflx_vis(isc:iec,jsc:jec) = Grd%tmask(isc:iec,jsc:jec,1)               &
        *(Ice_ocean_boundary%sw_flux_vis_dir(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) &
       +  Ice_ocean_boundary%sw_flux_vis_dif(isc_bnd:iec_bnd,jsc_bnd:jec_bnd))
  endif 

  call ocean_tpm_sbc(Dom, Grd, T_prog(:), Time, Ice_ocean_boundary%fluxes, runoff, &
                     isc_bnd, iec_bnd, jsc_bnd, jec_bnd)


  ! send fields to diagnostic manager 

  ! i-directed wind stress (N/m^2)
  if (id_tau_x > 0) used =  send_data(id_tau_x, Velocity%smf(:,:,1), &
                            Time%model_time, rmask=Grd%umask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! j-directed wind stress (N/m^2)
  if (id_tau_y > 0) used =  send_data(id_tau_y, Velocity%smf(:,:,2), &
                            Time%model_time, rmask=Grd%umask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

  ! wind stress curl (N/m^3) averaged to U-point
  ! Ekman pumping velocity averaged to U-point 
  if (id_tau_curl > 0 .or. id_ekman_we > 0) then 
      wrk2_2d(:,:) = 0.0
      wrk3_2d(:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            umask_norm   = Grd%umask(i,j,1)/(epsln+Grd%umask(i+1,j,1)+Grd%umask(i-1,j,1))
            wrk2_2d(i,j) = umask_norm*(                                                             &
                  Grd%umask(i+1,j,1)*( Velocity%smf(i+1,j,2)-Velocity%smf(i,j,2))  /Grd%dxtn(i+1,j) & 
                 +Grd%umask(i-1,j,1)*( Velocity%smf(i,j,2)  -Velocity%smf(i-1,j,2))/Grd%dxtn(i,j)   &
                 )
            umask_norm   = Grd%umask(i,j,1)/(epsln+Grd%umask(i,j+1,1)+Grd%umask(i,j-1,1))
            wrk2_2d(i,j) = wrk2_2d(i,j)                                                             &
                 -umask_norm*(                                                                      &
                  Grd%umask(i,j+1,1)*( Velocity%smf(i,j+1,1)-Velocity%smf(i,j,1))  /Grd%dyte(i,j+1) & 
                 +Grd%umask(i,j-1,1)*( Velocity%smf(i,j,1)  -Velocity%smf(i,j-1,1))/Grd%dyte(i,j)   &
                 )
            if(abs(Grd%yu(i,j)) > 1.0) then 
              wrk3_2d(i,j) = rho0r*wrk2_2d(i,j)/Grd%f(i,j) 
            endif 
         enddo
      enddo
      used = send_data(id_tau_curl, wrk2_2d(:,:),               & 
                       Time%model_time, rmask=Grd%umask(:,:,1), &
                       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
      used = send_data(id_ekman_we, wrk3_2d(:,:),               & 
                       Time%model_time, rmask=Grd%umask(:,:,1), &
                       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

  endif 


  ! Ekman Heat Transport component
  ! Defined as: Hek = -rho_cp \Int[ Taux/(f rho0) (T1 - <T>) ]dx
  ! where T1=first layer Temp; < >=vertical mean.
  ! Zonal sum done off-line.
  ! TAUx = smf*rho0, we want TAUx/rho0, therefore we use smf.  
  ! Algorithm from riccardo.farneti@noaa.gov
  if (id_ekman_heat > 0) then
  
    wrk1_2d(:,:) = 0.0
    do j=jsc,jec
      do i=isc,iec

        factor = 0.0
        if(Grd%f(i,j)==0.0 .and. j>1) then
           factor = 4.*Grd%f(i,j-1)
        else
           factor = 4.*Grd%f(i,j)
        endif

        totz = 1.0
        tbrz = 0.0
        maskt= 0.0
        do k=1,nk 
           maskt = Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
           tbrz  = tbrz + min(Thickness%dzt(i,j,k),Thickness%dzt(i,j+1,k))*Thickness%dzt(i,j,k)*maskt &
                   *(T_prog(index_temp)%field(i,j,k,tau)+T_prog(index_temp)%field(i,j+1,k,tau))
           totz  = totz + min(Thickness%dzt(i,j,k),Thickness%dzt(i,j+1,k))*Thickness%dzt(i,j,k)*maskt 
        enddo    
        if (totz /= 0.0) then
           tbrz = tbrz/totz
           wrk1_2d(i,j) = wrk1_2d(i,j) - ( Velocity%smf(i,j,1)*Grd%dxu(i,j) + Velocity%smf(i-1,j,1)*Grd%dxu(i-1,j) )&
                          *( T_prog(index_temp)%field(i,j,1,tau) + T_prog(index_temp)%field(i,j+1,1,tau) -tbrz )    &
                          *cos(Grd%phiu(i,j)) / factor    
        endif

      enddo
    enddo
  
    if(id_ekman_heat > 0) then 
        used = send_data(id_ekman_heat, wrk1_2d(:,:)*rho_cp,&
             Time%model_time, rmask=Grd%tmask(:,:,1),       &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

  endif



  ! temp of runoff
  if (id_trunoff(index_temp) > 0) then

      wrk1_2d=0.0
      if(land_model_heat_fluxes) then 
          ! temp_runoff = heat_runoff/(cp_liquid_runoff*runoff).
          ! this is the temperature that the land model would 
          ! call its runoff temperature. 
          do j=jsc,jec
             do i=isc,iec
                if(runoff(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%trunoff(i,j)*cp_liquid_runoff_r*cp_ocean
                endif
             enddo
          enddo

      else 

          do j=jsc,jec
             do i=isc,iec
                if(runoff(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%trunoff(i,j)
                endif
             enddo
          enddo

      endif

      used = send_data(id_trunoff(index_temp),     &
           wrk1_2d(:,:),                           &
           Time%model_time, rmask=Grd%tmask(:,:,1),&
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


  ! salinity of runoff
  if (id_trunoff(index_salt) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(runoff(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      used = send_data(id_trunoff(index_salt),             &
             wrk1_2d(:,:)*T_prog(index_salt)%trunoff(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1),      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! temp of calving solid runoff 
  if (id_tcalving(index_temp) > 0) then

      wrk1_2d=0.0
      if(land_model_heat_fluxes) then 
          ! temp_calving = heat_calving/(cp_solid_runoff*calving).
          ! this is the temperature that the land model would 
          ! call its calving land ice temperature. 
          do j=jsc,jec
             do i=isc,iec
                if(calving(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%tcalving(i,j)*cp_solid_runoff_r*cp_ocean
                endif
             enddo
          enddo

      else 

          do j=jsc,jec
             do i=isc,iec
                if(calving(i,j) /= 0.0) then 
                    wrk1_2d(i,j) = T_prog(index_temp)%tcalving(i,j)
                endif
             enddo
          enddo

      endif

      used = send_data(id_tcalving(index_temp),    &
           wrk1_2d(:,:),                           &
           Time%model_time, rmask=Grd%tmask(:,:,1),&
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! effective temp of calving, computed as 
  ! temp_calve = heat_calve/(cp_ocean*calve)
  if (id_temp_calving_eff > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(calving(i,j) /= 0.0) then 
                wrk1_2d(i,j) = 1.0
            endif 
         enddo
      enddo
      used = send_data(id_temp_calving_eff,                &
             wrk1_2d(:,:)*T_prog(index_temp)%tcalving(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1),       &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! salinity of calving
  if (id_tcalving(index_salt) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(calving(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      used = send_data(id_tcalving(index_salt),            &
             wrk1_2d(:,:)*T_prog(index_salt)%tcalving(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! temp of river 
  if (id_triver(index_temp) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(river(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      used = send_data(id_triver(index_temp),            &
             wrk1_2d(:,:)*T_prog(index_temp)%triver(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! salinity of river 
  if (id_triver(index_salt) > 0) then
      wrk1_2d=0.0
      do j=jsc,jec
         do i=isc,iec
            if(river(i,j) /= 0.0) wrk1_2d(i,j) = 1.0
         enddo
      enddo
      used = send_data(id_triver(index_salt),            &
             wrk1_2d(:,:)*T_prog(index_salt)%triver(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! surface heat flux (W/m2) passed through the coupler 
  if (id_stf_coupler(index_temp) > 0) then
      used = send_data(id_stf_coupler(index_temp),                     &
             T_prog(index_temp)%stf(:,:)*T_prog(index_temp)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total surface heat flux (Watts) passed through coupler 
  if(id_total_ocean_stf_coupler(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                    *T_prog(index_temp)%stf(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_coupler(index_temp), total_stuff*1e-15, Time%model_time)
  endif

  ! heat input from liquid river runoff relative to 0 degrees C (W/m2)
  if (id_stf_runoff(index_temp) > 0) then
      used = send_data(id_stf_runoff(index_temp),                                     &
             T_prog(index_temp)%runoff_tracer_flux(:,:)*T_prog(index_temp)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                 &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total heat flux from liquid river runoff (Watts), relative to 0C. 
  if(id_total_ocean_stf_runoff(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                     *T_prog(index_temp)%runoff_tracer_flux(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_runoff(index_temp), total_stuff*1e-15, Time%model_time)
  endif

  ! heat input from solid calving land ice relative to 0 degrees C (W/m2)
  if (id_stf_calving(index_temp) > 0) then
      used = send_data(id_stf_calving(index_temp),                                     &
             T_prog(index_temp)%calving_tracer_flux(:,:)*T_prog(index_temp)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total heat flux from solid calving land ice (Watts), relative to 0C. 
  if(id_total_ocean_stf_calving(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
         *T_prog(index_temp)%calving_tracer_flux(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_calving(index_temp), total_stuff*1e-15, Time%model_time)
  endif

  ! heat input from liquid precip relative to 0 degrees C (W/m2).
  ! note that frozen precip arrives at 0C, so contributes no heat 
  ! relative to 0C.  Assume temp of liquid precip same as tpme
  if (id_stf_prec(index_temp) > 0) then
      used = send_data(id_stf_prec(index_temp),                                             &
             liquid_precip(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion, &
             Time%model_time, rmask=Grd%tmask(:,:,1),                                       &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total heat flux from liquid precip (Watts)
  if(id_total_ocean_stf_prec(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                     *liquid_precip(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_prec(index_temp), total_stuff*1e-15, Time%model_time)
  endif

  ! heat sent away from ocean due to water mass leaving ocean
  ! via evaporation, measured relative to 0 degrees C (W/m2).
  ! Assume temp of evaporating water is same as tpme
  if (id_stf_evap(index_temp) > 0) then
      used = send_data(id_stf_evap(index_temp),                                          &
             evaporation(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total heat flux from evaporating water carrying heat away from ocean (Watts)
  if(id_total_ocean_stf_evap(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                     *evaporation(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_evap(index_temp), total_stuff*1e-15, Time%model_time)
  endif


  ! salt flux (kg/(m2*sec)) passed through the coupler 
  if (id_stf_coupler(index_salt) > 0) then
      used = send_data(id_stf_coupler(index_salt),                      &
             T_prog(index_salt)%stf(:,:)*T_prog(index_salt)%conversion, &
             Time%model_time, rmask=Grd%tmask(:,:,1),                   &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total salt flux (kg/sec) passed through coupler 
  if(id_total_ocean_stf_coupler(index_salt) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                    *T_prog(index_salt)%stf(:,:)*T_prog(index_salt)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_coupler(index_salt), total_stuff*1e-15, Time%model_time)
  endif

  ! salt input from liquid river runoff (kg/(m2*sec))
  if (id_stf_runoff(index_salt) > 0) then
      used = send_data(id_stf_runoff(index_salt),                                     &
             T_prog(index_salt)%runoff_tracer_flux(:,:)*T_prog(index_salt)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                 &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! salt input from calving land ice (kg/(m2*sec))
  if (id_stf_runoff(index_salt) > 0) then
      used = send_data(id_stf_runoff(index_salt),                                      &
             T_prog(index_salt)%calving_tracer_flux(:,:)*T_prog(index_salt)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total salt flux from liquid river runoff (kg/sec)
  if(id_total_ocean_stf_runoff(index_salt) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                    *T_prog(index_salt)%runoff_tracer_flux(:,:)*T_prog(index_salt)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_runoff(index_salt), total_stuff*1e-15, Time%model_time)
  endif

  ! total salt flux from calving land ice (kg/sec)
  if(id_total_ocean_stf_calving(index_salt) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                    *T_prog(index_salt)%calving_tracer_flux(:,:)*T_prog(index_salt)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_calving(index_salt), total_stuff*1e-15, Time%model_time)
  endif

  ! salt input from ice (kg/(m2*sec))
  if (id_salt_flux_ice > 0) then
      used = send_data(id_salt_flux_ice,             &
             melt(:,:)*ice_salt_concentration,       &
             Time%model_time, rmask=Grd%tmask(:,:,1),&
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total salt flux from ice (kg/sec)
  if(id_total_salt_flux_ice > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                    *melt(:,:)*ice_salt_concentration
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_salt_flux_ice, total_stuff*1e-15, Time%model_time)
  endif


  ! total mass flux per area from pme and river (kg/(m2*sec))  
  if (id_pme_river > 0) then
      used = send_data(id_pme_river, pme(:,:)+river(:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,1),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total mass flux from pme+river (kg/sec)
  if(id_total_ocean_pme_river > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*(pme(:,:)+river(:,:))
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_pme_river, total_stuff*1e-15, Time%model_time)
  endif


  ! mass flux per area from pme_sbc (kg/(m2*sec))  
  if (id_pme_sbc > 0) then
      used = send_data(id_pme_sbc, pme(:,:),             &
             Time%model_time, rmask=Grd%tmask(:,:,1),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total mass flux from pme_sbc (kg/sec)
  if(id_total_ocean_pme_sbc > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*pme(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_pme_sbc, total_stuff*1e-15, Time%model_time)
  endif

  ! mass flux per area from ice melt (kg/(m2*sec))  
  if (id_melt > 0) then
      used = send_data(id_melt, melt(:,:),               & 
             Time%model_time, rmask=Grd%tmask(:,:,1),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total mass flux from ice melt (kg/sec)
  if(id_total_ocean_melt > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*melt(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_melt, total_stuff*1e-15, Time%model_time)
  endif

  ! evaporative mass flux (kg/(m2*sec))
  if (id_evap > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%q_flux(i,j)
         enddo
      enddo
      used = send_data(id_evap, tmp_flux(:,:),    &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total mass flux from evap (kg/sec)
  if(id_total_ocean_evap > 0) then 
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%q_flux(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_evap, total_stuff*1e-15, Time%model_time)
  endif

  ! frozen precip (kg/(m2*sec))
  if (id_fprec > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%fprec(i,j)
         enddo
      enddo
      used = send_data(id_fprec, tmp_flux(:,:),        &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total frozen precip (kg/sec)
  if (id_total_ocean_fprec > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%fprec(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_fprec, total_stuff*1e-15, Time%model_time)
  endif


  ! liquid precip (kg/(m2*sec))
  if (id_lprec > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%lprec(i,j)
         enddo
      enddo
      used = send_data(id_lprec, tmp_flux(:,:),        &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total liquid precip (kg/sec)
  if (id_total_ocean_lprec > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%lprec(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_lprec, total_stuff*1e-15, Time%model_time)
  endif

  ! river (mass flux of land water (liquid+solid) ) entering ocean (kg/m^3)*(m/s)
  if (id_river > 0) used =  send_data(id_river, river(:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,1),             &
       is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
  ! global sum of river input (kg/sec)
  if(id_total_ocean_river > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*river(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_river, total_stuff*1e-15, Time%model_time)
  endif

  ! calving land ice (kg/(m2*sec)) entering the ocean 
  if (id_calving > 0) then 
      used = send_data(id_calving, calving(:,:),   &
           Time%model_time, rmask=Grd%tmask(:,:,1),&
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
  endif
  ! total mass of calving (kg/sec)
  if (id_total_ocean_calving > 0) then
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*calving(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_calving, total_stuff*1e-15, Time%model_time)
  endif

  ! liquid river runoff entering the ocean (kg/m^3)*(m/s)
  if (id_runoff > 0) then 
      used = send_data(id_runoff, runoff(:,:),      &
           Time%model_time, rmask=Grd%tmask(:,:,1), &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
  endif
  ! total liquid river runoff (kg/sec)
  if (id_total_ocean_runoff > 0) then
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*runoff(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_runoff, total_stuff*1e-15, Time%model_time)
  endif

  ! shortwave flux (W/m2)
  if (id_swflx > 0) used =  send_data(id_swflx, swflx(:,:),          &
                            Time%model_time, rmask=Grd%tmask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  ! total shortwave flux (Watts)
  if (id_total_ocean_swflx > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*swflx(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_swflx, total_stuff*1e-15, Time%model_time)
  endif 


  ! visible shortwave flux (W/m2)
  if (id_swflx_vis > 0) used =  send_data(id_swflx_vis, swflx_vis(:,:),  &
                                Time%model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  ! total visible shortwave flux (Watts)
  if (id_total_ocean_swflx_vis > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*swflx_vis(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_swflx_vis, total_stuff*1e-15, Time%model_time)
  endif 


  ! evaporative heating (W/m2) (<0 cools ocean)
  if (id_evap_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -hlv*Ice_ocean_boundary%q_flux(i,j)
         enddo
      enddo
      used = send_data(id_evap_heat, tmp_flux(:,:),    &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total evaporative heating (Watts) 
  if (id_total_ocean_evap_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -hlv*Ice_ocean_boundary%q_flux(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_evap_heat, total_stuff*1e-15, Time%model_time)
  endif

  ! longwave heating (W/m2)
  if (id_lw_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%lw_flux(i,j)
         enddo
      enddo
      used = send_data(id_lw_heat, tmp_flux(:,:),    &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total longwave heating (Watts) 
  if (id_total_ocean_lw_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = Ice_ocean_boundary%lw_flux(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_lw_heat, total_stuff*1e-15, Time%model_time)
  endif

  ! heating from melting the frozen precip (W/m2)
  if (id_fprec_melt_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -Ice_ocean_boundary%fprec(i,j)*hlf
         enddo
      enddo
      used = send_data(id_fprec_melt_heat, tmp_flux(:,:),    &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total heating from melting the frozen precip (Watts) 
  if (id_total_ocean_fprec_melt_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -Ice_ocean_boundary%fprec(i,j)*hlf
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_fprec_melt_heat, total_stuff*1e-15, Time%model_time)
  endif

  ! heating from the melting of calved land ice (W/m2)
  if (id_calving_melt_heat > 0) then
      used = send_data(id_calving_melt_heat, -calving(:,:)*hlf, &
             Time%model_time, rmask=Grd%tmask(:,:,1),      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total heat from the melting of calved land ice (Watts)  
  if (id_total_ocean_calving_melt_heat > 0) then
      wrk1_2d(:,:) = -hlf*Grd%tmask(:,:,1)*Grd%dat(:,:)*calving(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_calving_melt_heat, total_stuff*1e-15, Time%model_time)
  endif

  ! sensible heat flux (W/m2)
  if (id_sens_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -Ice_ocean_boundary%t_flux(i,j)
         enddo
      enddo
      used = send_data(id_sens_heat, tmp_flux(:,:),    &
             Time%model_time, rmask=Grd%tmask(:,:,1),  &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total sensible heat flux (Watts) 
  if (id_total_ocean_sens_heat > 0) then
      do j=jsc_bnd,jec_bnd
         do i=isc_bnd,iec_bnd
            ii=i+i_shift
            jj=j+j_shift  
            tmp_flux(ii,jj) = -Ice_ocean_boundary%t_flux(i,j)
         enddo
      enddo
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*tmp_flux(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_sens_heat, total_stuff*1e-15, Time%model_time)
  endif


end subroutine get_ocean_sbc
! </SUBROUTINE> NAME="get_ocean_sbc"


!#######################################################################
! <SUBROUTINE NAME="flux_adjust">
!
! <DESCRIPTION>
! Subroutine to compute the surface fluxes derived from a 
! restoring condition and/or correction from an input file. 
!
! We use a convention whereby a positive 
! flux enters the ocean:  (+) down convention. 
!
! When restoring salinity, one may choose to convert this
! flux to an implied water flux, or keep it a salt flux.
! The default is to keep it as a salt flux.  Converting to 
! a water flux will alter the sea level, and so alter the 
! concentration of other tracers.  
!
! </DESCRIPTION>
!
subroutine flux_adjust(Time, T_diag, T_prog, Velocity, river, melt, pme)

  type(ocean_time_type),         intent(in)    :: Time
  type(ocean_diag_tracer_type),  intent(in)    :: T_diag(:)
  type(ocean_prog_tracer_type),  intent(inout) :: T_prog(:)
  type(ocean_velocity_type),     intent(inout) :: Velocity
  real, dimension(isd:,jsd:),    intent(in)    :: river
  real, dimension(isd:,jsd:),    intent(in)    :: melt 
  real, dimension(isd:,jsd:),    intent(inout) :: pme

  real, dimension(isd:ied,jsd:jed) :: open_ocean_mask
  real, dimension(isd:ied,jsd:jed) :: pme_restore, flx_restore
  real, dimension(isd:ied,jsd:jed) :: pme_correct, flx_correct
  real                             :: pme_river_total
  real                             :: pme_restore_total, flx_restore_total 
  real                             :: pme_correct_total, flx_correct_total 
  real                             :: total_stuff 
  real                             :: tmp_delta_salinity 
  integer                          :: i, j, tau, taum1
  logical                          :: used
  
  tau      = Time%tau
  taum1    = Time%taum1  

  pme_correct = 0.0
  pme_restore = 0.0
  flx_correct = 0.0 
  flx_restore = 0.0 

  open_ocean_mask(isd:ied,jsd:jed) = Grd%tmask(isd:ied,jsd:jed,1)

  ! add restoring to fluxes from coupled model, or
  ! add flux correction to fluxes from coupled model.
  ! NOTE: (+) down convention

  if (id_restore(index_salt) > 0 .and. salt_damp_factor > 0.0) then
      call time_interp_external(id_restore(index_salt), Time%model_time, data)

      ! initialization of restoring fields 
      pme_restore = 0.0
      flx_restore = 0.0 

      ! use near-frazil condition as a proxy for where sea-ice is present 
      if(.not. salt_restore_under_ice) then 
          do j=jsc,jec
             do i=isc,iec
                if(Grd%tmask(i,j,1) == 1.0) then 
                    if(T_prog(index_temp)%field(i,j,1,tau) <= -0.0539*T_prog(index_salt)%field(i,j,1,tau)) then 
                        open_ocean_mask(i,j) = 0.0
                    endif
                endif
             enddo
          enddo
      endif

      if (use_waterflux .and. .not. salt_restore_as_salt_flux) then

          ! put all water restore into pme_restore (no river_restore has been coded) 
          if(max_delta_salinity_restore < 0.0) then 
              do j=jsc,jec
                 do i=isc,iec
                    pme_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                         *(T_prog(index_salt)%field(i,j,1,taum1)-data(i,j))                    &
                         /(T_prog(index_salt)%field(i,j,1,taum1)+epsln)
                 enddo
              enddo
          else 
              do j=jsc,jec
                 do i=isc,iec
                    tmp_delta_salinity = T_prog(index_salt)%field(i,j,1,taum1)-data(i,j)
                    tmp_delta_salinity = sign(1.0,tmp_delta_salinity) &
                                        *min(abs(tmp_delta_salinity),max_delta_salinity_restore)
                    pme_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                         *tmp_delta_salinity                                                   &
                         /(T_prog(index_salt)%field(i,j,1,taum1)+epsln)
                 enddo
              enddo
          endif

          ! produce a zero area average so there is no net input 
          ! of pme mass to the ocean associated with the restoring 
          if(zero_net_water_restore) then 
             pme_restore_total = mpp_global_sum(Dom%domain2d,pme_restore(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                             global_sum_flag)/Grd%tcellsurf

             ! truncate out the lower order bits if using non_bitwise_reproducible sums
             if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_restore_total = real(pme_restore_total, kind=FLOAT_KIND)
             pme_restore(isc:iec,jsc:jec) = pme_restore(isc:iec,jsc:jec) - pme_restore_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add pme_restore to pme
          pme(isc:iec,jsc:jec) = pme(isc:iec,jsc:jec) + pme_restore(isc:iec,jsc:jec)

          ! produce a zero area average of pme_restore + pme + river.
          ! remove the melt term from sea ice since the coupler 
          ! has already included the ice melt within the pme field. 
          pme_river(:,:) = 0.0
          if(zero_net_water_couple_restore) then 
              do j=jsc,jec
                 do i=isc,iec
                    pme_river(i,j) = pme(i,j) + river(i,j) - melt(i,j)
                 enddo
              enddo
              pme_river_total = mpp_global_sum(Dom%domain2d,pme_river(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                                global_sum_flag)/Grd%tcellsurf

              ! truncate out the lower order bits if using non_bitwise_reproducible sums
              if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_river_total = real(pme_river_total, kind=FLOAT_KIND)
              do j=jsc,jec
                 do i=isc,iec
                    pme(i,j)  = pme(i,j) - pme_river_total*Grd%tmask(i,j,1)
                 enddo
              enddo
          endif
          call mpp_update_domains(pme(:,:), Dom%domain2d)


      else  ! for opposite of if (use_waterflux .and. .not. salt_restore_as_salt_flux) then

          if(max_delta_salinity_restore < 0.0) then 
              do j=jsc,jec
                 do i=isc,iec
                    flx_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                                       *(data(i,j) - T_prog(index_salt)%field(i,j,1,taum1))
                 enddo
              enddo
          else 
              do j=jsc,jec
                 do i=isc,iec
                    tmp_delta_salinity = data(i,j) - T_prog(index_salt)%field(i,j,1,taum1)
                    tmp_delta_salinity = sign(1.0,tmp_delta_salinity) &
                                        *min(abs(tmp_delta_salinity),max_delta_salinity_restore)
                    flx_restore(i,j) = salt_damp_factor*open_ocean_mask(i,j)*restore_mask(i,j) &
                                       *tmp_delta_salinity
                 enddo
              enddo
          endif

          ! produce a zero area average so there is no net input
          ! of salt to the ocean associated with the restoring 
          if(zero_net_salt_restore) then 
            flx_restore_total =  &
            mpp_global_sum(Dom%domain2d,flx_restore(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), global_sum_flag)/Grd%tcellsurf

            ! truncate out the lower order bits if using non_bitwise_reproducible sums
            if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) flx_restore_total = real(flx_restore_total, kind=FLOAT_KIND)
            flx_restore(isc:iec,jsc:jec) = flx_restore(isc:iec,jsc:jec) - flx_restore_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add flx_restore to stf
          do j = jsc,jec
             do i = isc,iec
                T_prog(index_salt)%stf(i,j) = T_prog(index_salt)%stf(i,j) + flx_restore(i,j)
             enddo
          enddo

      endif  ! endif for if (use_waterflux .and. .not. salt_restore_as_salt_flux) then

  endif  ! endif for if(id_restore(index_salt) > 0 .and. salt_damp_factor > 0.0)


  ! add salt fluxes from a data file to perform a flux correction
  if (id_correction(index_salt) > 0 ) then

      call time_interp_external(id_correction(index_salt), Time%model_time, data)

      ! initialization of correction fields 
      pme_correct = 0.0
      flx_correct = 0.0 

      ! use near-frazil condition as a proxy for where sea-ice is present 
      if(.not. salt_restore_under_ice) then 
          do j=jsc,jec
             do i=isc,iec
                if(Grd%tmask(i,j,1) == 1.0) then 
                    if(T_prog(index_temp)%field(i,j,1,tau) <= -0.0539*T_prog(index_salt)%field(i,j,1,tau)) then 
                        open_ocean_mask(i,j) = 0.0
                    endif
                endif
             enddo
          enddo
      endif

      ! salt flux correction is assumed to be a pme_mass field (units kg/(m^2*sec))
      if (use_waterflux) then
          do j = jsc,jec
             do i = isc,iec
                pme_correct(i,j) = open_ocean_mask(i,j)*Grd%tmask(i,j,1)*salt_correction_scale*data(i,j)
             enddo
          enddo     

          ! produce a zero area average so there is no net input 
          ! of water to the ocean associated with the flux correction 
          if(zero_net_water_correction) then 
             pme_correct_total = mpp_global_sum(Dom%domain2d,pme_correct(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                                 global_sum_flag)/Grd%tcellsurf

             ! truncate out the lower order bits if using non_bitwise_reproducible sums
             if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) pme_correct_total = real(pme_correct_total, kind=FLOAT_KIND)
             pme_correct(isc:iec,jsc:jec) = pme_correct(isc:iec,jsc:jec) - pme_correct_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add pme_correct to pme
          pme(isc:iec,jsc:jec) = pme(isc:iec,jsc:jec) + pme_correct(isc:iec,jsc:jec)
          call mpp_update_domains(pme(:,:), Dom%domain2d)

      else  ! salt fluxes 

          ! need to convert the pme correction to an implied salt flux  
          do j = jsc,jec
             do i = isc,iec
                flx_correct(i,j) = open_ocean_mask(i,j)*Grd%tmask(i,j,1)*salt_correction_scale*data(i,j)/(rho0*0.001) 
             enddo
          enddo
         
          ! produce a zero area average so there is no net input
          ! of salt to the ocean associated with the flux correction 
          if(zero_net_salt_correction) then 
             flx_correct_total =  &
             mpp_global_sum(Dom%domain2d,flx_correct(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), global_sum_flag)/Grd%tcellsurf

             ! truncate out the lower order bits if using non_bitwise_reproducible sums
             if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) flx_correct_total = real(flx_correct_total, kind=FLOAT_KIND)
             flx_correct(isc:iec,jsc:jec) = flx_correct(isc:iec,jsc:jec) - flx_correct_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add flx_correct to stf
          do j = jsc,jec
             do i = isc,iec
                T_prog(index_salt)%stf(i,j) = T_prog(index_salt)%stf(i,j) + flx_correct(i,j)
             enddo
          enddo

      endif 

  endif   ! endif for if (id_correction(index_salt) > 0 )


  ! send diagnostics 

  ! salt from restoring 
  if (id_stf_restore(index_salt) > 0) then
      used = send_data(id_stf_restore(index_salt), flx_restore(:,:)*T_prog(index_salt)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                             &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total salt from restoring
  if (id_total_ocean_stf_restore(index_salt) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*flx_restore(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_restore(index_salt), total_stuff*1e-15, Time%model_time)
  endif 

  ! salt from correction 
  if (id_stf_correct(index_salt) > 0) then
      used = send_data(id_stf_correct(index_salt), flx_correct(:,:)*T_prog(index_salt)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                             &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total salt from correction 
  if (id_total_ocean_stf_correct(index_salt) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*flx_correct(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_correct(index_salt), total_stuff*1e-15, Time%model_time)
  endif 

  ! salt from all surface fluxes 
  if (id_stf_total(index_salt) > 0) then
      used = send_data(id_stf_total(index_salt), T_prog(index_salt)%stf(:,:)*T_prog(index_salt)%conversion,&
             Time%model_time, rmask=Grd%tmask(:,:,1),                                                      &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total salt from all fluxes 
  if (id_total_ocean_stf_sum(index_salt) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                     *T_prog(index_salt)%stf(:,:)*T_prog(index_salt)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_sum(index_salt), total_stuff*1e-15, Time%model_time)
  endif

  ! pme from salt restoring
  if (id_pme_restore > 0) then
      used = send_data(id_pme_restore, pme_restore(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total pme from salt restoring
  if(id_total_ocean_pme_restore > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*pme_restore(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_pme_restore, total_stuff*1e-15, Time%model_time)
  endif 

  ! pme from salt correction 
  if (id_pme_correct > 0) then
      used = send_data(id_pme_correct, pme_correct(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1),    &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total pme from salt correction 
  if(id_total_ocean_pme_correct > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*pme_correct(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_pme_correct, total_stuff*1e-15, Time%model_time)
  endif 

  ! pme from all surface terms 
  if (id_pme_net > 0) then
      used = send_data(id_pme_net, pme(:,:),          &
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total pme from all surface terms 
  if(id_total_ocean_pme_net > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*pme(:,:)
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_pme_net, total_stuff*1e-15, Time%model_time)
  endif 

  ! heat input from net pme relative to 0 degrees C (W/m2)
  if (id_stf_pme(index_temp) > 0) then
      used = send_data(id_stf_pme(index_temp),                                    &
             pme(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion, &
             Time%model_time, rmask=Grd%tmask(:,:,1),                             &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  ! total heat flux from net pme (Watts)
  if(id_total_ocean_stf_pme(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                     *pme(:,:)*T_prog(index_temp)%tpme(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_pme(index_temp), total_stuff*1e-15, Time%model_time)
  endif

  if (id_ice_mask > 0) then
      used = send_data(id_ice_mask, (1.0-open_ocean_mask(:,:))*Grd%tmask(:,:,1),&
             Time%model_time,rmask=Grd%tmask(:,:,1),                            &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  if (id_open_ocean_mask > 0) then
      used = send_data(id_open_ocean_mask, open_ocean_mask(:,:),&
             Time%model_time,rmask=Grd%tmask(:,:,1),            &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


  ! initialization of flux fields for temperature 
  flx_correct = 0.0 
  flx_restore = 0.0 

  ! temperature restoring   
  if (id_restore(index_temp) > 0 .and. temp_damp_factor > 0.0) then

     call time_interp_external(id_restore(index_temp), Time%model_time, data)

     if(prog_temp_variable==CONSERVATIVE_TEMP) then 
         do j=jsc,jec
            do i=isc,iec
               flx_restore(i,j) = temp_damp_factor*restore_mask(i,j) &
                                  *(data(i,j)-T_diag(index_diag_temp)%field(i,j,1))
               T_prog(index_temp)%stf(i,j) = T_prog(index_temp)%stf(i,j) + flx_restore(i,j)
            enddo
         enddo
     else 
         do j=jsc,jec
            do i=isc,iec
               flx_restore(i,j) = temp_damp_factor*restore_mask(i,j) &
                                  *(data(i,j)-T_prog(index_temp)%field(i,j,1,taum1))
               T_prog(index_temp)%stf(i,j) = T_prog(index_temp)%stf(i,j) + flx_restore(i,j)
            enddo
         enddo
     endif

  endif ! (id_restore(index_temp) > 0 .and. temp_damp_factor > 0.0)


  ! temperature flux correction 
  if (id_correction(index_temp) > 0 ) then

     call time_interp_external(id_correction(index_temp), Time%model_time, data)

     do j=jsc,jec
        do i=isc,iec
           flx_correct(i,j) = Grd%tmask(i,j,1)*data(i,j)*temp_correction_scale/rho_cp
           T_prog(index_temp)%stf(i,j) = T_prog(index_temp)%stf(i,j) + flx_correct(i,j)
        enddo
     enddo

  endif


  ! restoring heat flux 
  if (id_stf_restore(index_temp) > 0) then
      used = send_data(id_stf_restore(index_temp), flx_restore(:,:)*T_prog(index_temp)%conversion, &
             Time%model_time,rmask=Grd%tmask(:,:,1),                                               &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total of restoring heat flux 
  if(id_total_ocean_stf_restore(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*flx_restore(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_restore(index_temp), total_stuff*1e-15, Time%model_time)
  endif 

  ! flux correction heat flux 
  if (id_stf_correct(index_temp) > 0) then
      used = send_data(id_stf_correct(index_temp), flx_correct(:,:)*T_prog(index_temp)%conversion, &
             Time%model_time,rmask=Grd%tmask(:,:,1),                                               &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total of flux correction heat flux 
  if(id_total_ocean_stf_correct(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:)*flx_correct(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_correct(index_temp), total_stuff*1e-15, Time%model_time)
  endif 

  ! net heat from stf   
  if (id_stf_total(index_temp) > 0) then
      used = send_data(id_stf_total(index_temp),                        &
             T_prog(index_temp)%stf(:,:)*T_prog(index_temp)%conversion, &
             Time%model_time, rmask=Grd%tmask(:,:,1),                   &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  ! total of net heat flux 
  if(id_total_ocean_stf_sum(index_temp) > 0) then 
      wrk1_2d(:,:) = Grd%tmask(:,:,1)*Grd%dat(:,:) &
                     *T_prog(index_temp)%stf(:,:)*T_prog(index_temp)%conversion
      total_stuff  = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:), NON_BITWISE_EXACT_SUM)
      used = send_data (id_total_ocean_stf_sum(index_temp), total_stuff*1e-15, Time%model_time)
  endif 



  ! tau_x flux correction
  if (id_tau_x_correction > 0 ) then
     flx_correct = 0.0
     call time_interp_external(id_tau_x_correction, Time%model_time, data)
     do j=jsc,jec
        do i=isc,iec
           flx_correct(i,j) = tau_x_correction_scale*Grd%umask(i,j,1)*data(i,j)/rho0
           Velocity%smf(i,j,1) = Velocity%smf(i,j,1) + flx_correct(i,j)
        enddo
     enddo
  endif

  ! diagnostics for i-directed wind stress correction (N/m^2)
  if (id_tau_x_flux_correction > 0) then 
       used =  send_data(id_tau_x_flux_correction, flx_correct(:,:)*rho0, &
               Time%model_time, rmask=Grd%umask(:,:,1),                   &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


  ! tau_y correction
  if (id_tau_y_correction > 0 ) then
     flx_correct = 0.0
     call time_interp_external(id_tau_y_correction, Time%model_time, data)
     do j=jsc,jec
        do i=isc,iec
           flx_correct(i,j) = tau_y_correction_scale*Grd%umask(i,j,1)*data(i,j)/rho0
           Velocity%smf(i,j,2) = Velocity%smf(i,j,2) + flx_correct(i,j)
        enddo
     enddo
  endif

  ! diagnostics for j-directed wind stress correction (N/m^2)
  if (id_tau_y_flux_correction > 0) then 
       used =  send_data(id_tau_y_flux_correction, flx_correct(:,:)*rho0, &
               Time%model_time, rmask=Grd%umask(:,:,1),                   &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


  return

end subroutine flux_adjust
! </SUBROUTINE> NAME="flux_adjust"


end module ocean_sbc_mod
