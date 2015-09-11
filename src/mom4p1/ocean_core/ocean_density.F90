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
module ocean_density_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Compute density and related quantities.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module computes the in-situ density and its partial derivatives with 
! respect to conservative temperature or potential temperature, and with 
! respect to salinity.  
!
! Based on Jackett, McDougall, Feistel Wright, and Griffies(2005).  
!
! This equation of state is valid over a "cone-shaped" range 
! corresponding to 
!
! 0psu <= salinity <= 40 psu
!
! -3C <= theta <= 40C    "theta" = either conservative or potential temp	 
!
! 0dbar <= pressure <= 8000dbar 
!
! with the cone getting smaller in the deeper ocean where 
! theta and salinity vary over a smaller range.  
!
!  Input variables are the following:
!
!  --salinity in psu
!  --conservative temperature or potential temperature (theta) in deg C
!  --pressure in dbars  (1bar = 10dbar = 10^5 Newton/m^2 = 10^5 Pascals). 
!
!  Note that in the ocean, pressure increases roughly by 1dbar for each meter depth.
!  Also note that pressure is the "sea pressure", which is the absolute pressure
!  minus the pressure of a standard atmosphere, which is 10.1325 dbars.
!
! check values                                                        <BR/>

!  for "theta" = conservative temperature 
!  rho(s=20psu,theta=20C,p=1000dbar)   = 1017.842890411975 (kg/m^3)   <BR/>
!  alpha(s=20psu,theta=20C,p=1000dbar) = 2.436057013634649e-4 (1/C)   <BR/>
!  beta(s=20psu,theta=20C,p=1000dbar)  = 7.314818108935248e-4 (1/psu) <BR/>
!
!  for "theta" = potential temperature 
!  rho(s=20psu,theta=20C,p=1000dbar)   = 1017.728868019642 (kg/m^3)   <BR/>
!  alpha(s=20psu,theta=20C,p=1000dbar) = 2.525481286927133e-4 (1/C)   <BR/>
!  beta(s=20psu,theta=20C,p=1000dbar)  = 7.379638527217575e-4 (1/psu) <BR/> 
!
! This equation of state should be suitable for purposes of realistic 
! global ocean climate modeling. 
!
!
! B. Linear equation for use in idealized Boussinesq studies
! 
! This equation renders density a linear function of potential 
! temperature and salinity.  All nonlinearities are ignored, as are  
! pressure effects. 
!
! The valid range for theta and salinity arbitrary for the 
! linear equation of state. 
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Feistel (2003), A new extended Gibbs thermodynamic potential 
! of seawater. Progress in Oceanography. vol 58, pages 43-114.
! </REFERENCE>
!
! <REFERENCE>
! Jackett, McDougall, Feistel, Wright, and Griffies (2006)
! Algorithms for density, potential temperature, conservative
! temperature, and freezing temperature of seawater.  
! Journal of Atmospheric and Oceanic Technology, 2006, 
! in press. 
! </REFERENCE>
!
! <REFERENCE>
! McDougall and Jackett (2005)
! The material derivative of neutral density
! Journal of Marine Research, vol 63, pages 159-185.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison,  R.C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, R.C. Pacanowski, R.M. Schmidt, and V. Balaji
! Tracer Conservation with an Explicit Free Surface Method for 
! Z-coordinate Ocean Models
! Monthly Weather Review (2001) vol 129 pages 1081--1098
! </REFERENCE>
!
! <REFERENCE>
! T. McDougall (1987)
! Cabbeling, Thermobaricity, and water mass conversion
! JGR vol 92, pages 5448-5464
! </REFERENCE>
!
! <NOTE>
!
! Density is computed as a function of conservative temperature (degC) 
! or potential temperature (degC), salinity (psu), and in-situ pressure (dbar).
! The pressure contribution includes that from the free surface height 
! and the applied atmospheric and/or sea ice pressure.  
!
! For vert_coordinate==GEOPOTENTIAL, ZSTAR, or ZSIGMA, baroclinic component of
! hydrostatic pressure is not known until the density is known.  In this case,
! the baroclinic pressure contribution to density is lagged by a time step.  
! rho(tau) = rho[theta(tau),s(tau), p_atm(tau) + p_fs(tau) + p_baroclinic(tau-1)].  
! This issue does not arise when using vert_coordinate=PRESSURE, PSTAR, or PSIGMA.
!
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_density_nml">
!
!  <DATA NAME="write_a_restart" TYPE="logical">
!  Set true to write a restart.  False setting only for rare 
!  cases where wish to benchmark model without measuring the cost
!  of writing restarts and associated chksums.  
!  Default is write_a_restart=.true. 
!  </DATA> 
!
!  <DATA NAME="press_standard" UNITS="dbar" TYPE="real">
!  Standard atmospheric pressure (dbar).  The realistic 
!  EOS used in mom4 requires "sea pressure" as an argument
!  rather than absolute pressure.  Sea pressure is 
!  absolute pressure minus a standard atmospheric pressure 
!  of 10.1325dbar.  
!
!  For models that do have a realistic atmospheric loading, then it
!  is appropriate to remove 10.1325dbar prior to computing the EOS.
!  For those cases with zero atmospheric pressure, then it is not
!  necessary to remove the standard atmosphere.  The default for the 
!  press_standard is 0.0dbar.   
!  </DATA> 
!
!  <DATA NAME="t_test" UNITS="C" TYPE="real"> 			
!  Conservative temperature or potential temperature for 
!  testing the EOS.
!  </DATA> 
!  <DATA NAME="s_test" UNITS="psu" TYPE="real">
!  Salinity for testing the EOS.
!  </DATA> 
!  <DATA NAME="p_test" UNITS="dbar" TYPE="real">
!  Sea pressure for testing the EOS.
!  </DATA> 
!
!  <DATA NAME="tn_test" UNITS="C" TYPE="real"> 			
!  Conservative temperature or potential temperature for 
!  testing the equation for neutral density.
!  </DATA> 
!  <DATA NAME="sn_test" UNITS="psu" TYPE="real">
!  Salinity the equation for neutral density.
!  </DATA> 
!
!  <DATA NAME="linear_eos" TYPE="logical">
!  Set to true if wish to use the linear equation of state.  
!  </DATA>
!  <DATA NAME="alpha_linear_eos" TYPE="real">
!  Constant "thermal expansion coefficient" for EOS 
!  rho = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity
!  </DATA>
!  <DATA NAME="beta_linear_eos" TYPE="real">
!  Constant "saline contraction coefficient" for EOS 
!  rho = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity
!  </DATA>
!
!  <DATA NAME="potrho_press" UNITS="dbar" TYPE="real">
!  Sea pressure for computing diagnostic potential density of use 
!  for computing diagnostics with potential density.  
!  </DATA> 
!  <DATA NAME="potrho_min" UNITS="kg/m^3" TYPE="real">
!  Minimum potential density used to partition vertical according 
!  to potential density.  
!  </DATA> 
!  <DATA NAME="potrho_max" UNITS="kg/m^3" TYPE="real">
!  Maximum potential density used to partition vertical according 
!  to potential density.
!  </DATA> 
!
!  <DATA NAME="neutralrho_min" UNITS="kg/m^3" TYPE="real">
!  Minimum neutral density used to partition vertical according 
!  to rational polynomial approximation to neutral density.  
!  </DATA> 
!  <DATA NAME="neutralrho_max" UNITS="kg/m^3" TYPE="real">
!  Maximum neutral density used to partition vertical according 
!  to rational polynomial approximation to neutral density.  
!  </DATA> 
!
!  <DATA NAME="theta_min" UNITS="C" TYPE="real">
!  Minimum conservative temperature or potential temperature used to 
!  partition vertical according to temperature.  
!  </DATA> 
!  <DATA NAME="theta_max" UNITS="C" TYPE="real">
!  Maximum conservative temperature or potential temperature used to
!  partition vertical according to temperature. 
!  </DATA> 
!
!  <DATA NAME="layer_nk" TYPE="integer">
!  Number of classes used to partition vertical according to potential density,
!  conservative temperature, or potential temperature. Used for diagnostics. 
!  </DATA> 
!
!  <DATA NAME="buoyfreq_smooth_vert" TYPE="logical">
!  To smooth the vertical temp and salt derivative for diagnosing 
!  the buoyancy frequency. Default buoyfreq_smooth_vert=.true.
!  </DATA>
!
!  <DATA NAME="epsln_drhodz" UNITS="kg/m4" TYPE="real">
!  To normalize the inverse vertical derivative of neutral density 
!  for computing the buoyancy frequency. Default epsln_drhodz=1e-10.
!  </DATA>
!
!  <DATA NAME="mask_domain_restart" TYPE="logical">
!  For cases where use the domain masking, it is necessary to initialize the field 
!  denominator_r to nonzero in order to avoid NaNs in the case when change processor
!  layout in between restarts.  Note that when use solid wall boundary conditions, 
!  this logical should remain false in order to bitwise reproduce across restarts.
!  Default mask_domain_restart=.false. 
!  </DATA>
!
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For debugging nonlinear equation of state 
!  </DATA>
!  <DATA NAME="rho0_density" TYPE="logical">
!  For debugging, it is often useful to have rho=rho0 uniform.
!  </DATA>
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase 
!  efficiency.
!  default: do_bitwise_exact_sum=.false.
!  </DATA>
!
!</NAMELIST>

#include <fms_platform.h>

use constants_mod,       only: epsln, grav, c2dbars
use diag_manager_mod,    only: register_diag_field, register_static_field, diag_axis_init
use diag_manager_mod,    only: need_data, send_data
use fms_mod,             only: write_version_number, mpp_error
use fms_mod,             only: field_exist, FATAL, WARNING
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, file_exist
use fms_io_mod,          only: register_restart_field, save_restart, restore_state
use fms_io_mod,          only: reset_field_pointer, restart_file_type
use mpp_domains_mod,     only: mpp_update_domains, mpp_global_sum
use mpp_domains_mod,     only: BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use mpp_mod,             only: stdout, stdlog, mpp_chksum
use platform_mod,        only: i8_kind
use time_manager_mod,    only: time_type, increment_time
use field_manager_mod,   only: fm_get_index

use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: GEOPOTENTIAL, ZSTAR, ZSIGMA, PRESSURE, PSTAR, PSIGMA 
use ocean_parameters_mod, only: CONSERVATIVE_TEMP, POTENTIAL_TEMP
use ocean_parameters_mod, only: missing_value, onefourth, onehalf, rho0r, rho0
use ocean_pressure_mod,   only: pressure_in_dbars
use ocean_types_mod,      only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,      only: ocean_time_type, ocean_time_steps_type, ocean_options_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_density_type, ocean_external_mode_type
use ocean_types_mod,      only: ocean_diag_tracer_type
use ocean_util_mod,       only: write_timestamp
use ocean_workspace_mod,  only: wrk1, wrk2, wrk3, wrk4, wrk1_2d, wrk2_2d, wrk3_2d, wrk4_2d

implicit none

private

real, dimension(:), allocatable     :: a, b                   ! polynomial coefficients in the realistic EOS
real, dimension(:,:), allocatable   :: rhodz_tau              ! vertical integral of rho*dzt at time tau
real, dimension(:,:), allocatable   :: rhodz_taup1            ! vertical integral of rho*dzt at time taup1
real, dimension(:,:,:), allocatable :: denominator_r          ! reciprocal of denominator for relastic EOS
logical                             :: linear_eos=.false.     ! for use of ideal linear equation of state      
integer                             :: vert_coordinate        ! for setting which vertical coordinate to use
integer                             :: temperature_variable   ! for whether using conservative or potential temp
integer                             :: global_sum_flag        ! for determing type of global sum BITWISE/NON

! polynomial coefficients in the realistic EOS
real :: a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11      
real :: b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12

! multiplied constants for density partial derivs
real :: two_a2, three_a3, two_a6, six_a3, two_a8, two_a10, two_a11, four_a11             
real :: two_b2, three_b3, six_b3, four_b4, twelve_b4, three_b7, six_b7             
real :: onep5_b8, onep5_b9, two_b9, three_b9, two_b11, three_b11, six_b11, three_b12   

! polynomial coefficients in the neutral density 
real :: a0n, a1n, a2n, a3n, a4n, a5n, a6n
real :: b0n, b1n, b2n, b3n, b4n, b5n, b6n, b7n, b8n, b9n

! some test settings for density 
real :: s_test=20.0
real :: t_test=20.0
real :: p_test=1000.0
real :: jmfwg_rho 
real :: jmfwg_alpha 
real :: jmfwg_beta 

! some test settings for neutral density 
real :: sn_test=35.0
real :: tn_test=20.0
real :: mj_neutralrho 

! for linear EOS 
real :: alpha_linear_eos=0.255
real :: beta_linear_eos =0.0

! time steps 
real :: dtts =0.0

! inverse area (1/m) of k=1 tracer cells 
real :: area_total_r  

! inverse gravity constant 
real :: grav_r 

! for diagnostics manager 
logical :: used
integer :: id_drhodtheta     =-1
integer :: id_drhodsalt      =-1
integer :: id_drhodpress     =-1
integer :: id_sound_speed2   =-1
integer :: id_press          =-1
integer :: id_rho            =-1
integer :: id_rho0           =-1
integer :: id_neutral_rho    =-1
integer :: id_pot_rho        =-1
integer :: id_pot_rho_0      =-1
integer :: id_pot_rho_1      =-1
integer :: id_pot_rho_2      =-1
integer :: id_pot_rho_3      =-1
integer :: id_pot_rho_4      =-1
integer :: id_pot_rho_et     =-1
integer :: id_pot_rho_nt     =-1
integer :: id_pot_rho_wt     =-1
integer :: id_rho_average    =-1
integer :: id_mass_level     =-1
integer :: id_volume_level   =-1
integer :: id_drhodz_zt      =-1  
integer :: id_drhodz_wt      =-1  
integer :: id_buoyfreq2_zt   =-1  
integer :: id_buoyfreq2_wt   =-1  
integer :: id_cabbeling      =-1
integer :: id_thermobaricity =-1
integer :: id_int_rhodz      =-1

integer :: id_rhoave              =-1
integer :: id_pbot_adjust         =-1
integer :: id_eta_adjust          =-1
integer :: id_eta_adjust_approx   =-1
integer :: id_eta_adjust_cstvolume=-1


!for restart
integer                       :: id_restart_rho = 0
type(restart_file_type), save :: Den_restart

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

character(len=128) :: version=&
       '$Id: ocean_density.F90,v 16.0.2.6.4.1.2.2.10.5.6.1.4.3 2009/12/04 14:39:46 smg Exp $'
character (len=128) :: tagname = &
       '$Name: mom4p1_pubrel_dec2009_nnz $'

public ocean_density_init
public ocean_density_end
public ocean_density_diag
public update_ocean_density
public density
public potential_density
public neutral_density
public density_sfc
public density_derivs
public density_delta_z
public density_delta_sfc
public ocean_density_restart

private ocean_density_chksum
private compute_buoyfreq
private cabbeling_thermobaricity

! interfaces for density 
interface density
   module procedure density_field
   module procedure density_level
   module procedure density_line
   module procedure density_point
end interface

! interfaces for neutral density 
interface neutral_density
   module procedure neutral_density_field
   module procedure neutral_density_point
end interface


! interfaces for density derivatives 
interface density_derivs
   module procedure density_derivs_field
   module procedure density_derivs_level
   module procedure density_derivs_point
end interface

integer :: index_temp=-1
integer :: index_salt=-1
integer :: num_prog_tracers

! indices for diagnostic tracer density elements
integer :: ind_rho = -1
integer :: ind_neutralrho = -1
integer :: ind_potrho_0 = -1
integer :: ind_potrho_1 = -1
integer :: ind_potrho_2 = -1
integer :: ind_potrho_3 = -1
integer :: ind_potrho_4 = -1

! nml settings 

! for diagnostic partitioning of vertical according to 
! potential density, neutral density, or theta classes 
integer :: layer_nk = 80       ! # of classes used to compute diagnostics associated with 
                               ! potential density, potential temperature, or conservative
                               ! temperature classes. 

! for diagnostic partitioning of vertical according to potential density classes 
real :: potrho_press = 2000.0  ! sea pressure (dbars) for potential density computed for diagnostics
real :: potrho_min   = 1028.0  ! (kg/m^3)         
real :: potrho_max   = 1038.0  ! (kg/m^3)           

! for diagnostic partitioning of vertical according to 
! rational polynomial approximation to neutral density
real :: neutralrho_min = 1020.0  ! (kg/m^3)         
real :: neutralrho_max = 1030.0  ! (kg/m^3)           

! for diagnostic partitioning of vertical according 
! to potential temperature or conservative temperature classes 
real :: theta_min = -2.0  ! (degrees C)         
real :: theta_max = 30.0  ! (degrees C)           

! standard atmospheric pressure (dbar).
! Should be set to 0.0 if assume zero pressure for the 
! overlying atmosphere.  But if have a realistic atmospheric
! pressure loading, then set press_standard=10.1325.  
real :: press_standard = 0.0 

! for smoothing the vertical density gradient 
! when diagnosing the buoyancy frequency. 
logical :: buoyfreq_smooth_vert=.true.
integer :: num_121_passes      =1
real    :: epsln_drhodz        =1.e-10

! for case when use domain masking, 
logical :: mask_domain_restart = .false. 

logical :: debug_this_module      = .false. 
logical :: rho0_density           = .false.
logical :: module_is_initialized  = .false.
logical :: write_a_restart        = .true. 
logical :: do_bitwise_exact_sum   = .false.

namelist /ocean_density_nml/ s_test, t_test, p_test, press_standard,             &
                             sn_test, tn_test,                                   &
                             linear_eos, alpha_linear_eos, beta_linear_eos,      & 
                             potrho_press, potrho_min, potrho_max,               &
                             neutralrho_min, neutralrho_max,                     &
                             layer_nk, theta_min, theta_max,                     &
                             debug_this_module, rho0_density, write_a_restart,   & 
                             buoyfreq_smooth_vert, num_121_passes, epsln_drhodz, &
                             mask_domain_restart, do_bitwise_exact_sum
contains

!#######################################################################
! <SUBROUTINE NAME="ocean_density_init">
!
! <DESCRIPTION>
! Initialize the density module
! </DESCRIPTION>
!
  subroutine ocean_density_init (Grid, Domain, Time, Time_steps, Thickness, T_prog, T_diag, &
                                 Ocean_options, Dens, ver_coordinate, debug)

    type(ocean_grid_type),   target, intent(in)    :: Grid
    type(ocean_domain_type), target, intent(in)    :: Domain
    type(ocean_time_type),           intent(in)    :: Time
    type(ocean_time_steps_type),     intent(in)    :: Time_steps
    type(ocean_thickness_type),      intent(in)    :: Thickness
    type(ocean_prog_tracer_type),    intent(in)    :: T_prog(:)
    type(ocean_diag_tracer_type),    intent(inout) :: T_diag(:)
    type(ocean_options_type),        intent(inout) :: Ocean_options
    type(ocean_density_type),        intent(inout) :: Dens

    integer,          intent(in)            :: ver_coordinate
    logical,          intent(in), optional  :: debug
    
    integer :: ioun, io_status, ierr
    integer :: k, n
    integer :: tau, taup1, taum1
    real    :: neutralrho_test, rho_test, alpha_test, beta_test, speed2_test
    real    :: drho_dtheta_test, drho_dsal_test, drho_dpress_test
    real    :: diff_test
    real    :: rho0_nk(nk)

    real, allocatable, dimension(:) :: rho_average
    real, allocatable, dimension(:) :: mass_level
    real, allocatable, dimension(:) :: volume_level 

    integer :: id_potrho_bounds, id_potrho_axis
    integer :: id_potrho_xt, id_potrho_yt

    integer :: id_neutralrho_bounds, id_neutralrho_axis
    integer :: id_neutralrho_xt, id_neutralrho_yt

    integer :: id_theta_bounds, id_theta_axis
    integer :: id_theta_xt, id_theta_yt
    integer :: id_restart(5)

    real, allocatable, dimension(:) :: potrho_bounds
    real, allocatable, dimension(:) :: neutralrho_bounds
    real, allocatable, dimension(:) :: theta_bounds

    real    :: potrho_interval 
    real    :: neutralrho_interval 
    real    :: theta_interval 
    integer :: temp_variable
    character*128 filename

    integer :: stdoutunit,stdlogunit 
    stdoutunit=stdout();stdlogunit=stdlog() 


    if ( module_is_initialized ) then 
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod (ocean_density_init): module already initialized.')
    endif 

    module_is_initialized = .TRUE.

    tau             = Time%tau
    taup1           = Time%taup1
    taum1           = Time%taum1
    vert_coordinate = ver_coordinate
 
    num_prog_tracers = size(T_prog(:))
    do n=1,num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
       if (T_prog(n)%longname == 'Conservative temperature') temp_variable = CONSERVATIVE_TEMP
       if (T_prog(n)%longname == 'Potential temperature')    temp_variable = POTENTIAL_TEMP
    enddo

    call write_version_number( version, tagname )

    ! provide for namelist override of defaults
    ioun = open_namelist_file()
    read (ioun,ocean_density_nml,IOSTAT=io_status)
    write (stdoutunit,'(/)')
    write (stdoutunit,ocean_density_nml)
    write (stdlogunit,ocean_density_nml)
    ierr = check_nml_error(io_status, 'ocean_density_nml')
    call close_file(ioun)

    if (PRESENT(debug) .and. .not. debug_this_module) then
      debug_this_module = debug
    endif 
    if(debug_this_module) then 
      write(stdoutunit,'(a)') '==>Note: running ocean_density with debug_this_module=.true.'  
    endif 

    if(.not. write_a_restart) then 
      write(stdoutunit,'(a)') '==>Note: running ocean_density with write_a_restart=.false.'
      write(stdoutunit,'(a)') '   Will NOT write restart file, and so cannot restart the run.'
    endif 

    if(do_bitwise_exact_sum) then
       global_sum_flag = BITWISE_EXACT_SUM
    else
       global_sum_flag = NON_BITWISE_EXACT_SUM
    endif

    if(linear_eos) then
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: USING linear EOS designed for idealized Boussinesq simulations.'
        write (stdoutunit,'(7x,a)') &
        'It is a linear function of potential temperature and salinity.'
        write (stdoutunit,'(7x,a)') &
        'There is no pressure dependence.'
        Ocean_options%equation_of_state = 'Linear equation of state used for density.'
    else
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: USING full EOS, as relevant for realistic ocean climate simulations.' 
        write (stdoutunit,'(1x,a,f12.6,a)') &
        ' Subtracting standard atmosphere of ',press_standard,' dbar for EOS calculation.'  
        Ocean_options%equation_of_state = 'Realistic equation of state used for density.'
    endif

    if(temp_variable==CONSERVATIVE_TEMP) then
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing EOS assuming mom4p1 temp = conservative temperature.'
    elseif(temp_variable==POTENTIAL_TEMP) then
        write (stdoutunit,'(/1x,a)') &
        ' ==> Note: Computing EOS assuming mom4p1 temp = potential temperature.'
    else
      call mpp_error(FATAL, &
       '==>Error in ocean_density_mod: model temperature variable remains unspecified')
    endif 

    write(stdoutunit,*) &
    '==>Note: The Boussinesq rho0 density has a value of (kg/m^3) ',rho0

    if(rho0_density) then 
        write(stdoutunit,*) &
        '==>Warning: rho0_density=.true, so rho=rho0 everywhere. Is this what you wish to do?'
    endif


    ! 25 coefficients in the realistic equation of state 
    if(temp_variable==CONSERVATIVE_TEMP) then

        jmfwg_rho   = 1017.842890411975d0
        jmfwg_alpha = 2.436057013634649d-4
        jmfwg_beta  = 7.314818108935248d-4

        a0  =  9.9983912878771446d+02
        a1  =  7.0687133522652896d+00
        a2  = -2.2746841916232965d-02
        a3  =  5.6569114861400121d-04
        a4  =  2.3849975952593345d+00
        a5  =  3.1761924314867009d-04
        a6  =  1.7459053010547962d-03
        a7  =  1.2192536310173776d-02
        a8  =  2.4643435731663949d-07
        a9  =  4.0525405332794888d-06
        a10 = -2.3890831309113187d-08
        a11 = -5.9016182471196891d-12

        b0  =  1.0000000000000000d+00 
        b1  =  7.0051665739672298d-03
        b2  = -1.5040804107377016d-05 
        b3  =  5.3943915288426715d-07
        b4  =  3.3811600427083414d-10
        b5  =  1.5599507046153769d-03
        b6  = -1.8137352466500517d-06
        b7  = -3.3580158763335367d-10
        b8  =  5.7149997597561099d-06
        b9  =  7.8025873978107375d-10
        b10 =  7.1038052872522844d-06
        b11 = -2.1692301739460094d-17
        b12 = -8.2564080016458560d-18

        mj_neutralrho=1024.43863927763d0
        
        a0n =  1.0022048243661291d+003
        a1n =  2.0634684367767725d-001
        a2n =  8.0483030880783291d-002
        a3n = -3.6670094757260206d-004
        a4n = -1.4602011474139313d-003
        a5n = -2.5860953752447594d-003
        a6n = -3.0498135030851449d-007

        b0n =  1.0000000000000000d+000 
        b1n =  4.4946117492521496d-005
        b2n =  7.9275128750339643d-005
        b3n = -1.2358702241599250d-007
        b4n = -4.1775515358142458d-009
        b5n = -4.3024523119324234d-004
        b6n =  6.3377762448794933d-006
        b7n = -7.2640466666916413d-010
        b8n = -5.1075068249838284d-005
        b9n = -5.8104725917890170d-009


    elseif(temp_variable==POTENTIAL_TEMP) then 
    
        jmfwg_rho   =  1017.728868019642d0
        jmfwg_alpha = 2.525481286927133d-4
        jmfwg_beta  = 7.379638527217575d-4

        a0  =  9.9984085444849347d+02
        a1  =  7.3471625860981584d+00
        a2  = -5.3211231792841769d-02
        a3  =  3.6492439109814549d-04
        a4  =  2.5880571023991390d+00
        a5  = -6.7168282786692355d-03
        a6  =  1.9203202055760151d-03
        a7  =  1.1798263740430364d-02
        a8  =  9.8920219266399117d-08
        a9  =  4.6996642771754730d-06
        a10 = -2.5862187075154352d-08
        a11 = -3.2921414007960662d-12

        b0  =  1.0000000000000000d+00 
        b1  =  7.2815210113327091d-03
        b2  = -4.4787265461983921d-05 
        b3  =  3.3851002965802430d-07
        b4  =  1.3651202389758572d-10
        b5  =  1.7632126669040377d-03
        b6  = -8.8066583251206474d-06
        b7  = -1.8832689434804897d-10
        b8  =  5.7463776745432097d-06
        b9  =  1.4716275472242334d-09
        b10 =  6.7103246285651894d-06
        b11 = -2.4461698007024582d-17
        b12 = -9.1534417604289062d-18

        mj_neutralrho=1024.59416751197d0

        a0n =  1.0023063688892480d+003
        a1n =  2.2280832068441331d-001
        a2n =  8.1157118782170051d-002
        a3n = -4.3159255086706703d-004
        a4n = -1.0304537539692924d-004
        a5n = -3.1710675488863952d-003
        a6n = -1.7052298331414675d-007

        b0n =  1.0000000000000000d+000 
        b1n =  4.3907692647825900d-005
        b2n =  7.8717799560577725d-005
        b3n = -1.6212552470310961d-007
        b4n = -2.3850178558212048d-009
        b5n = -5.1268124398160734d-004
        b6n =  6.0399864718597388d-006
        b7n = -2.2744455733317707d-009
        b8n = -3.6138532339703262d-005
        b9n = -1.3409379420216683d-009
       
    endif 

    ! save some multiples of the coefficients 
    two_a2   = 2.0*a2
    three_a3 = 3.0*a3
    six_a3   = 6.0*a3
    two_a6   = 2.0*a6
    two_a8   = 2.0*a8
    two_a10  = 2.0*a10
    two_a11  = 2.0*a11
    four_a11 = 4.0*a11

    two_b2    = 2.0*b2
    three_b3  = 3.0*b3
    six_b3    = 6.0*b3
    four_b4   = 4.0*b4
    twelve_b4 = 12.0*b4
    three_b7  = 3.0*b7
    six_b7    = 6.0*b7
    onep5_b8  = 1.5*b8
    onep5_b9  = 1.5*b9
    two_b9    = 2.0*b9
    three_b9  = 3.0*b9
    two_b11   = 2.0*b11
    three_b11 = 3.0*b11
    six_b11   = 6.0*b11
    three_b12 = 3.0*b12


#ifndef MOM4_STATIC_ARRAYS
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
    allocate(Dens%rho(isd:ied,jsd:jed,nk,3))
    allocate(Dens%rho_fresh(isd:ied,jsd:jed,nk))
    allocate(Dens%potrho(isd:ied,jsd:jed,nk))
    allocate(Dens%neutralrho(isd:ied,jsd:jed,nk))
    allocate(Dens%pressure_at_depth(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodT(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodS(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodP(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodz_wt(isd:ied,jsd:jed,nk))
    allocate(Dens%drhodz_zt(isd:ied,jsd:jed,nk))
    allocate(Dens%dTdz_zt(isd:ied,jsd:jed,nk))
    allocate(Dens%dSdz_zt(isd:ied,jsd:jed,nk))
#endif

    do k=1,nk 
       Dens%rho_fresh(:,:,k)         = 1000.0*Grid%tmask(:,:,k)
       Dens%potrho(:,:,k)            = rho0*Grid%tmask(:,:,k)
       Dens%neutralrho(:,:,k)        = rho0*Grid%tmask(:,:,k)
       Dens%pressure_at_depth(:,:,k) = rho0*grav*Thickness%depth_zt(:,:,k)*c2dbars
    enddo
    Dens%rho(:,:,:,1)    = Grid%tmask(:,:,:)*rho0 
    Dens%rho(:,:,:,2)    = Grid%tmask(:,:,:)*rho0 
    Dens%rho(:,:,:,3)    = Grid%tmask(:,:,:)*rho0 
    Dens%drhodT(:,:,:)   = 0.0
    Dens%drhodS(:,:,:)   = 0.0
    Dens%drhodP(:,:,:)   = 0.0
    Dens%drhodz_wt(:,:,:)= 0.0
    Dens%drhodz_zt(:,:,:)= 0.0
    Dens%dTdz_zt(:,:,:)  = 0.0
    Dens%dSdz_zt(:,:,:)  = 0.0

    Dom  => Domain
    Grd  => Grid
    dtts = Time_steps%dtts

    ! for diagnostics    
    area_total_r = mpp_global_sum(Dom%domain2d,Grd%dat(:,:)*Grd%tmask(:,:,1), global_sum_flag)
    if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) area_total_r = real(area_total_r,kind=FLOAT_KIND)
    area_total_r = 1.0/(epsln+area_total_r) 

    grav_r = 1.0/grav
    allocate(rhodz_tau(isd:ied,jsd:jed))
    allocate(rhodz_taup1(isd:ied,jsd:jed))
    rhodz_tau   = 0.0
    rhodz_taup1 = 0.0


   ! send as 1-dim field, due to problems with diag manager 
   id_rho0 = register_static_field('ocean_model','rho0', Grid%tracer_axes(3:3), &
             'reference density for Boussinesq approximation','kg/m^3',         &
             missing_value=missing_value, range=(/-1e0,1e10/),                  &
             standard_name='reference_sea_water_density_for_boussinesq_approximation')
   rho0_nk(:) = rho0 
   if (id_rho0 > 0) then 
      used = send_data(id_rho0, rho0_nk(:), Time%model_time)
   endif 


    allocate(denominator_r(isd:ied,jsd:jed,nk))
    denominator_r(:,:,:) = 0.0

    filename = 'ocean_density.res.nc'    
    id_restart_rho = register_restart_field(Den_restart, filename, 'rho', Dens%rho(:,:,:,taup1), &
                    domain=Dom%domain2d )
    id_restart(1)  = register_restart_field(Den_restart, filename, 'pressure_at_depth', Dens%pressure_at_depth(:,:,:), &
                    domain=Dom%domain2d )
    id_restart(2)  = register_restart_field(Den_restart, filename, 'denominator_r', denominator_r(:,:,:), &
                    domain=Dom%domain2d )
    id_restart(3)  = register_restart_field(Den_restart, filename, 'drhodT', Dens%drhodT(:,:,:), &
                    domain=Dom%domain2d )
    id_restart(4)  = register_restart_field(Den_restart, filename, 'drhodS', Dens%drhodS(:,:,:), &
                    domain=Dom%domain2d )
    id_restart(5)  = register_restart_field(Den_restart, filename, 'drhodz_zt', Dens%drhodz_zt(:,:,:), &
                    domain=Dom%domain2d )

    filename = 'INPUT/ocean_density.res.nc'
    if (.NOT.file_exist(trim(filename)) )then
       Dens%rho(:,:,:,taum1) = density(T_prog(index_salt)%field(:,:,:,tau),&
                                       T_prog(index_temp)%field(:,:,:,tau),&
                                       Dens%pressure_at_depth(:,:,:))
       Dens%rho(:,:,:,tau)   = Dens%rho(:,:,:,taum1)
       Dens%rho(:,:,:,taup1) = Dens%rho(:,:,:,taum1)
    else
       write (stdoutunit,'(/a)') '  Reading restart for density from INPUT/ocean_density.res.nc'
       call restore_state( Den_restart, id_restart_rho )
       call mpp_update_domains(Dens%rho(:,:,:,taup1),        Dom%domain2d)
       call restore_state( Den_restart, id_restart(1) )
       call mpp_update_domains(Dens%pressure_at_depth(:,:,:),Dom%domain2d)
       call restore_state( Den_restart, id_restart(2) )

       ! initialize to epsln to eliminate NaNs when mask 
       ! out processors and then change layout upon a restart 
       if(mask_domain_restart) then 
          where (denominator_r==0.)  denominator_r=rho0r
       endif 
       call mpp_update_domains(denominator_r(:,:,:),Dom%domain2d)

       ! determine whether to read density derivs at the initial time step of integration.
       ! early versions of MOM4p1 did not contain drhodT and drhodS in the restart file, so 
       ! to allow older restart files to be used with newer releases of MOM4p1, we 
       ! provide for the possibility that the restarts do not contain drhodT and drhodS.
       if(field_exist(filename, 'drhodT')) then
           call restore_state( Den_restart, id_restart(3) )
           call restore_state( Den_restart, id_restart(4) )
           call mpp_update_domains(Dens%drhodT(:,:,:),Dom%domain2d)
           call mpp_update_domains(Dens%drhodS(:,:,:),Dom%domain2d)
       else 
           write(stdoutunit,'(a)') '==>Note: ocean_density_mod: did not read density derivatives from restart.' 
       endif 

   endif

   ! compute buoyancy frequency for diagnostic purposes 
   call compute_buoyfreq(Time, Thickness, T_prog(index_salt)%field(:,:,:,tau), &
                         T_prog(index_temp)%field(:,:,:,tau), Dens) 

       ! over-write drhodz_zt from compute_buoyfreq with drhodz_zt from restart. 
       ! Note that we determine whether to read drhodz_zt at the initial time 
       ! step of integration. early versions of MOM4p1 did not contain drhodz_zt 
       ! in the restart file, so to allow older restart files to be used with 
       ! newer releases of MOM4p1, we provide for the possibility that the 
       ! restarts do not contain drhodz_zt.
   if (file_exist(trim(filename)) )then 
       if(field_exist(filename, 'drhodz_zt')) then
           call restore_state( Den_restart, id_restart(5) )
          call mpp_update_domains(Dens%drhodz_zt(:,:,:),Dom%domain2d)
       else 
          write(stdoutunit,'(a)') '==>Note: ocean_density_mod: did not read drhodz_zt from restart.' 
       endif
   endif
   ! compute neutral density for diagnostic purposes 
   Dens%neutralrho(:,:,:) = Grd%tmask(:,:,:)             &
   *neutral_density(T_prog(index_salt)%field(:,:,:,tau), &
                    T_prog(index_temp)%field(:,:,:,tau))
   ind_neutralrho = fm_get_index('/ocean_mod/diag_tracers/rho_neutral')
   if (ind_neutralrho .gt. 0) then
     T_diag(ind_neutralrho)%field(:,:,:) = Dens%neutralrho(:,:,:)
   endif

   ! compute potential density for diagnostic purposes 
   Dens%potrho(:,:,:) = Grd%tmask(:,:,:)                   &
   *potential_density(T_prog(index_salt)%field(:,:,:,tau), &
                      T_prog(index_temp)%field(:,:,:,tau), potrho_press)

   ! compute potential densities for diagnostic purposes 
   ind_potrho_0 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_0')
   if (ind_potrho_0 .gt. 0) then
     T_diag(ind_potrho_0)%field(:,:,:) = Grd%tmask(:,:,:)*                &
          potential_density(T_prog(index_salt)%field(:,:,:,tau),          &
                            T_prog(index_temp)%field(:,:,:,tau), 0.0)
   endif
   ind_potrho_1 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_1')
   if (ind_potrho_1 .gt. 0) then
     T_diag(ind_potrho_1)%field(:,:,:) = Grd%tmask(:,:,:)*                &
          potential_density(T_prog(index_salt)%field(:,:,:,tau),          &
                            T_prog(index_temp)%field(:,:,:,tau), 1000.0)
   endif
   ind_potrho_2 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_2')
   if (ind_potrho_2 .gt. 0) then
     T_diag(ind_potrho_2)%field(:,:,:) = Grd%tmask(:,:,:)*                &
          potential_density(T_prog(index_salt)%field(:,:,:,tau),          &
                            T_prog(index_temp)%field(:,:,:,tau), 2000.0)
   endif
   ind_potrho_3 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_3')
   if (ind_potrho_3 .gt. 0) then
     T_diag(ind_potrho_3)%field(:,:,:) = Grd%tmask(:,:,:)*                &
          potential_density(T_prog(index_salt)%field(:,:,:,tau),          &
                            T_prog(index_temp)%field(:,:,:,tau), 3000.0)
   endif
   ind_potrho_4 = fm_get_index('/ocean_mod/diag_tracers/rho_pot_4')
   if (ind_potrho_4 .gt. 0) then
     T_diag(ind_potrho_4)%field(:,:,:) = Grd%tmask(:,:,:)*                &
          potential_density(T_prog(index_salt)%field(:,:,:,tau),          &
                            T_prog(index_temp)%field(:,:,:,tau), 4000.0)
   endif

   if(rho0_density) then 
       Dens%rho(:,:,:,taum1) = Grd%tmask(:,:,:)*rho0
       Dens%rho(:,:,:,tau)   = Grd%tmask(:,:,:)*rho0
       Dens%rho(:,:,:,taup1) = Grd%tmask(:,:,:)*rho0
   endif
   ind_rho = fm_get_index('/ocean_mod/diag_tracers/rho_in_situ')
   if (ind_rho .gt. 0) then
     T_diag(ind_rho)%field(:,:,:) = Dens%rho(:,:,:,tau)
   endif

   ! compute initial level mass, volume, and density 
   allocate(rho_average(0:nk))
   allocate(mass_level(0:nk))
   allocate(volume_level(0:nk))
   rho_average(:)  = 0.0
   mass_level(:)   = 0.0
   volume_level(:) = 0.0
   wrk1_2d(:,:)    = 0.0
   wrk2_2d(:,:)    = 0.0
   wrk3_2d(:,:)    = 0.0
   do k=1,nk
      wrk1_2d(:,:)    = Grd%dat(:,:)*Grd%tmask(:,:,k)  
      wrk2_2d(:,:)    = wrk1_2d(:,:)*Thickness%dzt(:,:,k)
      wrk3_2d(:,:)    = wrk1_2d(:,:)*Thickness%dzt(:,:,k)*Dens%rho(:,:,k,tau)
      volume_level(k) = mpp_global_sum(Dom%domain2d,wrk2_2d(:,:),global_sum_flag)
      mass_level(k)   = mpp_global_sum(Dom%domain2d,wrk3_2d(:,:),global_sum_flag)
      if (global_sum_flag .eq. NON_BITWISE_EXACT_SUM) then
        volume_level(k) = real(volume_level(k),kind=FLOAT_KIND)
        mass_level(k)   = real(mass_level(k),kind=FLOAT_KIND)
      endif
      if(volume_level(k) > 0.0) rho_average(k) = mass_level(k)/volume_level(k)
   enddo
   do k=1,nk
      volume_level(0) = volume_level(0) + volume_level(k)
      mass_level(0)   = mass_level(0)   + mass_level(k)
   enddo
   if(volume_level(0) > 0.0) rho_average(0) = mass_level(0)/volume_level(0)

   write(stdoutunit,'(/a,e24.12)') &
   ' ==>Note: From ocean_density_mod: Boussinesq density rho0(kg/m3) = ',rho0
   write(stdoutunit,'(a,e24.12)')  &
   '          Initial rho_average(kg/m3) = ',rho_average(0)
   if(rho_average(0) /= rho0 .and. Time%init) then 
     write(stdoutunit,'(a)')  &
     '          Since rho0 .ne. rho_average, consider changing rho0 in' 
     write(stdoutunit,'(a/)') &
     '          ocean_paramters.h to be equal to rho_average for better accuracy.'
   endif 
 
   id_rho_average = register_static_field('ocean_model','rho_average', Grid%tracer_axes(3:3), &
                   'level average density at initial step of segment','kg/m^3',               &
                   missing_value=missing_value, range=(/-1e0,1e10/))
   id_mass_level  = register_static_field('ocean_model','mass_level', Grid%tracer_axes(3:3), &
                   'mass of levels at initial step of segment','kg',                         &
                   missing_value=missing_value, range=(/-1e0,1e10/))
   id_volume_level= register_static_field('ocean_model','volume_level', Grid%tracer_axes(3:3), &
                   'volume of levels at initial step of segment','m^3',                        &
                   missing_value=missing_value, range=(/-1e0,1e10/))

   if (id_rho_average > 0) then 
      used = send_data(id_rho_average, rho_average(1:nk), Time%model_time)
   endif 
   if (id_mass_level > 0) then 
      used = send_data(id_mass_level, mass_level(1:nk), Time%model_time)
   endif 
   if (id_volume_level > 0) then 
      used = send_data(id_volume_level, volume_level(1:nk), Time%model_time)
   endif 

 
   ! Test values for EOS 
   if(.not. linear_eos) then 

      write(stdoutunit,'(/a)') 'From ocean_density_mod: initial density chksums'
      call write_timestamp(Time%model_time)
      call ocean_density_chksum(Time, Dens)

      write (stdoutunit,'(/,a)') 'EQUATION OF STATE TEST VALUES'

      write (stdoutunit,'(a,f6.2,a,f6.2,a,f8.2)') &
      's_test(psu) = ',s_test,', t_test(C) = ',t_test,', p_test(dbar) = ',p_test  

      rho_test = density(s_test,t_test,p_test) 
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'rho  (',s_test,',',t_test,',',p_test,') = ',rho_test,' kg/m^3'

      diff_test = rho_test-jmfwg_rho
      write (stdoutunit,' (a,e22.16,a)') 'diff from JMFWG = ',diff_test,' kg/m^3'

      call density_derivs(rho_test, s_test, t_test, p_test, &
                          drho_dtheta_test, drho_dsal_test, drho_dpress_test)

      alpha_test = -drho_dtheta_test/(epsln+rho_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'alpha(',s_test,',',t_test,',',p_test,') = ',alpha_test,' 1/C'

      diff_test = alpha_test-jmfwg_alpha
      write (stdoutunit,' (a,e22.16,a)') 'diff from JMFWG = ',diff_test,' 1/C'

      beta_test  =  drho_dsal_test  /(epsln+rho_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'beta (',s_test,',',t_test,',',p_test,') = ',beta_test,' 1/psu'

      speed2_test  =  1.0/(epsln + drho_dpress_test)
      write (stdoutunit,' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') &
      'squared sound speed (',s_test,',',t_test,',',p_test,') = ',speed2_test,' (m/s)^2'

      diff_test = beta_test-jmfwg_beta
      write (stdoutunit,' (a,e22.16,a)') 'diff from JMFWG = ',diff_test,' 1/psu'


      write (stdoutunit,'(/,a)') 'NEUTRAL DENSITY EQUATION TEST VALUES'

      write (stdoutunit,'(a,f6.2,a,f6.2)') &
      'sn_test(psu) = ',sn_test,', tn_test(C) = ',tn_test

      neutralrho_test = neutral_density(sn_test,tn_test) 
      write (stdoutunit,' (a,f6.2,a,f6.2,a,e22.16,a)') &
      'rho  (',sn_test,',',tn_test,') = ',neutralrho_test,' kg/m^3'

      diff_test = neutralrho_test-mj_neutralrho
      write (stdoutunit,' (a,e22.16,a)') 'diff from McDougall and Jackett test = ',diff_test,' kg/m^3'

  endif 

  ! define vertical axes according to potential density classes  
  allocate ( Dens%potrho_ref(layer_nk))
  allocate ( potrho_bounds(layer_nk+1))
  potrho_bounds(1) = potrho_min
  potrho_interval  = (potrho_max-potrho_min)/(epsln+layer_nk)
  do k=2,layer_nk+1
    potrho_bounds(k)=potrho_bounds(k-1)+potrho_interval
  enddo 
  do k=1,layer_nk
    Dens%potrho_ref(k)=potrho_bounds(k)+0.5*potrho_interval
  enddo 

  ! define vertical axes according to neutral density classes  
  allocate ( Dens%neutralrho_ref(layer_nk))
  allocate ( neutralrho_bounds(layer_nk+1))
  neutralrho_bounds(1) = neutralrho_min
  neutralrho_interval  = (neutralrho_max-neutralrho_min)/(epsln+layer_nk)
  do k=2,layer_nk+1
    neutralrho_bounds(k)=neutralrho_bounds(k-1)+neutralrho_interval
  enddo 
  do k=1,layer_nk
    Dens%neutralrho_ref(k)=neutralrho_bounds(k)+0.5*neutralrho_interval
  enddo 

  ! define vertical axes according to 
  ! potential temperature or conservative temperature classes
  allocate ( Dens%theta_ref(layer_nk))
  allocate ( theta_bounds(layer_nk+1))
  theta_bounds(1)  = theta_min
  theta_interval   = (theta_max-theta_min)/(epsln+layer_nk)
  do k=2,layer_nk+1
    theta_bounds(k) =theta_bounds(k-1) +theta_interval
  enddo 
  do k=1,layer_nk
    Dens%theta_ref(k) =theta_bounds(k) +0.5*theta_interval
  enddo 

! assuming that if the Grd%name = 'ocean' then we will
! be using this diagnostic axis. Otherwise it is not used.
! This is done to prevent duplicate axis names for multiple
! grid objects.
  
  if (trim(Grd%name) == 'ocean') then
    id_potrho_xt = diag_axis_init ('rho_xt_'//trim(Grd%name),Grd%grid_x_t,'degrees_E','x', &
                                   'tcell longitude',set_name='ocean', Domain2=Dom%domain2d)

    id_potrho_yt = diag_axis_init ('rho_yt_'//trim(Grd%name),Grd%grid_y_t,'degrees_N','y', &
                                   'tcell latitude',set_name='ocean', Domain2=Dom%domain2d)

    id_potrho_bounds = diag_axis_init &
     ( 'potrho_edges',potrho_bounds, 'kg/m^3', 'z','potential density edges',&
       direction=-1, set_name='ocean')

    id_potrho_axis   = diag_axis_init &
     ( 'potrho',Dens%potrho_ref,&
       'kg/m^3', 'z','potential density',edges=id_potrho_bounds,&
        direction=-1,set_name='ocean')

    Dens%potrho_axes = (/ id_potrho_xt, id_potrho_yt, id_potrho_axis /)


    id_neutralrho_xt = diag_axis_init ('neutralrho_xt_'//trim(Grd%name),Grd%grid_x_t,'degrees_E','x', &
                                       'tcell longitude',set_name='ocean', Domain2=Dom%domain2d)

    id_neutralrho_yt = diag_axis_init ('neutralrho_yt_'//trim(Grd%name),Grd%grid_y_t,'degrees_N','y', &
                                       'tcell latitude',set_name='ocean', Domain2=Dom%domain2d)

    id_neutralrho_bounds = diag_axis_init &
     ( 'neutralrho_edges',neutralrho_bounds, 'kg/m^3', 'z','neutral density edges',&
       direction=-1, set_name='ocean')

    id_neutralrho_axis   = diag_axis_init &
     ( 'neutral',Dens%neutralrho_ref,&
       'kg/m^3', 'z','neutral density',edges=id_neutralrho_bounds,&
        direction=-1,set_name='ocean')

    Dens%neutralrho_axes = (/ id_neutralrho_xt, id_neutralrho_yt, id_neutralrho_axis /)


    id_theta_xt  = diag_axis_init ('theta_xt_'//trim(Grd%name),Grd%grid_x_t,'degrees_E','x', &
                                   'tcell longitude',set_name='ocean', Domain2=Dom%domain2d)

    id_theta_yt  = diag_axis_init ('theta_yt_'//trim(Grd%name),Grd%grid_y_t,'degrees_N','y', &
                                   'tcell latitude', set_name='ocean', Domain2=Dom%domain2d)

    id_theta_bounds = diag_axis_init &
     ( 'theta_edges',potrho_bounds, 'C', 'z','potential or conservative temperature edges',&
       direction=1, set_name='ocean')

    id_theta_axis   = diag_axis_init &
     ( 'theta',Dens%theta_ref, 'C', 'z','potential or conservative temperature',edges=id_theta_bounds, &
       direction=1,set_name='ocean')

    Dens%theta_axes = (/ id_theta_xt,  id_theta_yt,  id_theta_axis  /)

  endif


! register diagnostic fields for diag_manager 

  id_press = register_diag_field ('ocean_model', 'press', Grd%tracer_axes(1:3),&
             Time%model_time, 'absolute pressure', 'dbar',                     &
             missing_value=missing_value, range=(/-10.,6000./))
  id_rho = register_diag_field ('ocean_model', 'rho', Grd%tracer_axes(1:3), &
           Time%model_time, 'in situ density', 'kg/m^3',                    &
           missing_value=missing_value, range=(/-10.0,1e5/))
  id_pot_rho = register_diag_field ('ocean_model', 'pot_rho', Grd%tracer_axes(1:3), &
               Time%model_time, 'potential density', 'kg/m^3',                      &
               missing_value=missing_value, range=(/-10.0,1e5/))    
  id_neutral_rho = register_diag_field ('ocean_model', 'neutral_rho', Grd%tracer_axes(1:3), &
                   Time%model_time, 'polynomial estimate of neutral density', 'kg/m^3',     &
                   missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_0 = register_diag_field ('ocean_model', 'pot_rho_0', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 0 dbar', 'kg/m^3',   &
                 missing_value=missing_value, range=(/-10.0,1e5/),                      &
                 standard_name='sea_water_potential_density')    
  id_pot_rho_1 = register_diag_field ('ocean_model', 'pot_rho_1', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 1000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_2 = register_diag_field ('ocean_model', 'pot_rho_2', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 2000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_3 = register_diag_field ('ocean_model', 'pot_rho_3', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 3000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_4 = register_diag_field ('ocean_model', 'pot_rho_4', Grd%tracer_axes(1:3), &
                 Time%model_time, 'potential density referenced to 4000 dbar', 'kg/m^3',&
                 missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_et = register_diag_field ('ocean_model', 'pot_rho_et', Grd%tracer_axes_flux_x(1:3), &
                  Time%model_time, 'potential density at east face of tracer cell', 'kg/m^3',    &
                  missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_nt = register_diag_field ('ocean_model', 'pot_rho_nt', Grd%tracer_axes_flux_y(1:3), &
                  Time%model_time, 'potential density at north face of tracer cell', 'kg/m^3',   &
                  missing_value=missing_value, range=(/-10.0,1e5/))    
  id_pot_rho_wt = register_diag_field ('ocean_model', 'pot_rho_wt', Grd%tracer_axes_wt(1:3), &
                  Time%model_time, 'potential density at wt points', 'kg/m^3',               &
                  missing_value=missing_value, range=(/-10.0,1e5/))    
  id_drhodtheta = register_diag_field ('ocean_model', 'drhodtheta', Grd%tracer_axes(1:3), &
                  Time%model_time, 'd(rho)/d(theta)', 'kg/m^3/C',                         &
                  missing_value=-1e10, range=(/-1e10,1e10/))
  id_drhodsalt   = register_diag_field ('ocean_model', 'drhodsalinity', Grd%tracer_axes(1:3), &
                   Time%model_time, 'd(rho)/d(salinity)', 'kg/m^3/psu',                       &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_drhodpress  = register_diag_field ('ocean_model', 'drhodpress', Grd%tracer_axes(1:3), &
                   Time%model_time, 'd(rho)/d(press)', 'kg/m^3/Pa',                        &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_sound_speed2 = register_diag_field ('ocean_model', 'sound_speed2', Grd%tracer_axes(1:3),&
                   Time%model_time, 'squared sound speed = 1/[d(rho)/d(press)]', '(m/s)^2',  &
                   missing_value=missing_value, range=(/-1e9,1e9/))
  id_cabbeling   = register_diag_field ('ocean_model', 'cabbeling', Grd%tracer_axes(1:3), &
                  Time%model_time, 'cabbeling parameter', '(1/degC)^2',                   &
                  missing_value=-1e10, range=(/-1e10,1e10/))
  id_thermobaricity = register_diag_field ('ocean_model', 'thermobaricity', Grd%tracer_axes(1:3), &
                  Time%model_time, 'thermobaricity parameter', '1/(dbar*degC)',                   &
                  missing_value=-1e10, range=(/-1e10,1e10/))

  id_buoyfreq2_zt = register_diag_field ('ocean_model','buoyfreq2_zt',  &
                    Grd%tracer_axes(1:3), Time%model_time,              &
                    'Squared buoyancy frequency at T-point', '1/s^2',   &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_buoyfreq2_wt = register_diag_field ('ocean_model','buoyfreq2_wt',      &
                    Grd%tracer_axes_wt(1:3), Time%model_time,               &
                    'Squared buoyancy frequency at T-cell bottom', '1/s^2', &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_drhodz_zt = register_diag_field ('ocean_model','drhodz_zt',  &
                    Grd%tracer_axes(1:3), Time%model_time,        &
                    'd(neutral rho)/dz at T-point', 'kg/m^4',     &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_drhodz_wt = register_diag_field ('ocean_model','drhodz_wt',  &
                    Grd%tracer_axes_wt(1:3), Time%model_time,     &
                    'd(neutral rho)/dz at W-point', 'kg/m^4',     &
                    missing_value=missing_value, range=(/-1e6,1e6/))
  id_int_rhodz = register_diag_field ('ocean_model','int_rhodz',                        &
                    Grd%tracer_axes(1:2), Time%model_time,                              &
                    'vertical sum of in-situ density times cell thickness dz', 'kg/m^2',&
                    missing_value=missing_value, range=(/-1e1,1e20/))
  id_rhoave    = register_diag_field('ocean_model','rhoave',            &
                   Time%model_time, 'global mean ocean in-situ density from ocean_density_mod',&
                   'kg/m^3', missing_value=missing_value, range=(/-1e10,1e10/))
  id_pbot_adjust = register_diag_field('ocean_model','pbot_adjust',                                          &
                   Time%model_time, 'pbot adjustment to counteract spurious mass source in Boussinesq fluid',&
                   'dbar', missing_value=missing_value, range=(/-1e10,1e10/))
  id_eta_adjust = register_diag_field('ocean_model','eta_adjust',                                        &
                   Time%model_time, 'global eta adjustment to include steric effect in Boussinesq fluid',&
                   'm', missing_value=missing_value, range=(/-100.0,100.0/))
  id_eta_adjust_approx = register_diag_field('ocean_model','eta_adjust_approx',                                      &
                   Time%model_time, 'approximate global eta adjustment to include steric effect in Boussinesq fluid',&
                   'm', missing_value=missing_value, range=(/-100.0,100.0/))
  id_eta_adjust_cstvolume = register_diag_field('ocean_model','eta_adjust_cstvolume',                                               &
                   Time%model_time, 'global eta adjustment (assuming constant volume) to include steric effect in Boussinesq fluid',&
                   'm', missing_value=missing_value, range=(/-100.0,100.0/))


  end subroutine ocean_density_init
! </SUBROUTINE> NAME="ocean_density_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_diag">
!
! <DESCRIPTION>
! Diagnose pressure_at_depth and diagnostic ocean density fields.  
! Also send some diagnostics to diagnostic manager.  
! </DESCRIPTION>
!
  subroutine ocean_density_diag(Time, Temp, Salt, Dens, T_diag)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_prog_tracer_type),   intent(in)    :: Temp
  type(ocean_prog_tracer_type),   intent(in)    :: Salt
  type(ocean_density_type),       intent(inout) :: Dens
  type(ocean_diag_tracer_type),   intent(inout) :: T_diag(:)
  
  logical :: used
  integer :: tau

  tau   = Time%tau

  if (ind_rho .gt. 0) then
    T_diag(ind_rho)%field(:,:,:) = Dens%rho(:,:,:,tau)
  endif

  ! compute neutral density for diagnostic purposes 
  Dens%neutralrho(:,:,:) = Grd%tmask(:,:,:) &
   *neutral_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau))
  if (ind_neutralrho .gt. 0) then
    T_diag(ind_neutralrho)%field(:,:,:) = Dens%neutralrho(:,:,:)
  endif

  ! compute potential density for diagnostic purposes 
  Dens%potrho(:,:,:) = Grd%tmask(:,:,:) &
   *potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), potrho_press)


  ! compute potential densities for diagnostic purposes 
  if (ind_potrho_0 .gt. 0) then
    T_diag(ind_potrho_0)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 0.0)
  endif
  if (ind_potrho_1 .gt. 0) then
    T_diag(ind_potrho_1)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 1000.0)
  endif
  if (ind_potrho_2 .gt. 0) then
    T_diag(ind_potrho_2)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 2000.0)
  endif
  if (ind_potrho_3 .gt. 0) then
    T_diag(ind_potrho_3)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 3000.0)
  endif
  if (ind_potrho_4 .gt. 0) then
    T_diag(ind_potrho_4)%field(:,:,:) = Grd%tmask(:,:,:) *                &
         potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 4000.0)
  endif

  if (id_neutral_rho > 0) & 
       used = send_data (id_neutral_rho, Dens%neutralrho(:,:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,:),                  &
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_pot_rho > 0) & 
       used = send_data (id_pot_rho, Dens%potrho(:,:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,:), &
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  if (id_pot_rho_0 > 0) then 
    if (ind_potrho_0 .gt. 0) then
      used = send_data (id_pot_rho_0, T_diag(ind_potrho_0)%field(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                      &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    else
      wrk1(:,:,:) = potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 0.0)
      used = send_data (id_pot_rho_0, wrk1(:,:,:),  &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
  endif

  if (id_pot_rho_1 > 0) then 
    if (ind_potrho_1 .gt. 0) then
      used = send_data (id_pot_rho_1, T_diag(ind_potrho_1)%field(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                      &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    else
      wrk1(:,:,:) = potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 1000.0)
      used = send_data (id_pot_rho_1, wrk1(:,:,:),  &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
  endif

  if (id_pot_rho_2 > 0) then 
    if (ind_potrho_2 .gt. 0) then
      used = send_data (id_pot_rho_2, T_diag(ind_potrho_2)%field(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                      &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    else
      wrk1(:,:,:) = potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 2000.0)
      used = send_data (id_pot_rho_2, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),&
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
  endif

  if (id_pot_rho_3 > 0) then 
    if (ind_potrho_3 .gt. 0) then
      used = send_data (id_pot_rho_3, T_diag(ind_potrho_3)%field(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                      &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    else
      wrk1(:,:,:) = potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 3000.0)
      used = send_data (id_pot_rho_3, wrk1(:,:,:),  &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
  endif

  if (id_pot_rho_4 > 0) then 
    if (ind_potrho_4 .gt. 0) then
      used = send_data (id_pot_rho_4, T_diag(ind_potrho_4)%field(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),                      &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    else
      wrk1(:,:,:) = potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 4000.0)
      used = send_data (id_pot_rho_4, wrk1(:,:,:),  &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
  endif


end subroutine ocean_density_diag
! </SUBROUTINE> NAME="ocean_density_diag"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_density">
!
! <DESCRIPTION>
! Diagnose pressure_at_depth and ocean density.  
! Also send some diagnostics to diagnostic manager.  
! </DESCRIPTION>
!
  subroutine update_ocean_density(Time, Thickness, Temp, Salt, Ext_mode, Dens)

  type(ocean_time_type),          intent(in)    :: Time
  type(ocean_thickness_type),     intent(in)    :: Thickness
  type(ocean_prog_tracer_type),   intent(in)    :: Temp
  type(ocean_prog_tracer_type),   intent(in)    :: Salt
  type(ocean_external_mode_type), intent(in)    :: Ext_mode
  type(ocean_density_type),       intent(inout) :: Dens
  type(time_type)                               :: next_time 
  
  logical :: used
  integer :: i, j, k
  integer :: tau, taup1
  real    :: mass, massip1, massjp1
  real    :: density_tau         =0.0
  real    :: density_taup1       =0.0
  real    :: volume_tau          =0.0
  real    :: volume_taup1        =0.0
  real    :: volume_taup12       =0.0
  real    :: mass_tau            =0.0
  real    :: mass_taup1          =0.0
  real    :: dArea               =0.0 
  real    :: pbot_adjust         =0.0
  real    :: rhoave              =0.0
  real    :: eta_adjust          =0.0
  real    :: eta_adjust_approx   =0.0
  real    :: eta_adjust_cstvolume=0.0

  integer :: stdoutunit 
  stdoutunit=stdout() 

  tau   = Time%tau
  taup1 = Time%taup1  

  ! compute pressure at depth (dbars)
  if(vert_coordinate==GEOPOTENTIAL) then 
      Dens%pressure_at_depth(:,:,:) = pressure_in_dbars(Thickness, Dens%rho(:,:,:,tau), &
                                      Ext_mode%patm_t(:,:,tau)+grav*Dens%rho(:,:,1,tau)*Ext_mode%eta_t(:,:,tau))
  elseif(vert_coordinate==ZSTAR .or. vert_coordinate==ZSIGMA) then 
      Dens%pressure_at_depth(:,:,:) = pressure_in_dbars(Thickness, Dens%rho(:,:,:,tau), &
                                      Ext_mode%patm_t(:,:,tau))
  elseif(vert_coordinate==PRESSURE) then 
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%pressure_at_depth(i,j,k) = c2dbars*Thickness%depth_st(i,j,k)
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PSTAR) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%pressure_at_depth(i,j,k) = c2dbars*(                   &
                      Ext_mode%patm_t(i,j,tau)                             &
                    + Thickness%depth_st(i,j,k)*Thickness%pbot0r(i,j)      &
                      *(Ext_mode%pbot_t(i,j,tau)-Ext_mode%patm_t(i,j,tau)) &
                    )
            enddo
         enddo
      enddo
  elseif(vert_coordinate==PSIGMA) then
      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Dens%pressure_at_depth(i,j,k) = c2dbars*(                   &
                      Ext_mode%patm_t(i,j,tau)                             &
                    + Thickness%depth_st(i,j,k)                            &
                      *(Ext_mode%pbot_t(i,j,tau)-Ext_mode%patm_t(i,j,tau)) &
                    )
            enddo
         enddo
      enddo
  endif
  
  ! diagnose in situ density 
  Dens%rho(:,:,:,taup1) = Grd%tmask(:,:,:) &
   *density(Salt%field(:,:,:,taup1), Temp%field(:,:,:,taup1), Dens%pressure_at_depth(:,:,:))

  if(rho0_density) then 
     Dens%rho(:,:,:,taup1) = Grd%tmask(:,:,:)*rho0
  endif 

  ! compute drhodT and drhodS and drhodP
  call density_derivs(Dens%rho(:,:,:,tau), Salt%field(:,:,:,tau),           &
                      Temp%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:), &
                      Time, Dens%drhodT(:,:,:), Dens%drhodS(:,:,:), Dens%drhodP(:,:,:)) 


  ! compute cabbeling and thermobaricity parameters for diagnostics 
  if(id_cabbeling > 0 .or. id_thermobaricity > 0) then 
      call cabbeling_thermobaricity(Dens%rho(:,:,:,tau), Salt%field(:,:,:,tau),           &
                                    Temp%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:), &
                                    Time, Dens%drhodT(:,:,:), Dens%drhodS(:,:,:), Dens%drhodP(:,:,:))
  endif 

  ! compute buoyancy frequency for diagnostic purposes 
  call compute_buoyfreq(Time, Thickness, Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), Dens) 

  if (id_press > 0) &
       used = send_data (id_press, Dens%pressure_at_depth(:,:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,:), &
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_rho > 0) &
       used = send_data (id_rho, Dens%rho(:,:,:,tau), &
       Time%model_time, rmask=Grd%tmask(:,:,:), &
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if(id_int_rhodz > 0) then 
      wrk1_2d(:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_2d(i,j) = wrk1_2d(i,j) + Thickness%dzt(i,j,k)*Dens%rho(i,j,k,tau) 
            enddo
         enddo
      enddo
      used = send_data (id_int_rhodz, wrk1_2d(:,:), &
           Time%model_time,rmask=Grd%tmask(:,:,1),  &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  next_time = increment_time(Time%model_time, int(dtts), 0)

  if (need_data(id_pot_rho_wt,next_time)) then 
      wrk1(:,:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,nk) = Dens%potrho(i,j,nk)
         enddo
      enddo
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               if(Grd%tmask(i,j,k) > 0.0) then   
                   wrk1(i,j,k) = ( Grd%tmask(i,j,k+1)*Dens%potrho(i,j,k+1)*Thickness%rho_dzt(i,j,k+1,tau) &
                        +Dens%potrho(i,j,k)  *Thickness%rho_dzt(i,j,k,tau))  &
                        /( epsln + Thickness%rho_dzt(i,j,k,tau) + Grd%tmask(i,j,k+1)*Thickness%rho_dzt(i,j,k+1,tau) ) 
               endif
            enddo
         enddo
      enddo
      used = send_data (id_pot_rho_wt, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  if (need_data(id_pot_rho_et,next_time)) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               mass      = Thickness%rho_dzt(i,j,k,tau)  *Grd%dat(i,j)  *Grd%tmask(i,j,k) 
               massip1   = Thickness%rho_dzt(i+1,j,k,tau)*Grd%dat(i+1,j)*Grd%tmask(i+1,j,k) 
               wrk1(i,j,k) = ( Dens%potrho(i,j,k)*mass + Dens%potrho(i+1,j,k)*massip1 ) &
                    /( epsln + mass + massip1 ) 
            enddo
         enddo
      enddo
      used = send_data (id_pot_rho_et, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  if (need_data(id_pot_rho_nt,next_time)) then 
      wrk1(:,:,:) = 0.0
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               mass      = Thickness%rho_dzt(i,j,k,tau)  *Grd%dat(i,j)  *Grd%tmask(i,j,k) 
               massjp1   = Thickness%rho_dzt(i,j+1,k,tau)*Grd%dat(i,j+1)*Grd%tmask(i,j+1,k) 
               wrk1(i,j,k) = ( Dens%potrho(i,j,k)*mass + Dens%potrho(i,j+1,k)*massjp1 ) &
                    /( epsln + mass + massjp1 ) 
            enddo
         enddo
      enddo
      used = send_data (id_pot_rho_nt, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


  if (id_pbot_adjust>0       .or. id_eta_adjust>0           .or. &
      id_eta_adjust_approx>0 .or. id_eta_adjust_cstvolume>0 .or. id_rhoave>0) then 
    ! These diagnostics are used to adjust the bottom pressure and surface height computed
    ! from a Boussinesq model.  The adjustments correct, globally, for the spurious 
    ! mass sources appearing in the Boussinesq fluid and the missing steric effect.  
    ! The adjustments are NOT needed with non-Boussinesq models (pressure based vert coord),
    ! since the non-Boussinesq fluid correctly computes the bottom pressure and 
    ! surface height according to mass conserving kinematics.
    !
    ! In a Bouss model (depth based vert coordinates), when global mean density increases, 
    ! so does the Boussinesq mass. So the pbot_adjust is negative in this case, in 
    ! order to globally counteract effects from the spurious mass source.
    !
    ! In a non-Bouss model, when global mean density increases, global mean sea level decreases. 
    ! This global steric effect is absent from the Boussinesq fluid.  To incorporate this 
    ! this missing contribution to sea level, eta_adjust will be negative in this case. 
    !
    ! note that dzt has already been updated to its taup1 value. that is why we use  
    ! rho0r*rho_dzt = dzt (when depth-based vertical coordinates).

      mass_tau         = 0.0
      mass_taup1       = 0.0
      volume_tau       = 0.0
      volume_taup1     = 0.0
      volume_taup12    = 0.0
      rhodz_tau(:,:)   = 0.0
      rhodz_taup1(:,:) = 0.0
      wrk1_2d(:,:)     = 0.0
      wrk2_2d(:,:)     = 0.0
      wrk3_2d(:,:)     = 0.0
      wrk4_2d(:,:)     = 0.0

      ! do the vertical sum of rho*dzt 
      if(vert_coordinate==GEOPOTENTIAL .or. vert_coordinate==ZSTAR .or. vert_coordinate==ZSIGMA) then 
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   rhodz_tau(i,j)  = rhodz_tau(i,j)  + Grd%tmask(i,j,k)*rho0r*Thickness%rho_dzt(i,j,k,tau)   &
                                                       *Dens%rho(i,j,k,tau)
                   rhodz_taup1(i,j)= rhodz_taup1(i,j)+ Grd%tmask(i,j,k)*rho0r*Thickness%rho_dzt(i,j,k,taup1) &
                                                       *Dens%rho(i,j,k,taup1)
                enddo
             enddo
          enddo
      else 
          k=1
          do j=jsc,jec
             do i=isc,iec
                rhodz_tau(i,j)   = Grd%tmask(i,j,k)*grav_r*(Ext_mode%pbot_t(i,j,tau)  -Ext_mode%patm_t(i,j,tau))
                rhodz_taup1(i,j) = Grd%tmask(i,j,k)*grav_r*(Ext_mode%pbot_t(i,j,taup1)-Ext_mode%patm_t(i,j,taup1))
             enddo
          enddo
      endif

      ! compute the mass and volume of a seawater column
      do j=jsd,jed
         do i=isd,ied
            wrk1_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*rhodz_tau(i,j)
            wrk2_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*rhodz_taup1(i,j)
            wrk3_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*(Grd%ht(i,j) + Ext_mode%eta_t(i,j,tau))
            wrk4_2d(i,j) = Grd%tmask(i,j,1)*Grd%dat(i,j)*(Grd%ht(i,j) + Ext_mode%eta_t(i,j,taup1))
         enddo
      enddo

      ! perform global sums to get global mass and volume for ocean 
      mass_tau    = mpp_global_sum(Dom%domain2d,wrk1_2d(:,:),NON_BITWISE_EXACT_SUM)
      volume_tau  = mpp_global_sum(Dom%domain2d,wrk3_2d(:,:),NON_BITWISE_EXACT_SUM)
      density_tau = mass_tau/(epsln+volume_tau)

      mass_taup1    = mpp_global_sum(Dom%domain2d,wrk2_2d(:,:),NON_BITWISE_EXACT_SUM)
      volume_taup1  = mpp_global_sum(Dom%domain2d,wrk4_2d(:,:),NON_BITWISE_EXACT_SUM)
      density_taup1 = mass_taup1/(epsln+volume_taup1)

      volume_taup12 = 0.5*(volume_taup1+volume_tau)

      if(id_pbot_adjust > 0) then 
         pbot_adjust = -grav*c2dbars*area_total_r*volume_taup1*(density_taup1-density_tau) 
         used = send_data(id_pbot_adjust, pbot_adjust, Time%model_time)
      endif 
      if(id_eta_adjust > 0) then 
         eta_adjust = area_total_r*volume_taup12*log((density_tau+epsln)/(density_taup1+epsln))
         used = send_data(id_eta_adjust, eta_adjust, Time%model_time)
      endif 
      if(id_eta_adjust_approx > 0) then 
         eta_adjust_approx = area_total_r*volume_taup1*(density_tau/(epsln+density_taup1)-1.0)
         used = send_data(id_eta_adjust_approx, eta_adjust_approx, Time%model_time)
      endif 
      if(id_eta_adjust_cstvolume > 0) then 
         eta_adjust_cstvolume = area_total_r*Grd%tcellv*(density_tau/(epsln+density_taup1)-1.0)
         used = send_data(id_eta_adjust_cstvolume, eta_adjust_cstvolume, Time%model_time)
      endif 
      if(id_rhoave > 0) then 
         used = send_data(id_rhoave, density_tau, Time%model_time)
      endif    

  endif 

  if(debug_this_module) then 
      write(stdoutunit,'(a)') ' ' 
      write(stdoutunit,*) 'From ocean_density_mod: density chksums'
      call write_timestamp(Time%model_time)
      call ocean_density_chksum(Time, Dens)
  endif


end subroutine update_ocean_density
! </SUBROUTINE> NAME="update_ocean_density"


!#######################################################################
! <FUNCTION NAME="density_field">
!
! <DESCRIPTION>
! Compute density for all grid points.  
!
! Note that pressure here is 
!
! sea pressure = absolute pressure - press_standard (dbars) 
!
! and salinity is in model units (psu).  
!
! </DESCRIPTION>
!
  function density_field (salinity, theta, press)

    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta
    real, dimension(isd:,jsd:,:), intent(in) :: press

    real, dimension(isd:ied,jsd:jed,nk) :: density_field

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_field): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_field(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
             enddo
          enddo
       enddo

    else 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,k)
                t2  = t1*t1

                s1  = salinity(i,j,k)
                sp5 = sqrt(s1) 

                p1   = press(i,j,k) - press_standard
                p1t1 = p1*t1

                num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                     + s1*(a4 + a5*t1  + a6*s1)        &
                     + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                     + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                     + p1*(b10 + p1t1*(b11*t2 + b12*p1))
 
                denominator_r(i,j,k) = 1.0/(epsln+den) 

                density_field(i,j,k) = num*denominator_r(i,j,k)

             enddo
          enddo
       enddo
    endif

  end function density_field
! </FUNCTION> NAME="density_field"


!#######################################################################
! <FUNCTION NAME="density_level">
!
! <DESCRIPTION>
! Compute density at a particular k-level. 
!
! Note that pressure here is 
!
! sea pressure = absolute pressure - press_standard (dbars)
!
! </DESCRIPTION>
!
  function density_level(salinity, theta, press)

    real, dimension(isd:,jsd:), intent(in) :: salinity
    real, dimension(isd:,jsd:), intent(in) :: theta
    real, dimension(isd:,jsd:), intent(in) :: press

    real, dimension(isd:ied,jsd:jed) :: density_level 

    integer :: i, j
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_level): module must be initialized')
    endif 

    if(linear_eos) then

       do j=jsd,jed
          do i=isd,ied
             density_level(i,j) = rho0 - alpha_linear_eos*theta(i,j) + beta_linear_eos*salinity(i,j)
          enddo
       enddo

    else 

       do j=jsd,jed
          do i=isd,ied

             t1  = theta(i,j)
             t2  = t1*t1

             s1  = salinity(i,j)
             sp5 = sqrt(s1) 

             p1   = press(i,j) - press_standard
             p1t1 = p1*t1

             num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                  + s1*(a4 + a5*t1  + a6*s1)        &
                  + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

             den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                  + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                  + p1*(b10 + p1t1*(b11*t2 + b12*p1))

             density_level(i,j) = num/(epsln+den)

          enddo
       enddo
    endif

  end function density_level
! </FUNCTION> NAME="density_level"


!#######################################################################
! <FUNCTION NAME="density_line">
!
! <DESCRIPTION>
! Compute density at a particular k-level and j index.  This scheme
! is used in the vectorized version of the full convection scheme. 
!
! Note that pressure here is
!
! sea pressure = absolute pressure - press_standard
!
! </DESCRIPTION>
!
  function density_line(salinity, theta, press)

    real, dimension(isd:), intent(in) :: salinity
    real, dimension(isd:), intent(in) :: theta
    real, dimension(isd:), intent(in) :: press

    real, dimension(isd:ied) :: density_line 

    integer :: i
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_line): module must be initialized')
    endif 

    if(linear_eos) then

        do i=isd,ied
           density_line(i) = rho0 - alpha_linear_eos*theta(i) + beta_linear_eos*salinity(i)
        enddo

    else 

        do i=isd,ied

           t1  = theta(i)
           t2  = t1*t1

           s1  = salinity(i)
           sp5 = sqrt(s1) 

           p1   = press(i) - press_standard
           p1t1 = p1*t1

           num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                + s1*(a4 + a5*t1  + a6*s1)        &
                + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

           den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                  + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                  + p1*(b10 + p1t1*(b11*t2 + b12*p1))

           density_line(i) = num/(epsln+den)
       enddo
    endif

  end function density_line
! </FUNCTION> NAME="density_line"


!#######################################################################
! <FUNCTION NAME="neutral_density_field">
!
! <DESCRIPTION>
! Compute neutral density according to a rational polynomial 
! approximation given by McDougall and Jackett (2005).
! </DESCRIPTION>
!
  function neutral_density_field (salinity, theta)

    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta

    real, dimension(isd:ied,jsd:jed,nk) :: neutral_density_field

    integer :: i, j, k
    real    :: t1, t2, s1, sp5
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (neutral_density_field): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                neutral_density_field(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
             enddo
          enddo
       enddo

    else 

          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   t1  = theta(i,j,k)
                   t2  = t1*t1

                   s1  = salinity(i,j,k)
                   sp5 = sqrt(s1) 

                   num = a0n + t1*(a1n + t1*(a2n+a3n*t1) )    &
                        + s1*(a4n + a5n*t1  + a6n*s1)       

                   den = b0n + t1*(b1n + t1*(b2n + t1*(b3n + t1*b4n)))      &
                        + s1*(b5n + t1*(b6n + b7n*t2) + sp5*(b8n + b9n*t2)) 

                   neutral_density_field(i,j,k) = num/(epsln+den)
                enddo
             enddo
          enddo

    endif  ! endif for linear_eos 

  end function neutral_density_field
! </FUNCTION> NAME="neutral_density_field"


!#######################################################################
! <FUNCTION NAME="neutral_density_point">
!
! <DESCRIPTION>
! Compute neutral density according to a rational polynomial 
! approximation given by McDougall and Jackett (2005).
! </DESCRIPTION>
!
  function neutral_density_point (salinity, theta)

    real, intent(in) :: salinity
    real, intent(in) :: theta

    real :: neutral_density_point
    real :: t1, t2, s1, sp5
    real :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (neutral_density_point): module must be initialized')
    endif 

    if(linear_eos) then

        neutral_density_point = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity

    else 

        t1  = theta
        t2  = t1*t1

        s1  = salinity
        sp5 = sqrt(s1) 

        num = a0n + t1*(a1n + t1*(a2n+a3n*t1) )    &
             + s1*(a4n + a5n*t1  + a6n*s1)       

        den = b0n + t1*(b1n + t1*(b2n + t1*(b3n + t1*b4n)))      &
             + s1*(b5n + t1*(b6n + b7n*t2) + sp5*(b8n + b9n*t2)) 

        neutral_density_point = num/(epsln+den)

    endif  ! endif for linear_eos 

  end function neutral_density_point
! </FUNCTION> NAME="neutral_density_point"


!#######################################################################
! <FUNCTION NAME="potential_density">
!
! <DESCRIPTION>
! Compute potential density referenced to some given sea pressure. 
!
! Note that potential density referenced to the surface (i.e., sigma_0)
! has a zero sea pressure, so pressure=0.0 should be the argument
! to the function. 
!
! Note that pressure here is 
! sea pressure = absolute pressure - press_standard  (dbars)
!
! </DESCRIPTION>
!
  function potential_density (salinity, theta, pressure)

    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta
    real,                         intent(in) :: pressure

    real, dimension(isd:ied,jsd:jed,nk) :: potential_density

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (potential_density): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                potential_density(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
             enddo
          enddo
       enddo

    else 

       p1 = pressure

       if(p1 > 0.0) then 

          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   t1  = theta(i,j,k)
                   t2  = t1*t1

                   s1  = salinity(i,j,k)
                   sp5 = sqrt(s1) 

                   p1t1 = p1*t1

                   num = a0 + t1*(a1 + t1*(a2+a3*t1) )    &
                        + s1*(a4 + a5*t1  + a6*s1)        &
                        + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                   den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                        + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                        + p1*(b10 + p1t1*(b11*t2 + b12*p1))

                   potential_density(i,j,k) = num/(epsln+den)
                enddo
             enddo
          enddo

       elseif(p1==0.0) then ! for sigma_0

          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   t1  = theta(i,j,k)
                   t2  = t1*t1

                   s1  = salinity(i,j,k)
                   sp5 = sqrt(s1) 

                   num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                        + s1*(a4 + a5*t1  + a6*s1)        
                        
                   den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                        + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) 

                   potential_density(i,j,k) = num/(epsln+den)
                enddo
             enddo
          enddo

       elseif(p1 < 0.0) then 
          call mpp_error(FATAL, &
          '==>Error in ocean_density_mod: potential density at pressure < 0 is not defined')

       endif  ! endif for value of pressure 

    endif  ! endif for linear_eos 

  end function potential_density
! </FUNCTION> NAME="potential_density"


!#######################################################################
! <FUNCTION NAME="density_sfc">
!
! <DESCRIPTION>
! Compute density as a function of surface salinity, surface theta, 
! and insitu sea pressure. 
!
! Note that pressure here is 
! sea pressure = absolute pressure - press_standard  (dbars)
! 
! For use in KPP mixed layer scheme 
! </DESCRIPTION>
!
  function density_sfc (salinity, theta, press)

    real, intent(in), dimension(isd:,jsd:,:) :: salinity
    real, intent(in), dimension(isd:,jsd:,:) :: theta
    real, intent(in), dimension(isd:,jsd:,:) :: press

    real, dimension(isd:ied,jsd:jed,nk) :: density_sfc

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    if ( .not. module_is_initialized ) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_density_mod (density_sfc): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_sfc(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,1) + beta_linear_eos*salinity(i,j,1)
             enddo
          enddo
       enddo

    else 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,1)
                t2  = t1*t1

                s1  = salinity(i,j,1)
                sp5 = sqrt(s1) 

                p1   = press(i,j,k) - press_standard 
                p1t1 = p1*t1

                num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                     + s1*(a4 + a5*t1  + a6*s1)        &
                     + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                     + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                     + p1*(b10 + p1t1*(b11*t2 + b12*p1))

                density_sfc(i,j,k) = num/(epsln + den)

             enddo
          enddo
       enddo

    endif

  end function density_sfc
! </FUNCTION> NAME="density_sfc"


!#######################################################################
! <FUNCTION NAME="density_point">
!
! <DESCRIPTION>
! Compute density at a single model grid point. 
!
! Note that pressure here is 
!
! sea pressure = absolute pressure - press_standard  (dbars)
!
! </DESCRIPTION>
!
  function density_point (s1, t1, p1_dbars)

    real, intent(in) :: s1, t1, p1_dbars
    real             :: t2, sp5, p1, p1t1
    real             :: num, den
    real             :: density_point

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_point): module must be initialized')
    endif 

    if(linear_eos) then

       density_point = rho0 - alpha_linear_eos*t1 + beta_linear_eos*s1

    else 

       t2  = t1*t1
       sp5 = sqrt(s1) 

       p1   = p1_dbars - press_standard 
       p1t1 = p1*t1

       num = a0 + t1*(a1 + t1*(a2+a3*t1))                   &
            + s1*(a4 + a5*t1  + a6*s1)                      &
            + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

       den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4 )))      &
            + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2))  &  
            + p1*(b10 + p1t1*(b11*t2 + b12*p1))

       density_point = num/(epsln+den)

    endif

  end function density_point
! </FUNCTION> NAME="density_point"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_field">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to potential 
! temperature and with respect to salinity.  Hold pressure constant.  
!
! Pressure here is 
!
! sea pressure = absolute press - press_standard (dbars)
!
! </DESCRIPTION>
!
  subroutine density_derivs_field (rho, salinity, theta, press, Time, &
             density_theta, density_salinity, density_press)

    type(ocean_time_type),        intent(in)    :: Time
    real, dimension(isd:,jsd:,:), intent(in)    :: rho
    real, dimension(isd:,jsd:,:), intent(in)    :: salinity
    real, dimension(isd:,jsd:,:), intent(in)    :: theta
    real, dimension(isd:,jsd:,:), intent(in)    :: press
    real, dimension(isd:,jsd:,:), intent(inout) :: density_theta
    real, dimension(isd:,jsd:,:), intent(inout) :: density_salinity
    real, dimension(isd:,jsd:,:), intent(inout) :: density_press

    integer :: i, j, k
    real :: t1, t2, s1, sp5, p1, p1t1
    real :: dnum_dtheta,    dden_dtheta 
    real :: dnum_dsalinity, dden_dsalinity
    real :: dnum_dpress, dden_dpress

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_derivs_field): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_theta(i,j,k)    = -alpha_linear_eos
                density_salinity(i,j,k) =  beta_linear_eos
                density_press(i,j,k)    =  0.0
             enddo
          enddo
       enddo

    else 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,k)
                t2  = t1*t1
                s1  = salinity(i,j,k)
                sp5 = sqrt(s1) 

                p1   = press(i,j,k) - press_standard 
                p1t1 = p1*t1

                dnum_dtheta = a1 + t1*(two_a2 + three_a3*t1) &
                     + a5*s1                                 &
                     + p1t1*(two_a8 + two_a11*p1)    
                dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                     + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                     + p1*p1*(three_b11*t2 + b12*p1)

                dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
                dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

                dnum_dpress = a7 + a9*s1 + two_a10*p1 + t2*(a8 + two_a11*p1)
                dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 

                density_theta(i,j,k)    = denominator_r(i,j,k)*(dnum_dtheta    - rho(i,j,k)*dden_dtheta)
                density_salinity(i,j,k) = denominator_r(i,j,k)*(dnum_dsalinity - rho(i,j,k)*dden_dsalinity)
                density_press(i,j,k)    = denominator_r(i,j,k)*(dnum_dpress    - rho(i,j,k)*dden_dpress)*c2dbars

             enddo
          enddo
       enddo

    endif

    if (id_drhodtheta > 0) used = send_data (id_drhodtheta, density_theta(:,:,:), &
                                  Time%model_time, rmask=Grd%tmask(:,:,:),        &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    if (id_drhodsalt  > 0) used = send_data (id_drhodsalt, density_salinity(:,:,:), &
                                  Time%model_time, rmask=Grd%tmask(:,:,:),          &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    if (id_drhodpress > 0) used = send_data (id_drhodpress, density_press(:,:,:),&
                                  Time%model_time, rmask=Grd%tmask(:,:,:),       &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    if (id_sound_speed2 > 0) then 
        wrk1(:,:,:) = 0.0  
        do k=1,nk
           do j=jsd,jed
              do i=isd,ied
                 if(density_press(i,j,k) /= 0.0) then 
                     wrk1(i,j,k) = 1.0/density_press(i,j,k)
                 endif
              enddo
           enddo
        enddo
        used = send_data (id_sound_speed2, wrk1(:,:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,:),  &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif


  end subroutine density_derivs_field
! </SUBROUTINE> NAME="density_derivs_field"



!#######################################################################
! <SUBROUTINE NAME="cabbeling_thermobaricity">
!
! <DESCRIPTION>
! Compute cabbeling and thermobaricity parameters, as defined in 
! McDougall (1987).
!
! Pressure here is 
! sea pressure = absolute press - press_standard (dbars)
!
! </DESCRIPTION>
!
  subroutine cabbeling_thermobaricity(rho, salinity, theta, press, Time, &
             density_theta, density_salinity, density_press)

    real, dimension(isd:,jsd:,:), intent(in) :: rho
    real, dimension(isd:,jsd:,:), intent(in) :: salinity
    real, dimension(isd:,jsd:,:), intent(in) :: theta
    real, dimension(isd:,jsd:,:), intent(in) :: press
    type(ocean_time_type),        intent(in) :: Time
    real, dimension(isd:,jsd:,:), intent(in) :: density_theta
    real, dimension(isd:,jsd:,:), intent(in) :: density_salinity
    real, dimension(isd:,jsd:,:), intent(in) :: density_press

    integer :: i, j, k
    real :: t1, t2, s1, sp5, spm5, p1, p2, p1t1

    real :: rhotheta_rhosalinity
    real :: rho_inv
    real :: d2rho_dtheta2, d2rho_dsalin2, d2rho_dsalin_dtheta
    real :: d2rho_dsalin_dpress, d2rho_dtheta_dpress

    real :: dden_dtheta, dden_dsalinity, dden_dpress

    real :: d2num_dtheta2, d2num_dsalin2, d2num_dtheta_dsalin
    real :: d2num_dsalin_dpress, d2num_dtheta_dpress
    real :: d2den_dtheta2, d2den_dsalin2, d2den_dsalin_dtheta 
    real :: d2den_dsalin_dpress, d2den_dtheta_dpress

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (cabbeling_thermobaricity): module must be initialized')
    endif 

    wrk1(:,:,:) = 0.0   ! cabbeling parameter 
    wrk2(:,:,:) = 0.0   ! thermobaricity parameter 

    if(.not. linear_eos) then
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1   = theta(i,j,k)
                t2   = t1*t1
                s1   = salinity(i,j,k)
                sp5  = sqrt(s1) 
                p1   = press(i,j,k) - press_standard 
                p2   = p1*p1
                p1t1 = p1*t1

                if(s1 > 0.0) then 
                   spm5 = 1.0/sp5
                else 
                   spm5 = 0.0
                endif
                rhotheta_rhosalinity = Grd%tmask(i,j,k)*density_theta(i,j,k)/(epsln+density_salinity(i,j,k))

                dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                     + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                     + p1*p1*(three_b11*t2 + b12*p1)
                dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)
                dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 

                d2num_dtheta2 = two_a2 + six_a3*t1 + two_a8*p1 + two_a11*p2
                d2num_dsalin2 = two_a6
                d2num_dtheta_dsalin = a5
                d2num_dsalin_dpress = a9
                d2num_dtheta_dpress = two_a8*t1 + four_a11*p1*t1

                d2den_dtheta2 = two_b2 + six_b3*t1 + twelve_b4*t2  &
                              + six_b7*s1*t1 + two_b9*s1*sp5 + six_b11*p2*t1
                d2den_dsalin2 = .75*spm5*(b8 + b9*t2)
                d2den_dsalin_dtheta = b6 + three_b7*t2 + three_b9*sp5*t1
                d2den_dsalin_dpress = 0.0
                d2den_dtheta_dpress = six_b11*p1*t2 + three_b12*p2

                d2rho_dtheta2 = denominator_r(i,j,k) &
                 *(d2num_dtheta2-2.0*density_theta(i,j,k)*dden_dtheta-rho(i,j,k)*d2den_dtheta2)

                d2rho_dsalin2 = denominator_r(i,j,k) &
                 *(d2num_dsalin2-2.0*density_salinity(i,j,k)*dden_dsalinity-rho(i,j,k)*d2den_dsalin2)

                d2rho_dsalin_dtheta = denominator_r(i,j,k)  &
                 *( d2num_dtheta_dsalin                     &
                   -density_salinity(i,j,k)*dden_dtheta     &
                   -density_theta(i,j,k)*dden_dsalinity     &
                   -rho(i,j,k)*d2den_dsalin_dtheta)

                d2rho_dsalin_dpress = denominator_r(i,j,k)  &
                 *( d2den_dsalin_dpress                     &
                   -density_press(i,j,k)*dden_dsalinity     &
                   -density_salinity(i,j,k)*dden_dpress     &
                   -rho(i,j,k)*d2den_dsalin_dpress)

                d2rho_dtheta_dpress = denominator_r(i,j,k)  &
                 *( d2den_dtheta_dpress                     &
                   -density_press(i,j,k)*dden_dtheta        &
                   -density_theta(i,j,k)*dden_dpress        &
                   -rho(i,j,k)*d2den_dtheta_dpress)
 
                rho_inv = Grd%tmask(i,j,k)/(rho(i,j,k) + epsln)

                wrk1(i,j,k) = -rho_inv*                          &
                  ( d2rho_dtheta2                                &
                   -2.0*d2rho_dsalin_dtheta*rhotheta_rhosalinity &
                   +d2rho_dsalin2*rhotheta_rhosalinity**2 )
                   
                wrk2(i,j,k) = -rho_inv*(d2rho_dtheta_dpress-d2rho_dsalin_dpress*rhotheta_rhosalinity)

             enddo
          enddo
       enddo

    endif

    if(id_cabbeling > 0) then  
      used = send_data (id_cabbeling, wrk1(:,:,:),   &
             Time%model_time, rmask=Grd%tmask(:,:,:),&
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif 
    if(id_thermobaricity > 0) then 
      used = send_data (id_thermobaricity, wrk2(:,:,:),&
             Time%model_time, rmask=Grd%tmask(:,:,:),  &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif 



  end subroutine cabbeling_thermobaricity
! </SUBROUTINE> NAME="cabbeling_thermobaricity"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_level">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to potential 
! temperature and with respect to salinity.  Hold pressure constant.  
!
! Pressure here is sea pressure = absolute press - press_standard
!
! </DESCRIPTION>
!
  subroutine density_derivs_level (rho, salinity, theta, press, Time, klevel, &
                                   density_theta, density_salinity, density_press)

    type(ocean_time_type),      intent(in)    :: Time
    real, dimension(isd:,jsd:), intent(in)    :: rho
    real, dimension(isd:,jsd:), intent(in)    :: salinity
    real, dimension(isd:,jsd:), intent(in)    :: theta
    real, dimension(isd:,jsd:), intent(in)    :: press
    integer,                    intent(in)    :: klevel 
    real, dimension(isd:,jsd:), intent(inout) :: density_theta
    real, dimension(isd:,jsd:), intent(inout) :: density_salinity
    real, dimension(isd:,jsd:), intent(inout) :: density_press

    integer :: i, j
    real :: t1, t2, s1, sp5, p1, p1t1
    real :: dnum_dtheta,    dden_dtheta 
    real :: dnum_dsalinity, dden_dsalinity
    real :: dnum_dpress, dden_dpress

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_density_mod (density_derivs_level): module must be initialized')
    endif 

    if(linear_eos) then

          do j=jsd,jed
             do i=isd,ied
                density_theta(i,j)    = -alpha_linear_eos 
                density_salinity(i,j) =  beta_linear_eos
                density_press(i,j)    =  0.0
             enddo
          enddo

    else 

          do j=jsd,jed
             do i=isd,ied


                t1  = theta(i,j)
                t2  = t1*t1
                s1  = salinity(i,j)
                sp5 = sqrt(s1) 

                p1   = press(i,j) - press_standard 
                p1t1 = p1*t1

                dnum_dtheta = a1 + t1*(two_a2 + three_a3*t1) &
                     + a5*s1                                 &
                     + p1t1*(two_a8 + two_a11*p1)    
                dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                     + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                     + p1*p1*(three_b11*t2 + b12*p1)

                dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
                dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

                dnum_dpress = a7 + a9*s1 + two_a10*p1 + t2*(a8 + two_a11*p1)
                dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 

                density_theta(i,j)    = denominator_r(i,j,klevel)*(dnum_dtheta    - rho(i,j)*dden_dtheta)
                density_salinity(i,j) = denominator_r(i,j,klevel)*(dnum_dsalinity - rho(i,j)*dden_dsalinity)
                density_press(i,j)    = denominator_r(i,j,klevel)*(dnum_dpress    - rho(i,j)*dden_dpress)

             enddo
          enddo

    endif

  end subroutine density_derivs_level
! </SUBROUTINE> NAME="density_derivs_level"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_point">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to potential 
! temperature and with respect to salinity.  Do so here for a point. 
!
! Pressure here is 
!
! sea pressure = absolute pressure - press_standard  (dbars)
!
! </DESCRIPTION>
!
  subroutine density_derivs_point (rho, salinity, theta, press,&
                                   density_theta, density_salinity, density_press)

    real, intent(in)  :: rho, salinity, theta, press
    real, intent(out) :: density_theta, density_salinity,density_press
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: dnum_dtheta,    dden_dtheta 
    real    :: dnum_dsalinity, dden_dsalinity
    real    :: dnum_dpress, dden_dpress
    real    :: density_point, denrecip

    if ( .not. module_is_initialized ) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_density_mod (density_derivs_point): module must be initialized')
    endif 

    if(linear_eos) then

       density_theta    = -alpha_linear_eos 
       density_salinity =  beta_linear_eos
       density_press    =  0.0

    else 

       t1  = theta
       t2  = t1*t1
       s1  = salinity
       sp5 = sqrt(s1) 

       p1   = press - press_standard 
       p1t1 = p1*t1

       density_point = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4 ))) &
            + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2))       &  
            + p1*(b10 + p1t1*(b11*t2 + b12*p1))

       denrecip = 1.0/(density_point + epsln)

       dnum_dtheta = a1 + t1*(2.0*a2 + 3.0*a3*t1) &
            + a5*s1                               &
            + 2.0*p1t1*(a8 + a11*p1)    
       dden_dtheta = b1 + t1*(2.0*b2 + t1*(3.0*b3 + 4.0*b4*t1)) &
            + s1*(b6 + t1*(3.0*b7*t1 + 2.0*b9*sp5))             &
            + p1*p1*(3.0*b11*t2 + b12*p1)

       dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
       dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

       dnum_dpress = a7 + a9*s1 + two_a10*p1 + t2*(a8 + two_a11*p1)
       dden_dpress = b10 + p1*t1*(two_b11*t2 + three_b12*p1) 
 
       density_theta    = denrecip*(dnum_dtheta    - rho*dden_dtheta)
       density_salinity = denrecip*(dnum_dsalinity - rho*dden_dsalinity)
       density_press    = denrecip*(dnum_dpress    - rho*dden_dpress)

    endif

  end subroutine density_derivs_point
! </SUBROUTINE> NAME="density_derivs_point"


!#######################################################################
! <FUNCTION NAME="density_delta_z">
!
! <DESCRIPTION>
! rho(k)-rho(k+1) for all i,j with both temperatures referenced to the 
! deeper pressure depth.  
!
! Of use for KPP scheme. 
! </DESCRIPTION>
!
function density_delta_z (rho_initial, salinity, theta, press)

  real, intent(in), dimension(isd:,jsd:,:) :: rho_initial
  real, intent(in), dimension(isd:,jsd:,:) :: salinity
  real, intent(in), dimension(isd:,jsd:,:) :: theta
  real, intent(in), dimension(isd:,jsd:,:) :: press

  real, dimension(isd:ied,jsd:jed,nk)      :: density_delta_z
  integer :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_density_mod (density_delta_z): module must be initialized')
  endif 

  do k=1,nk-1
    wrk1(:,:,k) = press(:,:,k+1)
  enddo
  wrk1(:,:,nk) = press(:,:,nk)

  wrk2(:,:,:) = density(salinity(:,:,:), theta(:,:,:), wrk1(:,:,:))

  do k=1,nk-1
    density_delta_z(:,:,k) = (wrk2(:,:,k) - rho_initial(:,:,k+1))*Grd%tmask(:,:,k+1)
  enddo
  density_delta_z(:,:,nk) = (wrk2(:,:,nk) - rho_initial(:,:,nk))*Grd%tmask(:,:,nk)

end function density_delta_z
! </FUNCTION> NAME="density_delta_z"



!#######################################################################
! <FUNCTION NAME="density_delta_sfc">
!
! <DESCRIPTION>
! rho(1)-rho(k+1) for all i,j. 
!
! Of use for KPP scheme. 
! </DESCRIPTION>
!
function density_delta_sfc (rho_initial, salinity, theta, press)

  real, intent(in), dimension(isd:,jsd:,:) :: rho_initial
  real, intent(in), dimension(isd:,jsd:,:) :: salinity
  real, intent(in), dimension(isd:,jsd:,:) :: theta
  real, intent(in), dimension(isd:,jsd:,:) :: press
  real, dimension(isd:ied,jsd:jed,nk)      :: density_delta_sfc
  integer :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, &
    '==>Error in ocean_density_mod (density_delta_sfc): module must be initialized')
  endif 

  do k=1,nk-1
    wrk1(:,:,k) = press(:,:,k+1)
  enddo
  wrk1(:,:,nk) = press(:,:,nk)

  wrk2(:,:,:) = density_sfc(salinity, theta, wrk1)

  do k=1,nk-1
    density_delta_sfc(:,:,k) = (wrk2(:,:,k) - rho_initial(:,:,k+1))*Grd%tmask(:,:,k+1)
  enddo
  density_delta_sfc(:,:,nk) = (wrk2(:,:,nk) - rho_initial(:,:,nk))*Grd%tmask(:,:,nk)

end function density_delta_sfc
! </FUNCTION> NAME="density_delta_sfc"


!#######################################################################
! <SUBROUTINE NAME="compute_buoyfreq">
! <DESCRIPTION>
!
! Diagnose the buoyancy frequency, both at T-points and at 
! vertical interfaces of T-cells. The algorithm follows that
! used in subroutine diagnose_wdianeutral.  
!
! Author: Stephen.Griffies@noaa.gov
!
! </DESCRIPTION>
!
subroutine compute_buoyfreq(Time, Thickness, salinity, theta, Dens)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  real, dimension(isd:,jsd:,:), intent(in)    :: salinity
  real, dimension(isd:,jsd:,:), intent(in)    :: theta
  type(ocean_density_type),     intent(inout) :: Dens

  integer   :: i,j,k,kp1,kbot,m,n
  integer   :: tau
  real      :: drhodz_prev, tmpdrhodz, drhodzr

  if (.not. module_is_initialized) then 
      call mpp_error(FATAL, &
      '==>Error from ocean_tracer_diag (compute_buoyfreq): module needs initialization')
  endif

  tau = Time%tau
  wrk1(:,:,:) = 0.0
  wrk2(:,:,:) = 0.0
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0

  ! vertical derivative of temp and salt
  do k=1,nk-1
     kp1 = k+1
     do j=jsd,jed
        do i=isd,ied
           wrk1(i,j,k) = (theta(i,j,k)-theta(i,j,kp1))/Thickness%dzwt(i,j,k)
           wrk2(i,j,k) = (salinity(i,j,k)-salinity(i,j,kp1))/Thickness%dzwt(i,j,k)
        enddo
     enddo
  enddo

  ! vertical derivative of neutral density.
  ! vanishes at the bottom of a column.
  do k=1,nk-1
     kp1=k+1
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = Dens%drhodT(i,j,k)*wrk1(i,j,k) + Dens%drhodS(i,j,k)*wrk2(i,j,k)  
           if(abs(wrk3(i,j,k)) < epsln_drhodz) then 
              wrk3(i,j,k) = sign(1.0,wrk3(i,j,k))*epsln_drhodz
           endif 
           wrk3(i,j,k) = Grd%tmask(i,j,kp1)*wrk3(i,j,k)
        enddo
     enddo
  enddo

  ! vertically smooth drhodz; otherwise can get noisy results 
  if(buoyfreq_smooth_vert) then
      do m=1,num_121_passes
         do j=jsd,jed
            do i=isd,ied
               drhodz_prev = onefourth*wrk3(i,j,1) 
               kbot=Grd%kmt(i,j)
               if(kbot>3) then
                   do k=2,kbot-2
                      tmpdrhodz   = wrk3(i,j,k)
                      wrk3(i,j,k) = drhodz_prev + onehalf*wrk3(i,j,k) + onefourth*wrk3(i,j,k+1)
                      drhodz_prev = onefourth*tmpdrhodz
                   enddo
               endif
            enddo
         enddo
      enddo
  endif

  ! drhodz at bottom of T-cell
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Dens%drhodz_wt(i,j,k) = wrk3(i,j,k)
        enddo
     enddo
  enddo
  if(id_buoyfreq2_wt > 0) then 
      used = send_data (id_buoyfreq2_wt, -grav*rho0r*wrk3(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:),               & 
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if(id_drhodz_wt > 0) then 
      used = send_data (id_drhodz_wt, wrk3(:,:,:),  &
           Time%model_time, rmask=Grd%tmask(:,:,:), & 
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  ! drhodz and squared buoyancy frequency at T-cell point 
  wrk4(:,:,:) = 0.0
  k=1
  do j=jsd,jed
     do i=isd,ied
        Dens%drhodz_zt(i,j,k) = Dens%drhodz_wt(i,j,k)
        wrk4(i,j,k)           = Dens%drhodz_wt(i,j,k)
     enddo
  enddo
  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           wrk4(i,j,k) = wrk3(i,j,k-1)*(Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k))   &
                        +wrk3(i,j,k)  *(Thickness%depth_zt(i,j,k) -Thickness%depth_zwt(i,j,k-1)) 
           wrk4(i,j,k) = wrk4(i,j,k)/(epsln+Thickness%dzt(i,j,k))
           Dens%drhodz_zt(i,j,k) = wrk4(i,j,k)
        enddo
     enddo
  enddo
  if(id_buoyfreq2_zt > 0) then 
       used = send_data (id_buoyfreq2_zt, -grav*rho0r*wrk4(:,:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,:),                    & 
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if(id_drhodz_zt > 0) then 
       used = send_data (id_drhodz_zt, wrk4(:,:,:), &
       Time%model_time, rmask=Grd%tmask(:,:,:),     & 
       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 

  ! dTdz and dSdz  at T-cell point 
  wrk3(:,:,:) = 0.0
  wrk4(:,:,:) = 0.0
  k=1
  do j=jsd,jed
     do i=isd,ied
        Dens%dTdz_zt(i,j,k) = wrk1(i,j,k)
        Dens%dSdz_zt(i,j,k) = wrk2(i,j,k)
     enddo
  enddo

  do k=2,nk
     do j=jsd,jed
        do i=isd,ied
           wrk3(i,j,k) = wrk1(i,j,k-1)*(Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k))   &
                        +wrk1(i,j,k)  *(Thickness%depth_zt(i,j,k) -Thickness%depth_zwt(i,j,k-1)) 
           wrk4(i,j,k) = wrk2(i,j,k-1)*(Thickness%depth_zwt(i,j,k)-Thickness%depth_zt(i,j,k))   &
                        +wrk2(i,j,k)  *(Thickness%depth_zt(i,j,k) -Thickness%depth_zwt(i,j,k-1)) 
           wrk3(i,j,k) = wrk3(i,j,k)/(epsln+Thickness%dzt(i,j,k))
           wrk4(i,j,k) = wrk4(i,j,k)/(epsln+Thickness%dzt(i,j,k))
           Dens%dTdz_zt(i,j,k) = wrk3(i,j,k)
           Dens%dSdz_zt(i,j,k) = wrk4(i,j,k)
        enddo
     enddo
  enddo


end subroutine compute_buoyfreq
! </SUBROUTINE>  NAME="compute_buoyfreq"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_end">
!
! <DESCRIPTION>
!
! Write density and pressure_at_depth to a restart.
!
! </DESCRIPTION>
subroutine ocean_density_end(Time, Dens)

  type(ocean_time_type),    intent(in)           :: Time
  type(ocean_density_type), intent(in)           :: Dens
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error in ocean_density_mod (ocean_density_end): module must be initialized')
  endif 

  if(.not. write_a_restart) then
    write(stdoutunit,'(/a)') '==>Warning from ocean_density_mod (ocean_density_end): NO restart written.'
    call mpp_error(WARNING,'==>Warning from ocean_density_mod (ocean_density_end): NO restart written.')
    return
  endif 

  call ocean_density_restart(Time, Dens)

  write(stdoutunit,'(/a)') ' From ocean_density_mod: ending density chksums'
  call write_timestamp(Time%model_time)
  call ocean_density_chksum(Time, Dens)

  module_is_initialized = .FALSE.

end subroutine ocean_density_end
! </SUBROUTINE> NAME="ocean_density_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_density_restart(Time, Dens, time_stamp)
  type(ocean_time_type),    intent(in)           :: Time
  type(ocean_density_type), intent(in)           :: Dens
  character(len=*),         intent(in), optional :: time_stamp
   integer :: tau, taup1

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error in ocean_density_mod (ocean_density_end): module must be initialized')
  endif 

  taup1  = Time%taup1

  call reset_field_pointer(Den_restart, id_restart_rho, Dens%rho(:,:,:,taup1) )

  call save_restart(Den_restart, time_stamp)

end subroutine ocean_density_restart
! </SUBROUTINE> NAME="ocean_density_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_chksum">
!
! <DESCRIPTION>
! Compute checksums for density. 
! </DESCRIPTION>
subroutine ocean_density_chksum(Time, Dens)

  type(ocean_time_type),    intent(in) :: Time
  type(ocean_density_type), intent(in) :: Dens
  integer(i8_kind)                     :: chk_sum

  integer :: taup1

  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_density_mod (ocean_density_chksum): module not yet initialized ')
  endif 

  taup1 = Time%taup1 

  wrk1(isc:iec,jsc:jec,:) = Dens%rho(isc:iec,jsc:jec,:,taup1)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum                 = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'chksum for rho(taup1)        = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Dens%pressure_at_depth(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum                 = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'chksum for pressure_at_depth = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = denominator_r(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum                 = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'chksum for denominator_r     = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Dens%drhodT(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum                 = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'chksum for drhodT            = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Dens%drhodS(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum                 = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'chksum for drhodS            = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Dens%drhodz_zt(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum                 = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'chksum for drhodz_zt         = ',  chk_sum

end subroutine ocean_density_chksum
! </SUBROUTINE>  NAME="ocean_density_chksum"


end module ocean_density_mod





