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
module land_model_mod

  use time_manager_mod, only: time_type

  use mpp_domains_mod,  only: domain1D, domain2d, mpp_get_layout, mpp_define_layout, &
                              mpp_define_domains, mpp_get_compute_domain,            &
                              CYCLIC_GLOBAL_DOMAIN, mpp_get_compute_domains,         &
                              mpp_get_domain_components, mpp_get_pelist

  use fms_mod,          only: write_version_number, error_mesg, FATAL, NOTE, &
                              mpp_pe, mpp_npes, mpp_root_pe, field_exist

  use mpp_mod,          only: mpp_error, FATAL

  use diag_manager_mod, only: send_data

  use mosaic_mod,       only: get_mosaic_ntiles, calc_mosaic_grid_area, &
                              get_mosaic_xgrid_size, get_mosaic_xgrid
  
  use constants_mod,    only: PI, RADIUS

implicit none

private

type land_data_type
   logical :: pe
   type(domain2d) :: domain  ! our computation domain
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size =>NULL(),       & ! fractional coverage of cell by tile, dimensionless
        t_surf =>NULL(),          & ! ground surface temperature, degK
        t_ca =>NULL(),            & ! canopy air temperature, degK
        albedo =>NULL(),          & ! snow-adjusted land albedo
        albedo_vis_dir =>NULL(),  & ! albedo for direct visible-band radiation
        albedo_nir_dir =>NULL(),  & ! albedo for direct nir-band radiation 
        albedo_vis_dif =>NULL(),  & ! albedo for diffuse visible-band radiation 
        albedo_nir_dif =>NULL(),  & ! albedo for diffuse nir-band radiation
        rough_mom =>NULL(),       & ! momentum roughness length, m
        rough_heat =>NULL(),      & ! roughness length for tracers (heat and water), m
        rough_scale =>NULL()        ! roughness length for drag scaling
   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity, kg/kg

   real, pointer, dimension(:,:) :: &  ! (lon, lat)
        discharge =>NULL(),       & ! flux from surface drainage network out of land model
        discharge_snow =>NULL(),  & ! snow analogue of discharge
	discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
	discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

   logical, pointer, dimension(:,:,:):: &
        mask =>NULL()                ! true if land

   integer :: axes(2)      ! axes IDs for diagnostics  
   logical, pointer, dimension(:,:) :: maskmap =>NULL() ! A pointer to an array indicating which
                                                        ! logical processors are actually used for
                                                        ! the ocean code. The other logical
                                                        ! processors would be all land points and
                                                        ! are not assigned to actual processors.
                                                        ! This need not be assigned if all logical
                                                        ! processors are used.
end type land_data_type

type atmos_land_boundary_type
   ! data passed from the atmosphere to the surface
   real, dimension(:,:,:,:), pointer :: tr_flux => NULL()
   real, dimension(:,:,:,:), pointer :: dfdtr   => NULL()
   real, dimension(:,:,:), pointer :: &
        t_flux =>NULL(),  &
        lw_flux =>NULL(), &
        lwdn_flux =>NULL(), &
        sw_flux =>NULL(), &
        sw_flux_down_vis_dir =>NULL(), &
        sw_flux_down_total_dir =>NULL(), &
        sw_flux_down_vis_dif =>NULL(), &
        sw_flux_down_total_dif =>NULL(), &
        lprec =>NULL(),   &
        fprec =>NULL(),   &
        tprec =>NULL()
   ! derivatives of the fluxes
   real, dimension(:,:,:), pointer :: &
        dhdt =>NULL(),    &
        dedt =>NULL(),    &
        dedq =>NULL(),    &
        drdt =>NULL()
   real, dimension(:,:,:), pointer :: &
        cd_m => NULL(),      &   ! drag coefficient for momentum, dimensionless
        cd_t => NULL(),      &   ! drag coefficient for tracers, dimensionless
        ustar => NULL(),     &   ! turbulent wind scale, m/s
        bstar => NULL(),     &   ! turbulent bouyancy scale, m/s
        wind => NULL(),      &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot => NULL(),     &   ! height of the bottom atmos. layer above surface, m
        drag_q =>NULL(),  &
        p_surf =>NULL()

   real, dimension(:,:,:), pointer :: &
        data =>NULL() ! collective field for "named" fields above
   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type

! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart

public update_land_model_fast   ! fast time-scale update of the land
public update_land_model_slow   ! slow time-scale update of the land
public send_averaged_data       ! send tile-averaged data to diag manager
public Lnd_stock_pe
! ==== end of public interfaces ==============================================

! ==== public data type =====================================================
public land_data_type, atmos_land_boundary_type


! ==== some names, for information only ======================================
logical                       :: module_is_initialized = .FALSE.
character(len=*),   parameter :: module_name = 'land_mod'
character(len=128), parameter :: version     = ''
character(len=128), parameter :: tagname     = ''

! ==== local module variables ================================================
integer            :: n_tiles = 1  ! number of tiles

! ---- interface to tile-averaged diagnostic routines ------------------------
interface send_averaged_data
   module procedure send_averaged_data2d
   module procedure send_averaged_data3d
end interface

contains

! ============================================================================
subroutine land_model_init &
     ( atmos2land, Land_bnd, time_init, time, dt_fast, dt_slow, atmos_domain)

  use fms_mod, only : field_size, read_data
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains,&
                              CYCLIC_GLOBAL_DOMAIN, mpp_get_data_domain
  use diag_manager_mod, only : diag_axis_init
! ============================================================================
! initialiaze land model using grid description file as an input. This routine
! reads land grid boundaries and area of land from a grid description file

! NOTES: theoretically, the grid description file can specify any regular
! rectangular grid for land, not jast lon/lat grid. Therefore the variables
! "xbl" and "ybl" in NetCDF grid spec file are not necessarily lon and lat
! boundaries of the grid.
!   However, at this time the module land_properties assumes that grid _is_
! lon/lat and therefore the entire module also have to assume that the land 
! grid is lon/lat.
!   lon/lat grid is also assumed for the diagnostics, but this is probably not
! so critical. 

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout) :: atmos2land ! land boundary data
  type(land_data_type), intent(inout) :: Land_bnd ! land boundary data
  type(time_type), intent(in)    :: time_init ! initial time of simulation (?)
  type(time_type), intent(in)    :: time      ! current time
  type(time_type), intent(in)    :: dt_fast   ! fast time step
  type(time_type), intent(in)    :: dt_slow   ! slow time step
  type(domain2d),  intent(in), optional :: atmos_domain ! domain of computations

  real, dimension(:),   allocatable :: glonb, glatb, glon, glat
  real, dimension(:),   allocatable :: xgrid_area
  integer, dimension(:),allocatable :: i1, j1, i2, j2
  real, dimension(:,:), allocatable :: gfrac, garea, tmpx, tmpy, geo_lonv, geo_latv
  integer                           :: siz(4), layout(2), is, ie, js, je, i, j, k, n, m
  integer                           :: nlon, nlat, nlonb, nlatb, ntiles, nfile_axl, nxgrid 
  character(len=256)                :: grid_file   = "INPUT/grid_spec.nc"
  character(len=256)                :: land_mosaic, tile_file, axl_file

  if( field_exist(grid_file, 'AREA_LND') ) then
     call field_size( grid_file, 'AREA_LND', siz)
     nlon = siz(1)
     nlat = siz(2)
     nlonb = nlon + 1
     nlatb = nlat + 1
     ! allocate data for longitude and latitude bounds
     allocate( glonb(nlonb), glatb(nlatb), glon(nlon), glat(nlat))
     allocate( garea(nlon,nlat), gfrac(nlon,nlat) )
     ! read coordinates of grid cell vertices
     call read_data(grid_file, "xbl", glonb, no_domain=.true.)
     call read_data(grid_file, "ybl", glatb, no_domain=.true.)

     ! read coordinates of grid cell centers
 !    call read_data(grid_file, "xtl", glon, no_domain=.true.)
 !    call read_data(grid_file, "ytl", glat, no_domain=.true.)
     glon(1:siz(1)) = (glonb(1:siz(1))+glonb(2:siz(1)+1))/2.0
     glat(1:siz(2)) = (glatb(1:siz(2))+glatb(2:siz(2)+1))/2.0

     ! read land area, land cell area and calculate the fractional land coverage 
     call read_data(grid_file, "AREA_LND_CELL", garea, no_domain=.true.)
     call read_data(grid_file, "AREA_LND",      gfrac,     no_domain=.true.)
     ! calculate land fraction
     gfrac     = gfrac/garea
  else if( field_exist(grid_file, 'lnd_mosaic_file') ) then ! read from mosaic file
     call read_data(grid_file, "lnd_mosaic_file", land_mosaic)     
     land_mosaic = "INPUT/"//trim(land_mosaic)
     ntiles = get_mosaic_ntiles(land_mosaic)
     if(ntiles .NE. 1) call error_mesg('land_model_init',  &
         ' ntiles should be 1 for land mosaic, contact developer', FATAL)
     call read_data(land_mosaic, "gridfiles", tile_file )
     tile_file = 'INPUT/'//trim(tile_file)
     call field_size(tile_file, "x", siz)
     if( mod(siz(1)-1,2) .NE. 0) call error_mesg("land_model_init", "size(x,1) - 1 should be divided by 2", FATAL);
     if( mod(siz(2)-1,2) .NE. 0) call error_mesg("land_model_init", "size(x,2) - 1 should be divided by 2", FATAL);
     nlon = (siz(1)-1)/2
     nlat = (siz(2)-1)/2
     nlonb = nlon + 1
     nlatb = nlat + 1
     allocate( glonb(nlonb), glatb(nlatb), glon(nlon), glat(nlat) )
     allocate( garea(nlon,nlat), gfrac(nlon,nlat) )
     !--- read the grid information on supergrid.
     allocate( tmpx(2*nlon+1, 2*nlat+1), tmpy(2*nlon+1, 2*nlat+1) )
     call read_data(tile_file, "x", tmpx, no_domain=.TRUE.)
     call read_data(tile_file, "y", tmpy, no_domain=.TRUE.)
     !--- make sure the grid is regular lat-lon grid.
     do j = 1, 2*nlat+1
        do i = 2, 2*nlon+1
           if(tmpy(i,j) .NE. tmpy(1,j)) call mpp_error(FATAL, "land_model_init:longitude is not uniform")
        end do
     end do
     do i = 1, 2*nlon+1
        do j = 2, 2*nlat+1
           if(tmpx(i,j) .NE. tmpx(i,1)) call mpp_error(FATAL, "land_model_init:latitude is not uniform")
        end do
     end do
     !--- make sure the grid is regular lat-lon grid.

     do i = 1, nlonb
        glonb(i) = tmpx(2*i-1,1)
     end do     
     do j = 1, nlatb
        glatb(j) = tmpy(1, 2*j-1)
     end do     
     do i = 1, nlon
        glon(i) = tmpx(2*i,1)
     end do     
     do j = 1, nlat
        glat(j) = tmpy(1, 2*j)
     end do     
     !--- land_cell_area will be calculated using the same way to calculate the area of xgrid.
     allocate(geo_lonv(nlonb,nlatb), geo_latv(nlonb,nlatb))
     do j = 1, nlatb
        do i = 1, nlonb
           geo_lonv(i,j) = tmpx(2*i-1,2*j-1)*pi/180.0
           geo_latv(i,j) = tmpy(2*i-1,2*j-1)*pi/180.0
        end do
     end do

     call calc_mosaic_grid_area(geo_lonv, geo_latv, garea)
     deallocate(tmpx, tmpy, geo_lonv, geo_latv)

     !--- land area will be calculated based on exchange grid area.
     call field_size(grid_file, "aXl_file", siz)
     nfile_axl = siz(2)
     gfrac = 0
     do n = 1, nfile_axl
        call read_data(grid_file, "aXl_file", axl_file, level=n)
        axl_file = 'INPUT/'//trim(axl_file)
        nxgrid = get_mosaic_xgrid_size(axl_file)
        allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), xgrid_area(nxgrid))
        call get_mosaic_xgrid(aXl_file, i1, j1, i2, j2, xgrid_area)
        do m = 1, nxgrid
           i = i2(m); j = j2(m)
           gfrac(i,j) = gfrac(i,j) + xgrid_area(m)
        end do
        deallocate(i1, j1, i2, j2, xgrid_area)
     end do     
     gfrac = gfrac*4*PI*RADIUS*RADIUS/garea
  else
     call error_mesg('land_model_init','both AREA_LND and lnd_mosaic_file do not exist in file '//trim(grid_file),FATAL)
  end if  
  
  if( ASSOCIATED(Land_bnd%maskmap) ) then
     layout(1) = size(Land_bnd%maskmap,1)
     layout(2) = size(Land_bnd%maskmap,2)
     call mpp_define_domains((/1,nlon,1,nlat/), layout, Land_bnd%domain, &
        xflags = CYCLIC_GLOBAL_DOMAIN, maskmap = Land_bnd%maskmap, name='land model' )
  else
     call mpp_define_layout((/1,nlon,1,nlat/), mpp_npes(), layout)
     call mpp_define_domains((/1,nlon,1,nlat/), layout, Land_bnd%domain, &
        xflags = CYCLIC_GLOBAL_DOMAIN, name='land model')
  end if

  call mpp_get_data_domain(Land_bnd%domain,is,ie,js,je)

  allocate ( &
       Land_bnd % tile_size      (is:ie,js:je,n_tiles), & 
       Land_bnd % t_surf         (is:ie,js:je,n_tiles), &
       Land_bnd % t_ca           (is:ie,js:je,n_tiles), &
       Land_bnd % tr             (is:ie,js:je,n_tiles,1), &     ! one tracer for q?
       Land_bnd % albedo         (is:ie,js:je,n_tiles), & 
       Land_bnd % albedo_vis_dir (is:ie,js:je,n_tiles), &
       Land_bnd % albedo_nir_dir (is:ie,js:je,n_tiles), &
       Land_bnd % albedo_vis_dif (is:ie,js:je,n_tiles), &
       Land_bnd % albedo_nir_dif (is:ie,js:je,n_tiles), &
       Land_bnd % rough_mom      (is:ie,js:je,n_tiles), & 
       Land_bnd % rough_heat     (is:ie,js:je,n_tiles), & 
       Land_bnd % rough_scale    (is:ie,js:je,n_tiles), & 
       Land_bnd % discharge      (is:ie,js:je),   & 
       Land_bnd % discharge_snow (is:ie,js:je),   &
       Land_bnd % discharge_heat (is:ie,js:je),   &
       Land_bnd % discharge_snow_heat (is:ie,js:je),   &
       Land_bnd % mask       (is:ie,js:je,n_tiles)    )
  
  do i=is,ie
     do j=js,je
        do k= 1, n_tiles
           if (gfrac(i,j) > 0.0) then
               Land_bnd%tile_size(i,j,k)  = 1.0/n_tiles
               Land_bnd%mask(i,j,k) = .true.
           else
               Land_bnd%tile_size(i,j,k) = 0.0
               Land_bnd%mask(i,j,k) = .false.
           endif
        enddo
     enddo
  enddo

  Land_bnd%t_surf = 273.0
  Land_bnd%t_ca = 273.0
  Land_bnd%tr = 0.0
  Land_bnd%albedo = 0.0
  Land_bnd % albedo_vis_dir = 0.0
  Land_bnd % albedo_nir_dir = 0.0
  Land_bnd % albedo_vis_dif = 0.0
  Land_bnd % albedo_nir_dif = 0.0
  Land_bnd%rough_mom = 0.01
  Land_bnd%rough_heat = 0.01
  Land_bnd%rough_scale = 1.0
  Land_bnd%discharge = 0.0
  Land_bnd%discharge_snow = 0.0
  Land_bnd%discharge_heat = 0.0
  Land_bnd%discharge_snow_heat = 0.0
  Land_bnd%mask = .true.
  
  Land_bnd%axes(1) = diag_axis_init('lon',glon,'degrees_E','X','longitude',&
       set_name='land',domain2 = Land_bnd%domain)

  Land_bnd%axes(2) = diag_axis_init('lat',glon,'degrees_N','Y','latitude',&
       set_name='land',domain2 = Land_bnd%domain)  
  
  allocate( atmos2land % t_flux  (is:ie,js:je,n_tiles) )
  allocate( atmos2land % tr_flux  (is:ie,js:je,n_tiles,1) )     ! one tracer for q?
  allocate( atmos2land % dfdtr    (is:ie,js:je,n_tiles,1) )     ! one tracer for q?
  allocate( atmos2land % lw_flux (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux (is:ie,js:je,n_tiles) )
  allocate( atmos2land % lprec   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % fprec   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % tprec   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % dhdt    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % dedt    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % dedq    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % drdt    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % drag_q  (is:ie,js:je,n_tiles) )
  allocate( atmos2land % p_surf  (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_vis_dir   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_total_dir (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_vis_dif   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_total_dif (is:ie,js:je,n_tiles) )

  return

end subroutine land_model_init






! ============================================================================
subroutine land_model_end ( atmos2land, bnd )
! ============================================================================
! destruct the land model data

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd

  module_is_initialized = .FALSE.

!  deallocate boundary exchange data
!  call deallocate_boundary_data ( bnd )
  
end subroutine land_model_end

subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp ! timestamp to add to the file name

end subroutine land_model_restart

! ============================================================================
subroutine update_land_model_fast ( atmos2land, bnd )
! ============================================================================
! updates state of the land model on the fast time scale

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout)    :: atmos2land
  type(land_data_type),  intent(inout) :: bnd ! state to update

  return

end subroutine update_land_model_fast



! ============================================================================
subroutine update_land_model_slow ( atmos2land, bnd )
! ============================================================================
! updates land on slow time scale

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd


  return
  
end subroutine update_land_model_slow


! ============================================================================
subroutine update_land_bnd_fast ( bnd )
! ============================================================================
! updates land boundary data (the ones that atmosphere sees) on the fast time 
! scale this routine does not update tiling structure, because it is assumed
! that the tiling does not change on fast time scale

  ! ---- arguments -----------------------------------------------------------
  type(land_data_type), intent(inout) :: bnd

  return

end subroutine update_land_bnd_fast


! ============================================================================
subroutine update_land_bnd_slow ( bnd )
! ============================================================================
! updates land boundary data for the atmosphere on the slow time scale. This
! subroutine is responsible for the changing of tiling structure, if necessary,
! as well as for changing albedo, drag coefficients and such.

! NOTE: if tiling structure has been modified, then probably the distribution
! of other boundary values, such as temperature and surface humidity, should be
! modified too, not yet clear how.

  ! ---- arguments -----------------------------------------------------------
  type(land_data_type), intent(inout) :: bnd

  return

end subroutine update_land_bnd_slow

! ============================================================================
function send_averaged_data2d ( id, field, area, time, mask )
! ============================================================================
! average the data over tiles and send then to diagnostics

  ! --- return value ---------------------------------------------------------
  logical                      :: send_averaged_data2d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id             ! id od the diagnostic field 
  real,    intent(in)          :: field(:,:,:)   ! field to average and send
  real,    intent(in)          :: area (:,:,:)   ! area of tiles (== averaging 
                                                 ! weights), arbitrary units
  type(time_type), intent(in)  :: time           ! current time
  logical, intent(in),optional :: mask (:,:,:)   ! land mask

  ! --- local vars -----------------------------------------------------------
  real  :: out(size(field,1), size(field,2))

  call average_tiles( field, area, mask, out )
  send_averaged_data2d = send_data( id, out, time, mask=ANY(mask,DIM=3) )
end function send_averaged_data2d


! ============================================================================
function send_averaged_data3d( id, field, area, time, mask )
! ============================================================================
! average the data over tiles and send then to diagnostics

  ! --- return value ---------------------------------------------------------
  logical                      :: send_averaged_data3d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id              ! id of the diagnostic field
  real,    intent(in)          :: field(:,:,:,:)  ! (lon, lat, tile, lev) field 
                                                  ! to average and send
  real,    intent(in)          :: area (:,:,:)    ! (lon, lat, tile) tile areas 
                                                  ! ( == averaging weights), 
                                                  ! arbitrary units
  type(time_type), intent(in)  :: time            ! current time
  logical, intent(in),optional :: mask (:,:,:)    ! (lon, lat, tile) land mask

  ! --- local vars -----------------------------------------------------------
  real    :: out(size(field,1), size(field,2), size(field,4))
  logical :: mask3(size(field,1), size(field,2), size(field,4))
  integer :: it

  do it=1,size(field,4)
     call average_tiles( field(:,:,:,it), area, mask, out(:,:,it) )
  enddo

  mask3(:,:,1) = ANY(mask,DIM=3)
  do it = 2, size(field,4)
     mask3(:,:,it) = mask3(:,:,1)
  enddo

  send_averaged_data3d = send_data( id, out, time, mask=mask3 )
end function send_averaged_data3d

subroutine average_tiles ( x, area, mask, out )
! ============================================================================
! average 2-dimensional field over tiles
  ! --- arguments ------------------------------------------------------------
  real,    intent(in)  :: x   (:,:,:) ! (lon, lat, tile) field to average
  real,    intent(in)  :: area(:,:,:) ! (lon, lat, tile) fractional area
  logical, intent(in)  :: mask(:,:,:) ! (lon, lat, tile) land mask
  real,    intent(out) :: out (:,:)   ! (lon, lat)       result of averaging

  ! --- local vars -----------------------------------------------------------------
  integer  :: it                      ! iterator over tile number
  real     :: s(size(x,1),size(x,2))  ! area accumulator

  s(:,:)   = 0.0
  out(:,:) = 0.0

  do it = 1,size(area,3)
     where (mask(:,:,it)) 
        out(:,:) = out(:,:) + x(:,:,it)*area(:,:,it)
        s(:,:)   = s(:,:) + area(:,:,it)
     endwhere
  enddo

  where( s(:,:) > 0 ) &
       out(:,:) = out(:,:)/s(:,:)

end subroutine average_tiles

subroutine Lnd_stock_pe(bnd,index,value)
type(land_data_type), intent(in) :: bnd
integer, intent(in) :: index
real,    intent(out) :: value

value = 0.0

end subroutine Lnd_stock_pe

end module land_model_mod
   
