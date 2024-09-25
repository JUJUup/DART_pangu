! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

use        types_mod, only : r8, i8, MISSING_R8

use time_manager_mod, only : time_type, set_time

use     location_mod, only : location_type, get_close_type, get_location, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, IS_VERTICAL, &
                             VERTISLEVEL, VERTISPRESSURE, VERTISSURFACE

use    utilities_mod, only : error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode

use   state_structure_mod,only : add_domain, get_dart_vector_index, get_domain_size, &
                                 get_dim_name, get_kind_index, get_num_dims, &
                                 get_num_variables, &
                                 get_model_variable_indices, state_structure_info

use ensemble_manager_mod, only : ensemble_type

use          obs_kind_mod, only: QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, &
                                 QTY_SURFACE_PRESSURE, QTY_TEMPERATURE, &
                                 QTY_PRESSURE, QTY_SPECIFIC_HUMIDITY, &
                                 QTY_10M_U_WIND_COMPONENT, QTY_10M_V_WIND_COMPONENT, QTY_2M_TEMPERATURE, &
                                 get_index_for_quantity
! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, read_model_time, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step

use distributed_state_mod, only : get_state, get_state_array

           
implicit none
private

! routines required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          adv_1step,              &
          init_time,              &
          init_conditions,        &
          shortest_time_between_assimilations, &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
"$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=32 ), parameter :: version = "$Revision$"
character(len=128), parameter :: tag = "$Id$"

logical :: module_initialized = .false.
integer :: dom_id ! used to access the state structure
type(time_type) :: assimilation_time_step 

! Example Namelist
! Use the namelist for options to be set at runtime.
character(len=256) :: template_file = 'model_restart.nc'
integer  :: time_step_days      = 0
integer  :: time_step_seconds   = 3600
integer, parameter :: max_state_variables = 100
integer, parameter :: num_state_table_columns = 2

integer, allocatable :: state_qty_list(:)

! this defines what variables are in the state, and in what order.
character(len=256) :: state_variables(num_state_table_columns,max_state_variables) = 'NULL'

! this defines pangu grid 
integer   :: nlat   = 721
integer   :: nlon   = 1440
integer   :: nvert   = 13
real(r8) :: grid_size = 0.25
real(r8), allocatable :: lats(:),lons(:),verts(:)
real(r8), dimension(13) :: pressure_level = (/100000.0, 92500.0, 85000.0, 70000.0, 60000.0, 50000.0, 40000.0, 30000.0, 25000.0, 20000.0, 15000.0, 10000.0, 5000.0/)
! real(r8), dimension(:), pointer :: lats,lons,verts

! (/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0/)
type pangu_grid
    real(r8), allocatable  :: lats(:), lons(:), verts(:)
    integer                :: nlat, nlon, nvert
end type pangu_grid

type(pangu_grid)           :: grid_data

namelist /model_nml/ nlat, nlon, nvert, grid_size, &
                    template_file, state_variables
! namelist /model_nml/ nlat, nlon, nvert, grid_size, pressure_level, &
!                       template_file, state_variables

character(len=256) :: errstring

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

integer  :: iunit, io
integer :: maxrows, i, numrows, rows, this_qty


module_initialized = .true.
! write(*, *) 'stage before reading namelist'

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime 
! feel free to hardcode it and remove from the namelist.

assimilation_time_step = set_time(time_step_seconds, &
                                  time_step_days)


! init pangu grid to define lat/lon/vert grids.
! write(*, *) 'stage before init_pangu_grid'
call init_pangu_grid(nlat,nlon,nvert,grid_size,pressure_level)
! Define which variables are in the model state
! loop around the variables - the number of rows is the number of
! fields in the state vector.
maxrows = size(state_variables, 2)

numrows = 0
COUNTROWS: do i=1, maxrows

   !>@todo leave a comment here about what 2 is for, once we figure it out
   !> shouldn't it be either 1, or both 1 and 2?
   if (state_variables(2, i) == 'NULL') exit COUNTROWS

   numrows = i
enddo COUNTROWS

if (numrows == 0) then
   call error_handler(E_ERR,'fill_domain_structure', "bad namelist: no variables found in namelist item 'state_variables'", &
                      source, revision, revdate)
endif

allocate(state_qty_list(numrows))
do i = 1, numrows
   state_qty_list(i) = get_index_for_quantity(state_variables(2,i))
   if (state_qty_list(i) < 0) then
      call error_handler(E_ERR,'fill_domain_structure', "bad namelist: unknown kind: "//trim(state_variables(2,i)), &
                         source, revision, revdate)
   endif
end do

! need a template file.
if (template_file == 'null') then
    call error_handler(E_ERR,'template file', "need a template file for pangu model", &
                        source, revision, revdate)
else
    dom_id = add_domain(template_file, numrows, state_variables(1,1:numrows), state_qty_list(:))
endif

    
end subroutine static_init_model

subroutine init_pangu_grid(nlat,nlon,nvert,grid_size,pressure_level)

integer,             intent(in) :: nlat, nlon, nvert
real(r8),            intent(in) :: grid_size, pressure_level(nvert)
! real(r8), dimension(nlat):: lats
! real(r8), dimension(nlon):: lons
! real(r8), dimension(nvert):: verts
real(r8)                        :: iter=0.0 ! temp for iteration
integer                         :: i=1, j=1, k=1 ! for doing loops
real(r8)                        :: lat_lower=-90.0, lat_upper=90.0, lon_lower=0.0, lon_upper=359.99
! type(pangu_grid)                :: grid_data

allocate(grid_data % lats(1:nlat))
allocate(grid_data % lons(1:nlon))
allocate(grid_data % verts(1:nvert))

grid_data % nlat = nlat
grid_data % nlon = nlon
grid_data % nvert = nvert
DO i = 1,nlat
    iter = lat_lower + (i-1) * grid_size
    grid_data % lats(i) = iter
end DO

DO j = 1,nlon
    iter = lon_lower + (j-1) * grid_size
    grid_data % lons(j) = iter
end DO
DO k = 1,nvert
    grid_data % verts(k) = pressure_level(k)
end DO

! print *, grid_data % verts
allocate(lats(1:nlat))
allocate(lons(1:nlon))
allocate(verts(1:nvert))
lats = grid_data % lats
lons = grid_data % lons
verts = grid_data % verts
end subroutine init_pangu_grid

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = get_domain_size(dom_id)

end function get_model_size


!------------------------------------------------------------------
! Given a state handle, a location, and a state quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case a positive istatus should be returned.
!
! For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer :: num_lons, num_lats, lon_below, lon_above, lat_below, lat_above, i
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real(r8) :: lon_fract, lat_fract, temp_lon
real(r8) :: lon, lat, level, lon_lat_lev(3), pressure

integer :: tmp_status(ens_size,4), e
real(r8) :: val(2,2, ens_size), a(2, ens_size)

if ( .not. module_initialized ) call static_init_model

! Need to check that status of all four temp vals.
! Start out assuming that all pass.
tmp_status = 0

! Would it be better to pass state as prog_var_type (model state type) to here?
! As opposed to the stripped state vector. YES. This would give time interp.
! capability here; do we really want this or should it be pushed up?

! Get the position, determine if it is model level or pressure in vertical
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1); lat = lon_lat_lev(2); 
! write(*, *) "locatioin is: ", lon_lat_lev
if(is_vertical(location, "LEVEL")) then 
   level = lon_lat_lev(3)
else if(is_vertical(location, "PRESSURE")) then
   pressure = lon_lat_lev(3)
else if(is_vertical(location, "SURFACE")) then
   ! level is not used for surface pressure observations
   level = -1
else
   call error_handler(E_ERR,'model_interpolate', &
      'Pangu can only handle pressure or model level for obs vertical coordinate', &
      source, revision, revdate)
endif
! get appropriate lon and lat grid specs on mass grid
num_lons = size(lons)
num_lats = size(lats)
bot_lon = lons(1)
top_lon = lons(num_lons)
delta_lon = lons(2) - lons(1)
bot_lat = lats(1)
top_lat = lats(num_lats)
delta_lat = lats(2) - lats(1)

! Compute bracketing lon indices
if(lon >= bot_lon .and. lon <= top_lon) then
    lon_below = int((lon - bot_lon) / delta_lon) + 1
    lon_above = lon_below + 1
 !   write(*, *) 'lon, delta_lon, bot_lon', lon, delta_lon, bot_lon
 !   write(*, *) 'prod ', ((lon_below - 1) * delta_lon + bot_lon)
    lon_fract = (lon - ((lon_below - 1) * delta_lon + bot_lon)) / delta_lon
 else
 ! At wraparound point
    lon_below = num_lons
    lon_above = 1
    if(lon < bot_lon) then 
       temp_lon = lon + 360.0_r8
    else
       temp_lon = lon
    endif
    lon_fract = (temp_lon - top_lon) / delta_lon
 endif
 
! Next, compute neighboring lat rows
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE MAY BE WRONG
 if(lat >= bot_lat .and. lat <= top_lat) then
    lat_below = int((lat - bot_lat) / delta_lat) + 1
    lat_above = lat_below + 1
    lat_fract = (lat - ((lat_below - 1) * delta_lat + bot_lat)) / delta_lat
 else if(lat <= bot_lat) then
 ! South of bottom lat NEED TO DO BETTER: NOT REALLY REGULAR
    lat_below = 1
    lat_above = 2
    lat_fract = 0.0_r8
 else
 ! North of top lat NEED TO DO BETTER: NOT REALLY REGULAR
    lat_below = num_lats - 1
    lat_above = num_lats
    lat_fract = 1.0_r8
 endif

! Case 1: model level specified in vertical
if(is_vertical(location, "LEVEL") .or. is_vertical(location, "SURFACE")) then
    ! Now, need to find the values for the four corners
    val(1, 1,:) =  get_val(state_handle, ens_size, &
                        lon_below, lat_below, nint(level), qty)
    val(1, 2,:) =  get_val(state_handle, ens_size, &
                        lon_below, lat_above, nint(level), qty)
    val(2, 1,:) =  get_val(state_handle, ens_size, &
                        lon_above, lat_below, nint(level), qty)
    val(2, 2,:) =  get_val(state_handle, ens_size, &
                        lon_above, lat_above, nint(level), qty)
else
    ! Case of pressure specified in vertical
   ! call error_handler(E_MSG,'model_interpolate', &
   ! 'doing pressure interpolation. ')
    val(1, 1,:) =  get_val_pressure(state_handle, ens_size, &
                        lon_below, lat_below, pressure, qty, tmp_status(:,1))
    val(1, 2,:) =  get_val_pressure(state_handle, ens_size, &
                        lon_below, lat_above, pressure, qty, tmp_status(:,2))
    val(2, 1,:) =  get_val_pressure(state_handle, ens_size, &
                        lon_above, lat_below, pressure, qty, tmp_status(:,3))
    val(2, 2,:) =  get_val_pressure(state_handle, ens_size, &
                        lon_above, lat_above, pressure, qty, tmp_status(:,4))
endif
    
! Do the weighted average for interpolation
do i = 1, 2
    a(i,:) = lon_fract * val(2, i,:) + (1.0 - lon_fract) * val(1, i,:)
 end do
 
 expected_obs = lat_fract * a(2,:) + (1.0 - lat_fract) * a(1,:)
 
 istatus = maxval(tmp_status, 2)
 
 ! if the forward operater failed set the value to missing_r8
 do e = 1, ens_size
    if (istatus(e) /= 0) then
       expected_obs(e) = MISSING_R8
    endif
 enddo


! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
! expected_obs(:) = MISSING_R8

! istatus for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
! istatus(:) = 1

end subroutine model_interpolate

function get_val(state_handle, ens_size, lon_index, lat_index, level, itype)

    type(ensemble_type), intent(in) :: state_handle
    integer, intent(in) :: lon_index, lat_index, level, itype
    integer, intent(in) :: ens_size
    real(r8) :: get_val(ens_size)

    character(len = 129) :: msg_string
    integer :: model_type
    integer(i8) :: state_index
    ! Need to change the obs kind defined itype to the appropriate model type
    ! The itype_in variable uses types defined in the kinds module. The whole bgrid 
    ! model_mod should be modified to use this correctly. However, as a fast patch
    ! for the initial I-release, will just map the u, v, t, and ps kinds to the
    ! default expectations of the bgrid which are hard coded as u=1, v=2, ps =3,
    ! t = 4, and tracers for numbers larger than 4. For now, the tracer observations
    ! are not implemented.
    model_type = get_varid_from_kind(itype)

    ! Find the index into state array and return this value
    state_index = get_dart_vector_index(lon_index, lat_index, level, dom_id, model_type)
    get_val     = get_state(state_index, state_handle)

end function get_val
    

recursive function get_val_pressure(state_handle, ens_size, lon_index, lat_index, pressure, itype, istatus) result(val_pressure)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
real(r8),            intent(in) :: pressure
integer,             intent(in) :: lon_index
integer,             intent(in) :: lat_index
integer,             intent(in) :: itype
integer,            intent(out) :: istatus(:)

real(r8) :: val_pressure(ens_size)

type(location_type) :: ps_location
real(r8) :: ps(1, 1,ens_size), p_model(1, 1, ens_size), pfull(1, 1, grid_data%nvert), rfrac(ens_size)
integer :: top_lev(ens_size), bot_lev(ens_size), i, e, model_type
integer(i8) :: state_index_bottom(ens_size), state_index_top(ens_size)
real(r8) :: bot_val(ens_size), top_val(ens_size), ps_lon

! set all istatus to successful
istatus(:) = 0

! Need to get the surface pressure at this point.
! For t or tracers (on mass grid with ps) this is trivial
! For u or v (on velocity grid)

ps(1,1,:) = get_val(state_handle, ens_size, lon_index, lat_index, -1, QTY_SURFACE_PRESSURE)
! write(*, *) "get_val_pressure: stage0"
do e = 1, ens_size
   ! Next, get the values on the levels for this ps,
   ! haoxing: this need some thoughts.
   ! haoxing: lets just try predefine a value for pfull and see if it works.
   ! How to compute the pressure of pangu in each grid? see here for a reference: 
   ! https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
   pfull(1,1,:) = pressure_level
!    call compute_pres_full(Dynam%Vgrid, ps(:,:,e), pfull)
   ! pfull(1,1,:)=(/632.555794027661, 3860.37673931619, 9929.20528881130, 18230.0441531918, 28203.5500645921, 39291.6346158435, 50936.3558092015, 62579.7522262883, 73663.7797588585, 83630.2558959853, 91920.7637549028, 97976.5138133526, 101238.082113486/)

   ! Interpolate in vertical to get two bounding levels

   ! What to do about pressures above top??? Just use the top for now.
   ! Could extrapolate, but that would be very tricky. Might need to
   ! reject somehow.
   ! write(*, *) "pfull:", pfull
   ! write(*, *) "pressure:", pressure
   ! if(pressure > pfull(1, 1, 1)) then
   if(pressure > ps(1,1,e)) then
      top_lev(e) = 1
      bot_lev(e) = 2
      rfrac(e) = 1.0_r8
   ! Actually, just fail using istatus
      istatus(e) = 1
   else if(pressure < pfull(1, 1, grid_data%nvert)) then
   ! Same for bottom
      bot_lev(e) = grid_data%nvert
      top_lev(e) = bot_lev(e) - 1
      rfrac(e) = 0.0_r8
   ! Actually, just fail using istatus
      istatus(e) = 1
   else if (pressure < ps(1,1,e) .and. pressure > pfull(1,1,1)) then
      ! This case, pressure is between the surface pressure and first layer of upper pressure
      frac(e) = (ps(1,1,e) - pressure) / &
               (ps(1,1,e) - pfull(1, 1, i - top_lev(e)))
      ! write(*, *) "get_val_pressure: stage1"
      ! Search down through pressures

      do i = 2, grid_data%nvert
         if(pressure > pfull(1, 1, i)) then
               top_lev(e) = i + 1
               bot_lev(e) = i
               rfrac(e) = (pfull(1, 1, bot_lev(e)) - pressure) / &
               (pfull(1, 1, bot_lev(e)) - pfull(1, 1, top_lev(e)))
               ! write(*, *) "top_lev: ", top_lev
               ! write(*, *) "bot_lev: ", bot_lev
               ! write(*, *) "rfrac: ", rfrac

            goto 21
         endif
      end do
      ! write(*, *) "get_val_pressure: stage2"
   end if
   21 continue
enddo

model_type = get_varid_from_kind(itype)

do e = 1, ens_size
   ! Find the index into state array and return this value
   state_index_bottom(e) = get_dart_vector_index(lon_index, lat_index, bot_lev(e), dom_id, model_type)
   state_index_top(e)    = get_dart_vector_index(lon_index, lat_index, top_lev(e), dom_id, model_type)
   ! write(*, *) "state_index_bottom: ", state_index_bottom
   ! write(*, *) "state_index_top: ", state_index_top
enddo
! write(*, *) "get_val_pressure: stage3"
! write(*, *) "state_index_bottom: ", state_index_bottom
! write(*, *) "state_index_top: ", state_index_top

call get_state_array(bot_val, state_index_bottom, state_handle)
call get_state_array(top_val, state_index_top,    state_handle)
! write(*, *) "get_val_pressure: stage4"
do e = 1, ens_size
   val_pressure(e) = (1.0_r8 - rfrac(e)) * bot_val(e) + rfrac(e) * top_val(e)
enddo 

end function get_val_pressure


!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty

real(r8) :: lon, lat, lev
integer  :: i, j, k
integer  :: myvarid, myqty, nd, this_qty, vtype

if ( .not. module_initialized ) call static_init_model

! from dart index_in calculate local variable indicies i,j,k and variable id
! write(*, *) 'stage before get_model_variable_indices'
call get_model_variable_indices(index_in, i, j, k, var_id=myvarid, kind_index=myqty)
! should be set to the actual location using set_location()

! write(*, *) 'stage before judgement'
! ! write(*, *) lats,lons,verts
this_qty = state_qty_list(myvarid)
if (this_qty == QTY_TEMPERATURE) then
    lon = lons(i)
    lat = lats(j)
    lev = pressure_level(k)
    vtype = VERTISPRESSURE

 else if (this_qty == QTY_SURFACE_PRESSURE .or. this_qty == QTY_10M_U_WIND_COMPONENT .or. this_qty == QTY_10M_V_WIND_COMPONENT .or. this_qty == QTY_2M_TEMPERATURE) then
    lon = lons(i)
    lat = lats(j)
    lev = 1
    vtype = VERTISPRESSURE
   !  vtype = VERTISSURFACE
 else  ! the same as QTY_TEMPERATURE
    lon = lons(i)
    lat = lats(j)
    lev = k
    vtype = VERTISPRESSURE
 endif
 ! write(*, *) 'stage before set_location'
 location = set_location(lon, lat, lev, vtype)
 
! should be set to the physical quantity, e.g. QTY_TEMPERATURE
if (present(qty)) qty = myqty

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     ! location of interest
integer,                       intent(in)    :: base_type    ! observation TYPE
type(location_type),           intent(inout) :: locs(:)      ! state locations
integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
integer,                       intent(out)   :: num_close    ! how many are close
integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'


call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)


end subroutine get_close_state


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()


end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

function get_varid_from_kind(dart_kind)

    integer, intent(in) :: dart_kind
    integer             :: get_varid_from_kind
    
    integer :: i
    
    do i = 1, get_num_variables(dom_id)
       if (dart_kind == state_qty_list(i)) then
          get_varid_from_kind = i
          return
       endif
    end do
    
    write(errstring, *) 'Kind ', dart_kind, ' not found in state vector'
    call error_handler(E_ERR,'get_varid_from_kind', errstring, &
                       source, revision, revdate)
    
end function get_varid_from_kind 

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

