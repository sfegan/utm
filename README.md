# Latitude/Longitude to UTM/UPS

This repository provides a set of C++ routines to convert coordinates between latitude/longitude and the UTM/UPS grids.

__Prior to 2025-12-21__ these were based on the series expansion provided by the US Defense Mapping Agency Technical Manual 8358.2.

__After 2025-12-21__ the UTM routines use an implementation based on the expansions of [Karney 2011](https://arxiv.org/abs/1002.1417), [Kawase 2011](http://www.gsi.go.jp/common/000062452.pdf), and [Kawase 2013](http://www.gsi.go.jp/common/000065826.pdf). The original DMATM implementation is also present for reference.

## Building and testing

The code for the transformations is contained in the routines `utm.cpp` and `datum.cpp`. A test function is also included which can be compiled by setting a flag.

To build the UTM test suite on Linux or macOS use:

    g++ -o test_ellipse test_ellipse.cpp datum.cpp utm.cpp
    ./test_ellipse

This will output a series of forward and backward conversions between latitude/longitude and the UTM grid, reproducing table 2-11 of DMATM 8358.2 (which is included in this repository).

## Function reference

A first set of functions are provided to convert between geographic coordinates and transverse Mercator (TM) or polar stereographic projections. For the TM case these can be used with respect to any reference meridian and scale factor, and applying any desired false northing and false easting. 

The second set of functions apply the standard meridians, scale factor, false northing, and false easting as defined by DMATM 8358.2 to produce conversions between geographic coordinates and the UTM/UPS grids. In the forward conversion the user can allow the grid zone and hemisphere to be chosen automatically based on the geographic coordinates, or can impose any zone, allowing coordinates near the limits of one zone to be mapped into the coordinates of another, for example.

### Enums

- `GridZone`: Enumerates grid zone options:
    * `UTM_ZONE_1` to `UTM_ZONE_60`: Specific UTM zones (zone 1 covers -180° to -174°, zone 60 covers 174° to 180°).
    * `UPS_NORTH`: Polar Stereographic projection for the North Pole.
    * `UPS_SOUTH`: Polar Stereographic projection for the South Pole.
    * `GRID_AUTO`: Automatically selects the appropriate grid zone and projection based on latitude and longitude. For latitudes ≥84° or <-80°, selects UPS; otherwise selects the appropriate UTM zone (1-60) based on longitude, with special handling for Norway and Svalbard. __This is the recommended setting to convert geographic coordinates to the grid in conformance to DMATM 8358.2__.
    * `UTM_ZONE_AUTO`: Automatically selects the appropriate UTM zone (1-60) based on longitude alone, ignoring latitude. With this option the polar regions are mapped onto the various UTM zones. This option is primarily intended for internal use but can be selected by the user if desired.

- `Hemisphere`: Enumerates hemisphere options:
    * `HEMI_NORTH`: Northern Hemisphere.
    * `HEMI_SOUTH`: Southern Hemisphere.
    * `HEMI_AUTO`: Automatically selects based on latitude (NORTH for ≥0°, SOUTH for <0°).

### Spherical Projections

- `geographic_to_tm_sphere(double R, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic coordinates (latitude, longitude in radians) to Transverse Mercator (TM) coordinates on a sphere. The function takes the following parameters:
  * `R`: sphere radius in meters (input),
  * `k0`: scale factor (input),
  * `lon_mer`: central meridian (input),
  *  `FN` and `FE`: false northing and easting in meters (input),
  * `lat_rad`: latitude in radians (input),
  * `lon_rad`: longitude in radians (input),
  * `N`: calculated northing in meters (output),
  * `E`: calculated easting in meters (output).
- `geographic_to_tm_sphere_with_convergence_and_scale(double R, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E, double* grid_convergence_rad, double* scale)`: Converts geographic coordinates to TM coordinates on a sphere and computes grid convergence and scale. Parameters as above, with additional outputs:
  * `grid_convergence_rad`: angle between lines of constant easting and true North in radians (output),
  * `scale`: scale on the grid (output).
- `tm_to_geographic_sphere(double R, double k0, double lon_mer, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts TM coordinates to geographic coordinates on a sphere. The parameters are as above, but with `N` and `E` as inputs, and `lat_rad` and `lon_rad` as output.
- `geographic_to_ps_sphere(double R, double k0, Hemisphere hemi, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic coordinates to Polar Stereographic (PS) coordinates on a sphere. Parameters as `geographic_to_tm_sphere`, but without `lon_mer` and with the addition of:
  * `hemi`: the hemisphere with respect to which the projection will be done.
- `geographic_to_ps_sphere_with_convergence_and_scale(double R, double k0, Hemisphere hemi, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E, double* grid_convergence_rad, double* scale)`: Converts geographic coordinates to PS coordinates on a sphere and computes grid convergence and scale. Parameters as above, with additional outputs:
  * `grid_convergence_rad`: angle between lines of constant easting and true North in radians (output),
  * `scale`: scale on the grid (output).
- `ps_to_geographic_sphere(double R, double k0, Hemisphere hemi, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts PS coordinates to geographic coordinates on a sphere. Parameters as above but with `N` and `E` as inputs, and `lat_rad` and `lon_rad` as outputs.

### Ellipsoidal Projections

- `geographic_to_tm(double a, double e2, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic to TM on the ellipsoid using the Karney/Kawase expansions. Parameters as in the case of `geographic_to_tm_sphere` but with `R` replaced by: 
  * `a`: semi-major axis in meters, and 
  * `e2`: eccentricity squared.
- `geographic_to_tm_with_convergence_and_scale(double a, double e2, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E, double* grid_convergence_rad, double* scale)`: Converts geographic to TM on the ellipsoid using the Karney/Kawase expansions and computes the grid convergence (in radians) and scale factor. The parameters to the function are as above, with the addition of the following outputs:
  * `grid_convergence_rad`: angle between lines of constant easting and true North in radians at the chosen point (output), and
  * `scale`: scale on the grid at the chosen point (output).
- `tm_to_geographic(double a, double e2, double k0, double lon_mer, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts TM to geographic on the elliposid using the Karney/Kawase expansions. Parameters as in `geographic_to_tm` with `N` and `E` as inputs, and `lat_rad` and `lon_rad` as outputs.
- `geographic_to_ps(double a, double e2, double k0, Hemisphere hemi, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic to PS on the ellipsoid using the DMATM algorithm.
- `geographic_to_ps_with_convergence_and_scale(double a, double e2, double k0, Hemisphere hemi, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E, double* grid_convergence_rad, double* scale)`: Converts geographic to PS on the ellipsoid using the DMATM algorithm and computes grid convergence and scale. Parameters as above, with additional outputs:
  * `grid_convergence_rad`: angle between lines of constant easting and true North in radians (output),
  * `scale`: scale on the grid (output).
- `ps_to_geographic(double a, double e2, double k0, Hemisphere hemi, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts PS to geographic on an ellipsoid using the DMATM algorithm.

### Ellipsoidal Projections (obsoleted DMATM expansions)

- `dmatm_geographic_to_tm(double a, double e2, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic to TM on the ellipse, using DMATM series expansion. Parameters as in the case of `geographic_to_tm`.
- `dmatm_tm_to_geographic(double a, double e2, double k0, double lon_mer, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts TM to geographic using DMATM series expansion. Parameters as in `tm_to_geographic`.

### Grid Conversions

- `geographic_to_grid(double a, double e2, double lat_rad, double lon_rad, GridZone* zone, Hemisphere* hemi, double* N, double* E, double* grid_convergence_rad=nullptr, double* scale=nullptr)`: Converts geographic coordinates to UTM/UPS grid coordinates. Automatically selects zone/hemisphere if set to AUTO. Returns 1 on success. The function takes the following parameters:
  * `a`: semi-major axis of the ellipsoid in meters (input),
  * `e2`: eccentricity squared (input),
  * `lat_rad`: latitude in radians (input),
  * `lon_rad`: longitude in radians (input),
  * `zone`: pointer to the grid zone (input/output). Set to `GRID_AUTO` to automatically select the appropriate UTM zone based on longitude (1-60) or UPS (for latitudes >=84° or < -80°), and the selected zone will be returned in this variable. Can be set to a specific UTM zone (1-60) or `UPS_NORTH`/`UPS_SOUTH` to force conversion into that zone, allowing coordinates near zone boundaries to be mapped into a different zone,
  * `hemi`: pointer to the hemisphere (input/output). Set to `HEMI_AUTO` to automatically select based on latitude (NORTH for >=0°, SOUTH for <0°), and the selected hemisphere will be returned in this variable. Can be set to `HEMI_NORTH` or `HEMI_SOUTH` to force a specific hemisphere,
  * `N`: pointer to the calculated northing in meters (output),
  * `E`: pointer to the calculated easting in meters (output).
  * `grid_convergence_rad`: optional pointer to the grid convergence in radians (output). Both this and `scale` must be non-null for convergence and scale to be computed,
  * `scale`: optional pointer to the scale factor (output). Both this and `grid_convergence_rad` must be non-null for convergence and scale to be computed.
- `grid_to_geographic(double a, double e2, GridZone zone, Hemisphere hemi, double N, double E, double* lat_rad, double* lon_rad)`: Converts UTM/UPS grid coordinates to geographic. Returns 1 on success. The function takes the following parameters:
  * `a`: semi-major axis of the ellipsoid in meters (input),
  * `e2`: eccentricity squared (input),
  * `zone`: grid zone (input),
  * `hemi`: hemisphere (input),
  * `N`: northing in meters (input),
  * `E`: easting in meters (input),
  * `lat_rad`: pointer to the calculated latitude in radians (output),
  * `lon_rad`: pointer to the calculated longitude in radians (output).

