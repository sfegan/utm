# Latitude/Longitude to UTM/UPS

This repository provides a set of C++ routines to convert coordinates between latitude/longitude and the UTM/UPS grids.

__Prior to 2025-12-21__ these were based on the series expansion provided by the US Defense Mapping Agency Technical Manual 8358.2.

__After 2025-12-21__ the UTM routines use an implementation based on the expansions of [Karney 2011](https://arxiv.org/abs/1002.1417), [Kawase 2011](http://www.gsi.go.jp/common/000062452.pdf), and [Kawase 2013](http://www.gsi.go.jp/common/000065826.pdf). The original DMATM implementaion are also present for reference.

## Bulding and testing

The code for the transformations is contained in the routines `utm.cpp` and `datum.cpp`. A test function is also included which can be compiled by seetting a flag.

To build the UTM test suite on Linux or MacOS use:

    g++ -DELLIPSE_TEST_MAIN -o ellipse_test utm.cpp datum.cpp -lm
    ./ellipse_test

This will output a series of foward and backward conversions between latitude/longitude and the UTM grid, reproducing table 2-11 of DMATM 8358.2 (which is included in this repository).

## Function reference

### Enums

- `GridZone`: Enumerates UTM zones (1-60), UPS_NORTH, UPS_SOUTH, and GRID_AUTO for automatic selection.
- `Hemisphere`: Enumerates HEMI_NORTH, HEMI_SOUTH, and HEMI_AUTO for automatic selection.

### Spherical Projections

- `geographic_to_tm_sphere(double R, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic coordinates (latitude, longitude in radians) to Transverse Mercator (TM) coordinates on a sphere. Parameters: R (sphere radius), k0 (scale factor), lon_mer (central meridian), FN/FE (false northing/easting).
- `tm_to_geographic_sphere(double R, double k0, double lon_mer, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts TM coordinates to geographic coordinates on a sphere.
- `geographic_to_ps_sphere(double R, double k0, Hemisphere hemi, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic coordinates to Polar Stereographic (PS) coordinates on a sphere.
- `ps_to_geographic_sphere(double R, double k0, Hemisphere hemi, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts PS coordinates to geographic coordinates on a sphere.

### Ellipsoidal Projections (DMATM, obsoleted)

- `dmatm_geographic_to_tm(double a, double e2, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic to TM using DMATM series expansion. Parameters: a (semi-major axis), e2 (eccentricity squared).
- `dmatm_tm_to_geographic(double a, double e2, double k0, double lon_mer, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts TM to geographic using DMATM series expansion.

### Ellipsoidal Projections (Karney/Kawase, current)

- `geographic_to_tm(double a, double e2, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic to TM using Karney/Kawase expansions.
- `geographic_to_tm_with_convergence_and_scale(double a, double e2, double k0, double lon_mer, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E, double* grid_convergence_rad, double* scale)`: Converts geographic to TM and computes grid convergence (in radians) and scale factor.
- `tm_to_geographic(double a, double e2, double k0, double lon_mer, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts TM to geographic using Karney/Kawase expansions.
- `geographic_to_ps(double a, double e2, double k0, Hemisphere hemi, double FN, double FE, double lat_rad, double lon_rad, double* N, double* E)`: Converts geographic to PS on an ellipsoid.
- `ps_to_geographic(double a, double e2, double k0, Hemisphere hemi, double FN, double FE, double N, double E, double* lat_rad, double* lon_rad)`: Converts PS to geographic on an ellipsoid.

### Grid Conversions

- `geographic_to_grid(double a, double e2, double lat_rad, double lon_rad, GridZone* zone, Hemisphere* hemi, double* N, double* E)`: Converts geographic coordinates to UTM/UPS grid coordinates. Automatically selects zone/hemisphere if set to AUTO. Returns 1 on success.
- `grid_to_geographic(double a, double e2, GridZone zone, Hemisphere hemi, double N, double E, double* lat_rad, double* lon_rad)`: Converts UTM/UPS grid coordinates to geographic. Returns 1 on success.

