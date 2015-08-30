/* -*-mode:c++; mode:font-lock;-*- */

/******************************************************************************
 * 
 * CONVERT COORDINATES BETWEEN LATITUDE/LONGITUDE AND THE UTM/UPS GRIDS
 *
 * Stephen Fegan, July 2005, sfegan@gmail.com
 *
 * These conversion routines are a C/C++ implementation of the algorithms 
 * described in the Defense Mapping Agency Technical Manual (DMATM) 8358.2
 * which is available from the US National Geospatial Mapping Agency.
 * At time of writing it could be downloaded at:
 *
 * http://earth-info.nga.mil/GandG/coordsys/csat_pubs.html
 *
 * A number of alternative conversion routines are available on the Web. Those 
 * that I've seen are (ultimately) based on John Snyder's algorithm, presented
 * in "MAP PROJECTIONS; A WORKING MANUAL", USGS Professional Paper 1395. My
 * reading of that work suggests that Snyder's algorithms was derived as an
 * approximation to the original DMA algorithms.
 *
 * I have made every effort to make sure that this implementation is correct.
 * The "main" at the end of this file reproduces the worked examples in the
 * DMA document (section 2-11, page 2-7) to the accuracy of the equations
 * (0.01 meter on the grid and 0.001 arc second for geographic coordinates).
 * However, I DO NOT GIVE ANY ASSURANCES THAT THE OUTPUT OF THIS CODE IS 
 * CORRECT. If you need this code for mission critical applications it
 * is your responsibility to ensure its accuracy to the degree you require.
 *
 * The variable names may seem cryptic, but they were chosen to reflect the
 * DMA algorithm and to make the code writing easier.
 *
 * There only adjustable parameter is the tolerance to which the iteration 
 * to find the true meridional distance is performed in the reverse 
 * calculation. It is set at 0.001 meter which exceeds the accuracy the DMA 
 * claims across the range of applicability of the series. A lower value is
 * only meaningful if you are working close to the meridian.
 *
 * This code is in C++ but the main conversion routines should work without
 * modification under a C compiler. The "main" test code will need some
 * reworking to compile under C.
 *
 * I do not regard this work is anything other than a simple translation of 
 * the DMA algorithms, which were released as "Distribution Unlimited". 
 * Although I am a strong supporter of the GPL and free software, I do not 
 * think it would be appropriate to release this code under the GPL since 
 * all the hard work was done by the DMA and was unconditionally released to
 * the public. Therefore this code is released into the public domain.
 *
 * I request, however, that if you distribute the source, or a modification,
 * that you leave this header intact and that you leave the test "main" code
 * attached so that others can verify the code reproduces the DMA examples.
 *
 * $Id: utm.h,v 1.3 2006/11/30 21:37:33 sfegan Exp $
 *
 *****************************************************************************/

#ifndef UTM_H
#define UTM_H

enum GridZone
  {
    UTM_ZONE_AUTO = 0,
    UTM_ZONE_1=1,  UTM_ZONE_2,    UTM_ZONE_3,    UTM_ZONE_4,    UTM_ZONE_5,
    UTM_ZONE_6,    UTM_ZONE_7,    UTM_ZONE_8,    UTM_ZONE_9,    UTM_ZONE_10,
    UTM_ZONE_11,   UTM_ZONE_12,   UTM_ZONE_13,   UTM_ZONE_14,   UTM_ZONE_15,
    UTM_ZONE_16,   UTM_ZONE_17,   UTM_ZONE_18,   UTM_ZONE_19,   UTM_ZONE_20,
    UTM_ZONE_21,   UTM_ZONE_22,   UTM_ZONE_23,   UTM_ZONE_24,   UTM_ZONE_25,
    UTM_ZONE_26,   UTM_ZONE_27,   UTM_ZONE_28,   UTM_ZONE_29,   UTM_ZONE_30,
    UTM_ZONE_31,   UTM_ZONE_32,   UTM_ZONE_33,   UTM_ZONE_34,   UTM_ZONE_35,
    UTM_ZONE_36,   UTM_ZONE_37,   UTM_ZONE_38,   UTM_ZONE_39,   UTM_ZONE_40,
    UTM_ZONE_41,   UTM_ZONE_42,   UTM_ZONE_43,   UTM_ZONE_44,   UTM_ZONE_45,
    UTM_ZONE_46,   UTM_ZONE_47,   UTM_ZONE_48,   UTM_ZONE_49,   UTM_ZONE_50,
    UTM_ZONE_51,   UTM_ZONE_52,   UTM_ZONE_53,   UTM_ZONE_54,   UTM_ZONE_55,
    UTM_ZONE_56,   UTM_ZONE_57,   UTM_ZONE_58,   UTM_ZONE_59,   UTM_ZONE_60,
    UPS_NORTH, UPS_SOUTH, 
    GRID_AUTO
  };

enum Hemisphere
  {
    HEMI_AUTO = 0, HEMI_NORTH, HEMI_SOUTH
  };


#if !defined(__cplusplus)
typedef enum GridZone GridZone;
typedef enum Hemisphere Hemisphere;
#endif

/* FORWARD AND BACK TM/PS PROJECTIONS FOR A SPHERE */

void geographic_to_tm_sphere(double R, double k0, 
			     double lon_mer, double FN, double FE,
			     double lat_rad, double lon_rad,
			     double* N, double* E);

void tm_to_geographic_sphere(double R, double k0, 
			     double lon_mer, double FN, double FE,
			     double N, double E,
			     double* lat_rad, double* lon_rad);

void geographic_to_ps_sphere(double R, double k0, 
			     Hemisphere hemi, double FN, double FE,
			     double lat_rad, double lon_rad,
			     double* N, double* E);

void ps_to_geographic_sphere(double R, double k0, 
			     Hemisphere hemi, double FN, double FE,
			     double N, double E,
			     double* lat_rad, double* lon_rad);

/* FORWARD AND BACK TM/PS PROJECTIONS FOR AN ELLIPSOID */

#ifndef TM_TO_GEOGRAPHIC_TOLERANCE_M
#define TM_TO_GEOGRAPHIC_TOLERANCE_M 0.001
#endif

void geographic_to_tm(double a, double e2, double k0, 
		      double lon_mer, double FN, double FE,
		      double lat_rad, double lon_rad,
		      double* N, double* E);

void tm_to_geographic(double a, double e2, double k0, 
		      double lon_mer, double FN, double FE,
		      double N, double E,
		      double* lat_rad, double* lon_rad);

void geographic_to_ps(double a, double e2, double k0, 
		      Hemisphere hemi, double FN, double FE,
		      double lat_rad, double lon_rad,
		      double* N, double* E);

void ps_to_geographic(double a, double e2, double k0, 
		      Hemisphere hemi, double FN, double FE,
		      double N, double E,
		      double* lat_rad, double* lon_rad);

/* FORWARD AND BACK PROJECTIONS FOR AN ELLIPSOID ONTO THE UTM/UPS GRID */

int geographic_to_grid(double a, double e2,
		       double lat_rad, double lon_rad, 
		       GridZone* zone, Hemisphere* hemi, double* N, double* E);

int grid_to_geographic(double a, double e2,		       
		       GridZone zone, Hemisphere hemi, double N, double E,
		       double* lat_rad, double* lon_rad);

#endif /* UTM_H */
