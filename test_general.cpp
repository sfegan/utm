// Sample code: test_general.cpp
// Compile: g++ -o test_general test_general.cpp datum.cpp utm.cpp

// Calculate latitude and longitude of some location, given its surveyed
// position on the ground relative to a known reference point

#include<iostream>
#include<iomanip>
#include<cmath>
#include"datum.h"
#include"utm.h"

int main(int argc, char** argv)
{
  double lat_ref_rad = 32.0/180.0*M_PI; // Location of reference point
  double lon_ref_rad = -120/180.0*M_PI;

  double N = 150.22; // Location of surveyed point, relative to reference
  double E = 300.50; 

  const Ellipse* e = standard_ellipse(ELLIPSE_WGS84);

  // Calculate projection of reference point using TM from meridian
  // going through the reference point, with k0=1.0, so that scale is
  // accurately represented on the projection meridian. E_ref should
  // be zero.

  double N_ref;
  double E_ref;

  geographic_to_tm(e->a, e->e2, 1.0, lon_ref_rad, 0, 0, 
		   lat_ref_rad, lon_ref_rad,
		   &N_ref, &E_ref);

  std::cout << std::fixed << std::setprecision(2)
	    << E_ref << ' ' << N_ref << std::endl;

  // Calculate the position of the surveyed point on the ellipsoid
  // using inverse TM projection from the meridian through the
  // reference meridian with the scale factor set to 1.0. The position
  // of the reference point is negated and used as a false northing
  // and easting, thereby setting the origin of the E,N to be the
  // reference point.

  double lat_rad;
  double lon_rad;

  tm_to_geographic(e->a, e->e2, 1.0, lon_ref_rad, -N_ref, -E_ref,
		   N, E,
		   &lat_rad, &lon_rad);

  std::cout << std::fixed << std::setprecision(6)
            << lat_rad/M_PI*180.0 << ' ' << lon_rad/M_PI*180.0 << std::endl;
}
