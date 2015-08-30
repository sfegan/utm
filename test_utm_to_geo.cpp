// Sample code: test_utm_to_geo.cpp
// Compile: g++ -o test_utm_to_geo test_utm_to_geo.cpp datum.cpp utm.cpp

#include<iostream>
#include<iomanip>
#include<cmath>
#include"datum.h"
#include"utm.h"

int main(int argc, char** argv)
{
  double N = 3544404.13;
  double E = 216577.22;
  GridZone zone = UTM_ZONE_11;
  Hemisphere hemi = HEMI_NORTH;
  double lat;
  double lon;
  const Ellipse* e = ellipse(ELLIPSE_WGS84);
  grid_to_geographic(e->a, e->e2, zone, hemi, N, E, &lat, &lon);
  std::cout << std::fixed << std::setprecision(6)
            << lat/M_PI*180.0 << ' ' << lon/M_PI*180.0 << ' ' << std::endl;
}
