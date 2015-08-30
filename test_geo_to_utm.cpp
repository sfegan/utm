// Sample code: test_utm.cpp
// Compile: g++ -o test_geo_to_utm test_geo_to_utm.cpp datum.cpp utm.cpp

#include<iostream>
#include<iomanip>
#include<cmath>
#include"datum.h"
#include"utm.h"

int main(int argc, char** argv)
{
  double lat = 32.0/180.0*M_PI;
  double lon = -120/180.0*M_PI;
  double N;
  double E;
  GridZone zone = GRID_AUTO;
  Hemisphere hemi = HEMI_AUTO;
  const Ellipse* e = ellipse(ELLIPSE_WGS84);
  geographic_to_grid(e->a, e->e2, lat, lon, &zone, &hemi, &N, &E);
  std::cout << std::fixed << std::setprecision(2)
            << E << ' ' << N << ' '<< zone << ' ' << hemi << std::endl;
}
