// Sample code: test_sphere.cpp
// Compile: g++ -o test_input test_input.cpp datum.cpp utm.cpp

#include<iostream>
#include<iomanip>
#include<cmath>
#include"datum.h"
#include"utm.h"

// Convert latitude and longitude in degrees to UTM grid coordinates

#include<iostream>
#include<sstream>
#include<iomanip>

#define RAD(x) ((x)/180.0*M_PI)

int main(int argc, char** argv)
{
  const char* program = *argv;
  argv++,argc--;

  std::istream* stream = &std::cin;
  if(argc)
    {
      stream = new std::ifstream(*argv);
      argv++,argc--;
    }

  const double a = 6378137.0;
  const double e2 = 0.006694379990;
      
  double lon;
  double lat;
  
  *stream >> lon >> lat;
  while(*stream)
    {
      lat = RAD(lat);
      lon = RAD(lon);

      double N;
      double E;
      
      GridZone zone = GRID_AUTO;
      Hemisphere hemi = HEMI_AUTO;

      geographic_to_grid(a, e2, lat, lon, &zone, &hemi, &N, &E);
      
      std::cout << std::fixed 
		    << std::setw(10) << std::setprecision(2) << N << "   "
		    << std::setw(10) << std::setprecision(2) << E << "   "
    		<< std::setw(2) << (unsigned)zone << "   "
    		<< std::setw(2) << (unsigned)hemi 
    		<< std::endl;

      *stream >> lon >> lat;
    }

  if(stream != &std::cin)delete stream;
}
