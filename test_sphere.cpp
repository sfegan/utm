// Sample code: test_sphere.cpp
// Compile: g++ -o test_sphere test_sphere.cpp datum.cpp utm.cpp

// Reporoduce tables 2-11 AND 3-7 of DMA TM 8358.2

#include<iostream>
#include<iomanip>
#include<cmath>
#include"datum.h"
#include"utm.h"

// Test of the exact spherical TM conversion. Reproduce Table 10 from
// "MAP PROJECTIONS; A WORKING MANUAL", John Snyder, USGS Professional
// Paper 1395, 1987. Page 59-60.

#include<iostream>
#include<sstream>
#include<iomanip>

#define RAD(x) ((x)/180.0*M_PI)

int main()
{
  for(unsigned ilat=0;ilat<10;ilat++)
    {
      if(ilat!=0)std::cout << std::endl;

      double lat_rad=double(9-ilat)*RAD(10);
      double x[10];
      double y[10];
      for(unsigned ilon=0;ilon<10;ilon++)
	{
	  double lon_rad=double(ilon)*RAD(10);
	  geographic_to_tm_sphere(1.0, 1.0, 0.0, 0.0, 0.0, 
				  lat_rad, lon_rad, &y[ilon], &x[ilon]);

#if 0
	  // Test FORWARD followed by BACKWARD conversion
	  tm_to_geographic_sphere(1.0, 1.0, 0.0, 0.0, 0.0, 
				  y[ilon], x[ilon], &y[ilon], &x[ilon]);
	  y[ilon] *= 180/M_PI;
	  x[ilon] *= 180/M_PI;
#endif
	}

      for(unsigned ilon=0;ilon<10;ilon++)
	{
	  if(ilon)std::cout << ' ';
	  std::cout << std::fixed << std::setw(7) << std::setprecision(5)
		    << x[ilon];
	}
      std::cout << std::endl;
      for(unsigned ilon=0;ilon<10;ilon++)
	{
	  if(ilon)std::cout << ' ';
	  std::cout << std::fixed << std::setw(7) << std::setprecision(5)
		    << y[ilon];
	}
      std::cout << std::endl;
    }
}
