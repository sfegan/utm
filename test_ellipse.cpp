// Sample code: test_ellipse.cpp
// Compile: g++ -o test_ellipse test_ellipse.cpp datum.cpp utm.cpp

// Reporoduce tables 2-11 AND 3-7 of DMA TM 8358.2

#include<sstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include"datum.h"
#include"utm.h"

bool dmsStringToRad(const std::string& str, double& rad)
{
  unsigned degs=0;
  unsigned mins=0;
  unsigned secs=0;
  unsigned fracs=0;
  unsigned frac10s=1;
  unsigned i=0;

  bool negative=false;
  if(str[i]=='-') { negative = true; i++; }
  else if(str[i]=='+') { negative = false; i++; }

  while(i<str.length())
    {
      if((str[i]>='0')&&(str[i]<='9'))degs=degs*10+(str[i++]-'0');
      else if((str[i]==':')||(str[i]=='d')) { i++; break; }
      else return false;
    }

  while(i<str.length())
    {
      if((str[i]>='0')&&(str[i]<='9'))mins=mins*10+(str[i++]-'0');
      else if((str[i]==':')||(str[i]=='m')) { i++; break; }
      else return false;
    }

  while(i<str.length())
    {
      if((str[i]>='0')&&(str[i]<='9'))secs=secs*10+(str[i++]-'0');
      else if(str[i]=='.') { i++; break; }
      else if(str[i]=='s') { break; }
      else return false;
    }

  while(i<str.length())
    {
      if((str[i]>='0')&&(str[i]<='9'))
        { fracs=fracs*10+(str[i++]-'0'); frac10s=frac10s*10; }
      else if(str[i]=='s') { i++; break; }
      else return false;
    }

  if(i<str.length())return false;

  rad = ((negative?-1:1)*
	 (double(degs)+double(mins)/60+double(secs)/(60*60)+
	  double(fracs)/double(frac10s)/(60*60))/180.0*M_PI);
  
  return true;
}

std::string radToDMSString(double rad, 
			   unsigned sec_digits=3, bool dmsSep=false)
{
  const unsigned divisor10[] = { 1, 10, 100, 1000, 10000, 100000, 1000000,
				 10000000, 100000000, 1000000000 };

  double deg = fmod(fmod(rad/M_PI*180,360)+360,360);
  if(deg>=180)deg-=360;

  unsigned multiplier = divisor10[sec_digits];

  unsigned iangle = unsigned(floor(fabs(deg)*60*60*multiplier+0.5));

  bool negative =(deg<0);
  unsigned degs = iangle/(60*60*multiplier);
  unsigned mins = (iangle/(60*multiplier))%60;
  unsigned secs =(iangle/multiplier)%60;
  unsigned fsec = iangle%multiplier;

  std::ostringstream stream;
  stream << (negative?'-':'+')
         << std::setfill('0')
         << std::setw(3) << std::setprecision(3) << degs << (dmsSep?'d':':')
         << std::setw(2) << std::setprecision(2) << mins << (dmsSep?'m':':')
         << std::setw(2) << std::setprecision(2) << secs;
  if(sec_digits)stream << '.' << std::setw(sec_digits)
                       << std::setprecision(sec_digits) << fsec;
  if(dmsSep)stream << 's';

  return stream.str();
}

void write_entry(std::ostream& stream, bool fwd, 
		 double lat, double lon, GridZone zone, double N, double E)
{
  std::string zone_str;
  if(zone==UPS_NORTH)zone_str="NP";
  else if(zone==UPS_SOUTH)zone_str="SP";
  else zone_str=std::to_string((unsigned)zone);
  stream << "LAT: " << radToDMSString(lat) << "  LON: "
	 << radToDMSString(lon) << (fwd ? " --> " : " <-- ") << std::fixed
   << "Z: " << std::setw(2) << zone_str << "  "
	 << "N: " << std::setw(10) << std::setprecision(2) << N << "  " 
	 << "E: " << std::setw(9) << std::setprecision(2) << E 
	 << std::endl;
}

void write_entry(std::ostream& stream, bool fwd, 
		 double lat, double lon, GridZone zone, double N, double E, double gc_rad, double scale)
{
  std::string zone_str;
  if(zone==UPS_NORTH)zone_str="NP";
  else if(zone==UPS_SOUTH)zone_str="SP";
  else zone_str=std::to_string((unsigned)zone);
  stream << "LAT: " << radToDMSString(lat) << "  LON: "
	 << radToDMSString(lon) << (fwd ? " --> " : " <-- ") << std::fixed
   << "Z: " << std::setw(2) << zone_str << "  "
	 << "N: " << std::setw(10) << std::setprecision(2) << N << "  " 
	 << "E: " << std::setw(9) << std::setprecision(2) << E << "  "
   << "C: " << radToDMSString(gc_rad) << "  "
   << "S: " << std::setw(10) << std::setprecision(8) << scale
	 << std::endl;
}


// Test of the series approximate ellipsoidal TM conversion. Reproduce
// Section 2-11 from "The Universal Grids", Defense Mapping Agency
// Technical Manual (DMATM) 8358.2

int main()
{
  double a = 6378388.0;
  double e2 = 0.006722670022;

  double lon_rad;
  double lat_rad;
  GridZone zone;
  Hemisphere hemi;

  double E;
  double N;
  double gc_rad;
  double scale;

  std::cout << "Tests transformations to/from UTM grid (reproduces Table 2-11 of DMA TM 8358.2)\n\n";

  // ----------------------------------------------
  // TEST OF FORWARD GOING UTM ELLIPSOID CONVERSION
  // ----------------------------------------------

  // ------
  // ID = 1
  // ------

  dmsStringToRad("+045d00m00.000s",lon_rad);
  dmsStringToRad("+73d00m00.000s",lat_rad);
  zone = UTM_ZONE_38;
  hemi = HEMI_NORTH;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  // ------
  // ID = 2
  // ------

  dmsStringToRad("+102d00m00.000s",lon_rad);
  dmsStringToRad("+30d00m00.000s",lat_rad);
  zone = UTM_ZONE_47;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  zone = UTM_ZONE_48;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  // ------
  // ID = 3
  // ------

  dmsStringToRad("-113d54m43.321s",lon_rad);
  dmsStringToRad("+72d04m32.110",lat_rad);
  zone = UTM_ZONE_12;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  zone = UTM_ZONE_11;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  std::cout << std::endl;

  // -----------------------------------------------
  // TEST OF BACKWARD GOING UTM ELLIPSOID CONVERSION
  // -----------------------------------------------

  // ------
  // ID = 4
  // ------

  N = 3322824.35;
  E = 210577.93;

  zone = UTM_ZONE_48;
  grid_to_geographic(a, e2, zone, HEMI_NORTH, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);


  N = 3322824.08;
  E = 789411.59;

  zone = UTM_ZONE_47;
  grid_to_geographic(a, e2, zone, HEMI_NORTH, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  // ------
  // ID = 5
  // ------

  N = 1000000.00;
  E = 200000.00;

  zone = UTM_ZONE_31;
  grid_to_geographic(a, e2, zone, HEMI_NORTH, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  N = 1000491.75;
  E = 859739.88;

  zone = UTM_ZONE_30;
  grid_to_geographic(a, e2, zone, HEMI_NORTH, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  // ------
  // ID = 6
  // ------

  N = 9000000.00;
  E = 500000.00;

  zone = UTM_ZONE_43;
  grid_to_geographic(a, e2, zone, HEMI_NORTH, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  // ------
  // ID = 7
  // ------

  N = 4000000.00;
  E = 700000.00;

  zone = UTM_ZONE_30;
  grid_to_geographic(a, e2, zone, HEMI_SOUTH, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  N = 4000329.42;
  E = 307758.89;

  zone = UTM_ZONE_31;
  grid_to_geographic(a, e2, zone, HEMI_SOUTH, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  std::cout << std::endl;
  
  // ----------------------------------------------
  // TEST OF FORWARD GOING UPS ELLIPSOID CONVERSION
  // ----------------------------------------------

  std::cout << "Tests transformations to/from UPS grid (reproduces Table 3-7 of DMA TM 8358.2)\n\n";

  a = 6378137.0;
  e2 = 0.006694379990;

  // ------
  // ID = 1
  // ------

  dmsStringToRad("-132d14m52.761s",lon_rad);
  dmsStringToRad("+84d17m14.042s",lat_rad);
  zone = UPS_NORTH;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  // ------
  // ID = 2
  // ------

  dmsStringToRad("+044d00m00.000s",lon_rad);
  dmsStringToRad("+73d00m00.000s",lat_rad);
  zone = UPS_NORTH;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  // ------
  // ID = 3
  // ------

  dmsStringToRad("+132d14m52.303s",lon_rad);
  dmsStringToRad("-87d17m14.400s",lat_rad);
  zone = UPS_SOUTH;

  geographic_to_grid(a, e2, lat_rad, lon_rad, &zone, &hemi, &N, &E, &gc_rad, &scale);
  write_entry(std::cout, true, lat_rad, lon_rad, zone, N, E, gc_rad, scale);

  std::cout << std::endl;

  // -----------------------------------------------
  // TEST OF BACKWARD GOING UPS ELLIPSOID CONVERSION
  // -----------------------------------------------

  // ------
  // ID = 4
  // ------

  N = 2426773.60;
  E = 1530125.78;

  zone = UPS_NORTH;
  grid_to_geographic(a, e2, zone, HEMI_AUTO, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  // ------
  // ID = 5
  // ------

  N = 632668.43;
  E = 3320416.75;

  zone = UPS_NORTH;
  grid_to_geographic(a, e2, zone, HEMI_AUTO, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);

  // ------
  // ID = 6
  // ------

  N = 1500000.00;
  E = 2500000.00;

  zone = UPS_SOUTH;
  grid_to_geographic(a, e2, zone, HEMI_AUTO, N, E, &lat_rad, &lon_rad);
  write_entry(std::cout, false, lat_rad, lon_rad, zone, N, E);
}
