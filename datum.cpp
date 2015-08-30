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
 * $Id: datum.cpp,v 1.4 2006/12/11 20:47:08 sfegan Exp $
 *
 *****************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "datum.h"

#define ECC2(f) (2.0*(1.0-1.0/(f))/(f)) /* Caution: Uses f twice! */
#define NUMOF(x) (sizeof(x)/sizeof(*x))

static StandardEllipse reference_ellipse[] = {
  { "Airy 1830",               "AA", 6377563.396, ECC2(299.3249646) },
  { "Australian National",     "AN", 6378160,     ECC2(298.25) },
  { "Bessel 1841, Ethiopia, Indonesia, Japan and Korea",
                               "BR", 6377397.155, ECC2(299.1528128) },
  { "Bessel 1841, Namibia",    "BN", 6377483.865, ECC2(299.1528128) },
  { "Clarke 1866",             "CC", 6378206.4,   ECC2(294.9786982) },
  { "Clarke 1880",             "CD", 6378249.145, ECC2(293.465) },
  { "Everest, Brunei and E. Malaysia (Sabah and Sarawak)",
                               "EB", 6377298.556, ECC2(300.8017) },
  { "Everest, India 1830",     "EA", 6377276.345, ECC2(300.8017) },
  { "Everest, India 1956",     "EC", 6377301.243, ECC2(300.8017) },
  { "Everest, Pakistan",       "EF", 6377309.613, ECC2(300.8017) },
  { "Everest, W. Malaysia and Singapore 1948", 
                               "EE", 6377304.063, ECC2(300.8017) },  
  { "Everest, W. Malaysia 1969",
                               "ED", 6377295.664, ECC2(300.8017) },
  { "Geodetic Reference System 1980",
                               "RF", 6378137,     ECC2(298.257222101) },
  { "Helmert 1906",            "HE", 6378200,     ECC2(298.3) },
  { "Hough 1960",              "HO", 6378270,     ECC2(297) },
  { "Indonesian 1974",         "ID", 6378160,     ECC2(298.247) },
  { "International 1924",      "IN", 6378388,     ECC2(297) },
  { "Krassovsky 1940",         "KA", 6378245,     ECC2(298.3) },
  { "Modified Airy",           "AM", 6377340.189, ECC2(299.3249646) },
  { "Modified Fischer 1960",   "FA", 6378155,     ECC2(298.3) },
  { "South American 1969",     "SA", 6378160,     ECC2(298.25) },
  { "WGS 1972",                "WD", 6378135,     ECC2(298.26) },
  { "WGS 1984",                "WE", 6378137,     ECC2(298.257223563) }
};

static DatumTransformationParameters datum_transformation[] = {

  /* APPENDIX B */
  /* LOCAL GEODETIC DATUMS RELATED TO WGS84 THROUGH SATELLITE TIES */
  
  /* Continent: AFRICA */
  
  { "ADINDAN, Mean Solution (Ethiopia and Sudan)", "ADI-M",
                                                ELLIPSE_CD, -166,  -15,  204 },
  { "ADINDAN, Burkina Faso",           "ADI-E", ELLIPSE_CD, -118,  -14,  218 },
  { "ADINDAN, Cameroon",               "ADI-F", ELLIPSE_CD, -134,   -2,  210 },
  { "ADINDAN, Ethiopia",               "ADI-A", ELLIPSE_CD, -165,  -11,  206 },
  { "ADINDAN, Mali",                   "ADI-C", ELLIPSE_CD, -123,  -20,  220 },
  { "ADINDAN, Senegal",                "ADI-D", ELLIPSE_CD, -128,  -18,  224 },
  { "ADINDAN, Sudan",                  "ADI-B", ELLIPSE_CD, -161,  -14,  205 },
  { "AFGOOYE, Somalia",                "AFG",   ELLIPSE_KA,  -43, -163,   45 },
  { "ARC 1950, Mean Solution (Botswana, Lesotho,Malawi, Swaziland, Zaire, "
    "Zambia and Zimbabwe)",            "ARF-M", ELLIPSE_CD, -143,  -90, -294 },

  { "ARC 1950, Botswana",              "ARF-A", ELLIPSE_CD, -138, -105, -289 },
  { "ARC 1950, Burundi",               "ARF-H", ELLIPSE_CD, -153,   -5, -292 },
  { "ARC 1950, Lesotho",               "ARF-B", ELLIPSE_CD, -125, -108, -295 },
  { "ARC 1950, Malawi",                "ARF-C", ELLIPSE_CD, -161,  -73, -317 },
  { "ARC 1950, Swaziland",             "ARF-D", ELLIPSE_CD, -134, -105, -295 },
  { "ARC 1950, Zaire",                 "ARF-E", ELLIPSE_CD, -169,  -19, -278 },
  { "ARC 1950, Zambia",                "ARF-F", ELLIPSE_CD, -147,  -74, -283 },
  { "ARC 1950, Zimbabwe",              "ARF-G", ELLIPSE_CD, -142,  -96, -293 },

  { "ARC 1960, Mean Solution (Kenya and Tanzania)",
                                       "ARS-M", ELLIPSE_CD, -160,   -6, -302 },
  { "ARC 1960, Kenya",                 "ARS-A", ELLIPSE_CD, -157,   -2, -299 },
  { "ARC 1960, Tanzania",              "ARS-B", ELLIPSE_CD, -175,  -23, -303 },
  { "AYABELLE LIGHTHOUSE, Djibouti",   "PHA",   ELLIPSE_CD,  -79, -129,  145 },
  { "BISSAU, Guinea-Bissau",           "BID",   ELLIPSE_IN, -173,  253,   27 },
  { "CAPE, South Africa",              "CAP",   ELLIPSE_CD, -136, -108,  292 },

  { "CARTHAGE, Tunisia",               "CGE",   ELLIPSE_CD, -263,    6,  431 },
  { "DABOLA, Guinea",                  "DAL",   ELLIPSE_CD,  -83,   37,  124 },
  { "EUROPEAN 1950, Egypt",            "EUR-F", ELLIPSE_IN, -130, -117, -151 },
  { "EUROPEAN 1950, Tunisia",          "EUR-T", ELLIPSE_IN, -112,  -77, -145 },
  { "LEIGON, Ghana",                   "LEH",   ELLIPSE_CD, -130,   29,  364 },
  { "LIBERIA 1964, Liberia",           "LIB",   ELLIPSE_CD,  -90,   40,   88 },

  { "MASSAWA, Eritrea (Ethiopia)",     "MAS",   ELLIPSE_BR,  639,  405,   60 },
  { "MERCHICH, Morocco",               "MER",   ELLIPSE_CD,   31,  146,   47 },
  { "MINNA, Cameroon",                 "MIN-A", ELLIPSE_CD,  -81,  -84,  115 },
  { "MINNA, Nigeria",                  "MIN-B", ELLIPSE_CD,  -92,  -93,  122 },
  { "M'PORALOKO, Gabon",               "MPO",   ELLIPSE_CD,  -74, -130,   42 },
  { "NORTH SAHARA 1959, Algeria",      "NSD",   ELLIPSE_CD, -186,  -93,  310 },

  { "OLD EGYPTIAN 1907, Egypt",        "OEG",   ELLIPSE_HE, -130,  110,  -13 },
  { "POINT 58, Mean Solution (Burkina Faso and Niger)",
                                       "PTB",   ELLIPSE_CD, -106, -129,  165 },
  { "POINTE NOIRE 1948, Congo",        "PTN",   ELLIPSE_CD, -148,   51, -291 },
  { "SCHWARZECK, Namibia",             "SCK",   ELLIPSE_BN,  616,   97, -251 },
  { "SIERRA LEONE 1960, Sierra Leone", "SRL",   ELLIPSE_CD,  -88,    4,  101 },
  { "VOIROL 1960, Algeria",            "VOR",   ELLIPSE_CD, -123, -206,  219 },

  /* Continent: ASIA */
  
  { "AIN EL ABD 1970, Bahrain Island", "AIN-A", ELLIPSE_IN, -150, -250,   -1 },
  { "AIN EL ABD 1970, Saudi Arabia",   "AIN-B", ELLIPSE_IN, -143, -236,    7 },
  { "DJAKARTA (BATAVIA), Sumatra (Indonesia)",
                                       "BAT",   ELLIPSE_BR, -377,  681,  -50 },
  { "EUROPEAN 1950, Iran",             "EUR-H", ELLIPSE_IN, -117, -132, -164 },
  { "HONG KONG 1963, Hong Kong",       "HKD",   ELLIPSE_IN, -156, -271, -189 },
  { "HU-TZU-SHAN, Taiwan",             "HTN",   ELLIPSE_IN, -637, -549, -203 },

  { "INDIAN, Bangladesh",              "IND-B", ELLIPSE_EA,  282,  726,  254 },
  { "INDIAN, India and Nepal",         "IND-I", ELLIPSE_EC,  295,  736,  257 },
  { "INDIAN 1954, Thailand",           "INF-A", ELLIPSE_EA,  217,  823,  299 },
  { "INDIAN 1960, Vietnam (near 16°N)","ING-A", ELLIPSE_EA,  198,  881,  317 },
  { "INDIAN 1960, Con Son Island (Vietnam)",
                                       "ING-B", ELLIPSE_EA,  182,  915,  344 },
  { "INDIAN 1975, Thailand",           "INH-A", ELLIPSE_EA,  209,  818,  290 },
  { "INDIAN 1975, Thailand",           "INH-A1",ELLIPSE_EA,  210,  814,  289 },

  { "INDONESIAN 1974, Indonesia",      "IDN",   ELLIPSE_ID,  -24,  -15,    5 },
  { "KANDAWALA, Sri Lanka",            "KAN",   ELLIPSE_EA,  -97,  787,   86 },
  { "KERTAU 1948, West Malaysia and Singapore",  
                                       "KEA",   ELLIPSE_EE,  -11,  851,    5 },
  { "KOREAN GEODETIC SYSTEM 1995, South Korea",
                                       "KGS",   ELLIPSE_WE,    0,    0,    0 },

  { "NAHRWAN, Masirah Island (Oman)",  "NAH-A", ELLIPSE_CD, -247, -148,  369 },
  { "NAHRWAN, United Arab Emirates",   "NAH-B", ELLIPSE_CD, -249, -156,  381 },
  { "NAHRWAN, Saudi Arabia",           "NAH-C", ELLIPSE_CD, -243, -192,  477 },
  { "OMAN, Oman",                      "FAH",   ELLIPSE_CD, -346,   -1,  224 },
  { "QATAR NATIONAL, Qatar",           "QAT",   ELLIPSE_IN, -128, -283,   22 },
  { "SOUTH ASIA, Singapore",           "SOA",   ELLIPSE_FA,    7,  -10,  -26 },

  { "TIMBALAI 1948, Brunei and East Malaysia (Sarawak and Sabah)",
                                       "TIL",   ELLIPSE_EB, -679,  669,  -48 },
  { "TOKYO, Mean Solution (Japan, Okinawa and South Korea)",
                                       "TOY-M", ELLIPSE_BR, -148,  507,  685 },
  { "TOKYO, Japan",                    "TOY-A", ELLIPSE_BR, -148,  507,  685 },
  { "TOKYO, Okinawa",                  "TOY-C", ELLIPSE_BR, -158,  507,  676 },
  { "TOKYO, South Korea",              "TOY-B", ELLIPSE_BR, -146,  507,  687 },
  { "TOKYO, South Korea",              "TOY-B1",ELLIPSE_BR, -147,  506,  687 },

  /* Continent: AUSTRALIA */

  { "AUSTRALIAN GEODETIC 1966, Australia and Tasmania",
                                       "AUA",   ELLIPSE_AN, -133,  -48,  148 },
  { "AUSTRALIAN GEODETIC 1984, Australia and Tasmania",
                                       "AUG",   ELLIPSE_AN, -134,  -48,  149 },
  
  /* Continent: EUROPE */

  { "CO-ORDINATE SYSTEM 1937 OF ESTONIA, Estonia",
                                       "EST",   ELLIPSE_BR,  374,  150,  588 },
  { "EUROPEAN 1950, Mean Solution {Austria, Belgium, Denmark, Finland, "
    "France, FRG (Federal Republic of Germany), Gibraltar, Greece, Italy, "
    "Luxembourg, Netherlands, Norway, Portugal, Spain, Sweden and "
    "Switzerland}",                    "EUR-M", ELLIPSE_IN,  -87,  -98, -121 },

  { "EUROPEAN 1950, Western Europe {Limited to Austria, Denmark, France, FRG "
    "(Federal Republic of Germany), Netherlands and Switzerland}",
                                       "EUR-A", ELLIPSE_IN,  -87,  -96, -120 },
  { "EUROPEAN 1950, Cyprus",           "EUR-E", ELLIPSE_IN, -104, -101, -140 },
#if 0
  { "EUROPEAN 1950, Egypt",            "EUR-F", ELLIPSE_IN, -130, -117, -151 },
#endif
  { "EUROPEAN 1950, England, Channel Islands, Scotland and Shetland Islands",
                                       "EUR-G", ELLIPSE_IN,  -86,  -96, -120 },
  { "EUROPEAN 1950, England, Ireland, Scotland and Shetland Islands",
                                       "EUR-K", ELLIPSE_IN,  -86,  -96, -120 },

  { "EUROPEAN 1950, Greece",           "EUR-B", ELLIPSE_IN,  -84,  -95, -130 },
#if 0
  { "EUROPEAN 1950, Iran",             "EUR-H", ELLIPSE_IN, -117, -132, -164 },
#endif
  { "EUROPEAN 1950, Italy, Sardinia",  "EUR-I", ELLIPSE_IN,  -97, -103, -120 },
  { "EUROPEAN 1950, Italy, Sicily",    "EUR-J", ELLIPSE_IN,  -97,  -88, -135 },
  { "EUROPEAN 1950, Malta",            "EUR-L", ELLIPSE_IN, -107,  -88, -149 },
  { "EUROPEAN 1950, Norway and Finland",
                                       "EUR-C", ELLIPSE_IN,  -87,  -95, -120 },
  { "EUROPEAN 1950, Portugal and Spain",
                                       "EUR-D", ELLIPSE_IN,  -84, -107, -120 },
#if 0
  { "EUROPEAN 1950, Tunisia",          "EUR-T", ELLIPSE_IN, -112,  -77, -145 },
#endif

  { "EUROPEAN 1979, Mean Solution (Austria, Finland, Netherlands, Norway, "
    "Spain, Sweden and Switzerland)",  "EUS",   ELLIPSE_IN,  -86,  -98, -119 },
  { "HJORSEY 1955, Iceland",           "HJO",   ELLIPSE_IN,  -73,   46,  -86 },
  { "IRELAND 1965",                    "IRL",   ELLIPSE_AM,  506, -122,  611 },
  { "ORDNANCE SURVEY OF GREAT BRITAIN 1936, Mean Solution (England, "
    "Isle of Man, Scotland, Shetland Islands and Wales)",
                                       "OGB-M", ELLIPSE_AA,  375, -111,  431 },

  { "ORDNANCE SURVEY OF GREAT BRITAIN 1936, England",
                                       "OGB-A", ELLIPSE_AA,  371, -112,  434 },
  { "ORDNANCE SURVEY OF GREAT BRITAIN 1936, England, Isle of Man and Wales",
                                       "OGB-B", ELLIPSE_AA,  371, -111,  434 },
  { "ORDNANCE SURVEY OF GREAT BRITAIN 1936, Scotland and Shetland Islands",
                                       "OGB-C", ELLIPSE_AA,  384, -111,  425 },
  { "ORDNANCE SURVEY OF GREAT BRITAIN 1936, Wales",
                                       "OGB-D", ELLIPSE_AA,  370, -108,  434 },
  { "ROME 1940, Sardinia",             "MOD",   ELLIPSE_IN, -225,  -65,    9 },
  { "S-42 (PULKOVO 1942), Hungary",    "SPK-A", ELLIPSE_KA,   28, -121,  -77 },
  { "S-42 (PULKOVO 1942), Poland",     "SPK-B", ELLIPSE_KA,   23, -124,  -82 },
  
  { "S-42 (PULKOVO 1942), Czechoslovakia", 
                                       "SPK-C", ELLIPSE_KA,   26, -121,  -78 },
  { "S-42 (PULKOVO 1942), Latvia",     "SPK-D", ELLIPSE_KA,   24, -124,  -82 },
  { "S-42 (PULKOVO 1942), Kazakhstan", "SPK-E", ELLIPSE_KA,   15, -130,  -84 },
  { "S-42 (PULKOVO 1942), Albania",    "SPK-F", ELLIPSE_KA,   24, -130,  -92 },
  { "S-42 (PULKOVO 1942), Romania",    "SPK-G", ELLIPSE_KA,   28, -121,  -77 },

  { "S-JTSK Czechoslovakia",           "CCD",   ELLIPSE_BR,  589,   76,  480 },

  /* Continent: NORTH AMERICA */

  { "CAPE CANAVERAL, Mean Solution (Florida and Bahamas)",
                                       "CAC",   ELLIPSE_CC,   -2,  151,  181 },
  { "NORTH AMERICAN 1927, Mean Solution (CONUS)",
                                       "NAS-C", ELLIPSE_CC,   -8,  160,  176 },
  { "NORTH AMERICAN 1927, Western United States (Arizona, Arkansas, "
    "California, Colorado, Idaho, Iowa, Kansas, Montana, Nebraska, Nevada, "
    "New Mexico, North Dakota, Oklahoma, Oregon, South Dakota, Texas, Utah, "
    "Washington and Wyoming)",         "NAS-B", ELLIPSE_CC,   -8,  159,  175 },

  { "NORTH AMERICAN 1927, Eastern United States (Alabama, Connecticut, "
    "Delaware, District of Columbia, Florida, Georgia, Illinois, Indiana, "
    "Kentucky, Louisiana, Maine, Maryland, Massachusetts, Michigan, "
    "Minnesota, Mississippi, Missouri, New Hampshire, New Jersey, New York, "
    "North Carolina, Ohio, Pennsylvania, Rhode Island, South Carolina, "
    "Tennessee, Vermont, Virginia, West Virginia and Wisconsin)",
                                       "NAS-A", ELLIPSE_CC,   -9,  161,  179 },

  { "NORTH AMERICAN 1927, Alaska (Excluding Aleutian Islands)",
                                       "NAS-D", ELLIPSE_CC,   -5,  135,  172 },
  { "NORTH AMERICAN 1927, Aleutian Islands, East of 180°W",
                                       "NAS-V", ELLIPSE_CC,   -2,  152,  149 },
  { "NORTH AMERICAN 1927, Aleutian Islands, West of 180°W",
                                       "NAS-W", ELLIPSE_CC,    2,  204,  105 },
  { "NORTH AMERICAN 1927, Bahamas (Excluding San Salvador Island)",
                                       "NAS-Q", ELLIPSE_CC,   -4,  154,  178 },
  { "NORTH AMERICAN 1927, San Salvador Island",
                                       "NAS-R", ELLIPSE_CC,    1,  140,  165 },
  { "NORTH AMERICAN 1927, Canada Mean Solution (Including Newfoundland)",
                                       "NAS-E", ELLIPSE_CC,  -10,  158,  187 },
  { "NORTH AMERICAN 1927, Alberta and British Columbia",
                                       "NAS-F", ELLIPSE_CC,   -7,  162,  188 },

  { "NORTH AMERICAN 1927, Eastern Canada (Newfoundland, New Brunswick, "
    "Nova Scotia and Quebec)",         "NAS-G", ELLIPSE_CC,  -22,  160,  190 },
  { "NORTH AMERICAN 1927, Manitoba and Ontario",
                                       "NAS-H", ELLIPSE_CC,   -9,  157,  184 },
  { "NORTH AMERICAN 1927, Northwest Territories and Saskatchewan",
                                       "NAS-I", ELLIPSE_CC,    4,  159,  188 },
  { "NORTH AMERICAN 1927, Yukon",      "NAS-J", ELLIPSE_CC,   -7,  139,  181 },
  { "NORTH AMERICAN 1927, Canal Zone", "NAS-O", ELLIPSE_CC,    0,  125,  201 },
  { "NORTH AMERICAN 1927, Caribbean (Antigua Island, Barbados, Barbuda, "
    "Caicos Islands, Cuba, Dominican Republic, Grand Cayman, Jamaica and "
    "Turks Islands)",                  "NAS-P", ELLIPSE_CC,   -3,  142,  183 },
  
  { "NORTH AMERICAN 1927, Central America (Belize, Costa Rica, El Salvador, "
    "Guatemala, Honduras and Nicaragua)",
                                       "NAS-N", ELLIPSE_CC,    0,  125,  194 },
  { "NORTH AMERICAN 1927, Cuba",       "NAS-T", ELLIPSE_CC,   -9,  152,  178 },
  { "NORTH AMERICAN 1927, Greenland (Hayes Peninsula)",
                                       "NAS-U", ELLIPSE_CC,   11,  114,  195 },
  { "NORTH AMERICAN 1927, Mexico",     "NAS-L", ELLIPSE_CC,  -12,  130,  190 },
  { "NORTH AMERICAN 1983, Alaska (Excluding Aleutian Islands)",
                                       "NAR-A", ELLIPSE_RF,    0,    0,    0 },
  { "NORTH AMERICAN 1983, Aleutian Islands",
                                       "NAR-E", ELLIPSE_RF,   -2,    0,    4 },
  { "NORTH AMERICAN 1983, Canada",     "NAR-B", ELLIPSE_RF,    0,    0,    0 },

  { "NORTH AMERICAN 1983, CONUS",      "NAR-C", ELLIPSE_RF,    0,    0,    0 },
  { "NORTH AMERICAN 1983, Hawaii",     "NAR-H", ELLIPSE_RF,    1,    1,   -1 },
  { "NORTH AMERICAN 1983, Mexico and Central America",
                                       "NAR-D", ELLIPSE_RF,    0,    0,    0 },

  /* Continent: SOUTH AMERICA */

  { "BOGOTA OBSERVATORY, Colombia",    "BOO",   ELLIPSE_IN,  307,  304, -318 },
  { "CAMPO INCHAUSPE 1969, Argentina", "CAI",   ELLIPSE_IN, -148,  136,   90 },
  { "CHUA ASTRO, Paraguay",            "CHU",   ELLIPSE_IN, -134,  229,  -29 },
  { "CORREGO ALEGRE, Brazil",          "COA",   ELLIPSE_IN, -206,  172,   -6 },

  { "PROVISIONAL SOUTH AMERICAN 1956, Mean Solution (Bolivia, Chile, "
    "Colombia, Ecuador, Guyana, Peru and Venezuela)",
                                       "PRP-M", ELLIPSE_IN, -288,  175, -376 },
  { "PROVISIONAL SOUTH AMERICAN 1956, Bolivia",
                                       "PRP-A", ELLIPSE_IN, -270,  188, -388 },
  { "PROVISIONAL SOUTH AMERICAN 1956, Chile, Northern Chile (near 19°S)",
                                       "PRP-B", ELLIPSE_IN, -270,  183, -390 },
  { "PROVISIONAL SOUTH AMERICAN 1956, Southern Chile (near 43°S)",
                                       "PRP-C", ELLIPSE_IN, -305,  243, -442 },
  { "PROVISIONAL SOUTH AMERICAN 1956, Colombia",
                                       "PRP-D", ELLIPSE_IN, -282,  169, -371 },
  { "PROVISIONAL SOUTH AMERICAN 1956, Ecuador",
                                       "PRP-E", ELLIPSE_IN, -278,  171, -367 },

  { "PROVISIONAL SOUTH AMERICAN 1956, Guyana",
                                       "PRP-F", ELLIPSE_IN, -298,  159, -369 },
  { "PROVISIONAL SOUTH AMERICAN 1956, Peru",
                                       "PRP-G", ELLIPSE_IN, -279,  175, -379 },
  { "PROVISIONAL SOUTH AMERICAN 1956, Venezuela",
                                       "PRP-H", ELLIPSE_IN, -295,  173, -371 },
  { "PROVISIONAL SOUTH CHILEAN 1963, Southern Chile (near 53°S)",
                                       "HIT",   ELLIPSE_IN,   16,  196,   93 },

  { "SOUTH AMERICAN 1969, Mean Solution (Argentina, Bolivia, Brazil, Chile, "
    "Colombia, Ecuador, Guyana, Paraguay, Peru, Trinidad and Tobago "
    "and Venezuela)",                  "SAN-M", ELLIPSE_SA,  -57,    1,  -41 },
  { "SOUTH AMERICAN 1969, Argentina",  "SAN-A", ELLIPSE_SA,  -62,   -1,  -37 },
  { "SOUTH AMERICAN 1969, Bolivia",    "SAN-B", ELLIPSE_SA,  -61,    2,  -48 },
  { "SOUTH AMERICAN 1969, Brazil",     "SAN-C", ELLIPSE_SA,  -60,   -2,  -41 },
  { "SOUTH AMERICAN 1969, Chile",      "SAN-D", ELLIPSE_SA,  -75,   -1,  -44 },
  { "SOUTH AMERICAN 1969, Colombia",   "SAN-E", ELLIPSE_SA,  -44,    6,  -36 },
  
  { "SOUTH AMERICAN 1969, Ecuador (Excluding Galapagos Islands)",
                                       "SAN-F", ELLIPSE_SA,  -48,    3,  -44 },
  { "SOUTH AMERICAN 1969, Baltra and Galapagos Islands",
                                       "SAN-J", ELLIPSE_SA,  -47,   26,  -42 },
  { "SOUTH AMERICAN 1969, Guyana",     "SAN-G", ELLIPSE_SA,  -53,    3,  -47 },
  { "SOUTH AMERICAN 1969, Paraguay",   "SAN-H", ELLIPSE_SA,  -61,    2,  -33 },
  { "SOUTH AMERICAN 1969, Peru",       "SAN-I", ELLIPSE_SA,  -58,    0,  -44 },
  { "SOUTH AMERICAN 1969, Trinidad and Tobago", 
                                       "SAN-K", ELLIPSE_SA,  -45,   12,  -33 },
  { "SOUTH AMERICAN 1969, Venezuela",  "SAN-L", ELLIPSE_SA,  -45,    8,  -33 },

  { "SOUTH AMERICAN GEOCENTRIC REFERENCE SYSTEM (SIRGAS)",
                                       "SIR",   ELLIPSE_RF,    0,    0,    0 },
  { "ZANDERIJ, Suriname",              "ZAN",   ELLIPSE_IN, -265,  120, -358 },

  /* Continent: ATLANTIC OCEAN */
  
  { "ANTIGUA ISLAND ASTRO 1943, Antigua and Leeward Islands",
                                       "AIA",   ELLIPSE_CD, -270,   13,   62 },
  { "ASCENSION ISLAND 1958, Ascension Island",
                                       "ASC",   ELLIPSE_IN, -205,  107,   53 },
  { "ASTRO DOS 71/4, St. Helena Island",
                                       "SHB",   ELLIPSE_IN, -320,  550, -494 },
  { "BERMUDA 1957, Bermuda Islands",   "BER",   ELLIPSE_CC,  -73,  213,  296 },
#if 0
  { "CAPE CANAVERAL, Mean Solution (Bahamas and Florida)",
                                       "CAC",   ELLIPSE_CC,   -2,  151,  181 },
#endif

  { "DECEPTION ISLAND, Deception Island and Antarctica",
                                       "DID",   ELLIPSE_CD,  260,   12, -147 },
  { "FORT THOMAS 1955, Nevis, St. Kitts and Leeward Islands",
                                       "FOT",   ELLIPSE_CD,   -7,  215,  225 },
  { "GRACIOSA BASE SW 1948, Faial, Graciosa, Pico, Sao Jorge and Terceira"
    "Islands (Azores)",                "GRA",   ELLIPSE_IN, -104,  167,  -38 },
#if 0
  { "HJORSEY 1955, Iceland",           "HJO",   ELLIPSE_IN,  -73,   46,  -86 },
#endif
  { "ISTS 061 ASTRO 1968, South Georgia Island",
                                       "ISG",   ELLIPSE_IN, -794,  +25,  +25 },

  { "L. C. 5 ASTRO 1961, Cayman Brac Island",
                                       "LCF",   ELLIPSE_CC,   42,  124,  147 },
  { "MONTSERRAT ISLAND ASTRO 1958, Montserrat and Leeward Islands",
                                       "ASM",   ELLIPSE_CD,  174,  359,  365 },
  { "NAPARIMA BWI, Trinidad and Tobago", 
                                       "NAP",   ELLIPSE_IN,  -10,  375,  165 },
  { "OBSERVATORIO METEOROLOGICO 1939, Corvo and Flores Islands (Azores)",
                                       "FLO",   ELLIPSE_IN, -425, -169,   81 },
  { "PICO DE LAS NIEVES, Canary Islands", 
                                       "PLN",   ELLIPSE_IN, -307,  -92,  127 },

  { "PORTO SANTO, Porto Santo and Madeira Islands", 
                                       "POS",   ELLIPSE_IN, -499, -249,  314 },
  { "PUERTO RICO, Puerto Rico and Virgin Islands",
                                       "PUR",   ELLIPSE_CC,   11,   72, -101 },
  { "QORNOQ, South Greenland",         "QUO",   ELLIPSE_IN,  164,  138, -189 },
  { "SAO BRAZ, Sao Miguel and Santa Maria Islands (Azores)",
                                       "SAO",   ELLIPSE_IN, -203,  141,   53 },
  { "SAPPER HILL, East Falkland Island",
                                       "SAP",   ELLIPSE_IN, -355,   21,   72 },

  { "SELVAGEM GRANDE 1938, Salvage Islands",
                                       "SGM",   ELLIPSE_IN, -289, -124,   60 },
  { "TRISTAN ASTRO 1968, Tristan da Cunha",
                                       "TDC",   ELLIPSE_IN, -632,  438, -609 },

  /* Continent: INDIAN OCEAN */

  { "ANNA 1 ASTRO 1965, Cocos Islands","ANO",   ELLIPSE_AN, -491,  -22,  435 },
  { "GAN 1970, Republic of Maldives",  "GAA",   ELLIPSE_IN, -133, -321,   50 },
  { "ISTS 073 ASTRO 1969, Diego Garcia",
                                       "IST",   ELLIPSE_IN,  208, -435, -229 },
  { "KERGUELEN ISLAND 1949, Kerguelen Island",
                                       "KEG",   ELLIPSE_IN,  145, -187,  103 },
  { "MAHE 1971, Mahe Island",          "MIK",   ELLIPSE_CD,   41, -220, -134 },
  { "REUNION, Mascarene Islands",      "REU",   ELLIPSE_IN,   94, -948,-1262 },

  /* Continent: PACIFIC OCEAN */

  { "AMERICAN SAMOA 1962, American Samoa Islands",
                                       "AMA",   ELLIPSE_CC, -115,  118,  426 },
  { "ASTRO BEACON \"E\", Iwo Jima",    "ATF",   ELLIPSE_IN,  145,   75, -272 },
  { "ASTRO TERN ISLAND (FRIG) 1961, Tern Island",
                                       "TRN",   ELLIPSE_IN,  114, -116, -333 },
  { "ASTRONOMICAL STATION 1952, Marcus Island",
                                       "ASQ",   ELLIPSE_IN,  124, -234,  -25 },
  { "BELLEVUE (IGN),Efate and Erromango Islands",
                                       "IBE",   ELLIPSE_IN, -127, -769,  472 },

  { "CANTON ASTRO 1966, Phoenix Islands",
                                       "CAO",   ELLIPSE_IN,  298, -304, -375 },
  { "CHATHAM ISLAND ASTRO 1971, Chatham Island (New Zealand)",
                                       "CHI",   ELLIPSE_IN,  175,  -38,  113 },
  { "DOS 1968, Gizo Island (New Georgia Islands)",
                                       "GIZ",   ELLIPSE_IN,  230, -199, -752 },
  { "EASTER ISLAND 1967, Easter Island",
                                       "EAS",   ELLIPSE_IN,  211,  147,  111 },
  { "GEODETIC DATUM 1949, New Zealand","GEO",   ELLIPSE_IN,   84,  -22,  209 },
  { "GUAM 1963, Guam",                 "GUA",   ELLIPSE_CC, -100, -248,  259 },

  { "GUX l ASTRO, Guadalcanal Island", "DOB",   ELLIPSE_IN,  252, -209, -751 },
#if 0
  { "INDONESIAN 1974, Indonesia",      "IDN",   ELLIPSE_ID,  -24,  -15,    5 },
#endif
  { "JOHNSTON ISLAND 1961, Johnston Island",
                                       "JOH",   ELLIPSE_IN,  189,  -79, -202 },
  { "KUSAIE ASTRO 1951, Caroline Islands, Fed. States of Micronesia",
                                       "KUS",   ELLIPSE_IN,  647, 1777,-1124 },
  { "LUZON, Philippines (Excluding Mindanao Island)",
                                       "LUZ-A", ELLIPSE_CC, -133,  -77,  -51 },

  { "LUZON, Mindanao Island",          "LUZ-B", ELLIPSE_CC, -133,  -79,  -72 },
  { "MIDWAY ASTRO 1961, Midway Islands 2003",
                                       "MID",   ELLIPSE_IN,  403,  -81,  277 },
  { "MIDWAY ASTRO 1961, Midway Islands 1987",
                                       "MID-87",ELLIPSE_IN,  912,  -58, 1227 },
  { "OLD HAWAIIAN, Mean Solution",     "OHA-M", ELLIPSE_CC,   61, -285, -181 },
  { "OLD HAWAIIAN, Hawaii",            "OHA-A", ELLIPSE_CC,   89, -279, -183 },
  { "OLD HAWAIIAN, Kauai",             "OHA-B", ELLIPSE_CC,   45, -290, -172 },
  { "OLD HAWAIIAN, Maui",              "OHA-C", ELLIPSE_CC,   65, -290, -190 },
  { "OLD HAWAIIAN, Oahu",              "OHA-D", ELLIPSE_CC,   58, -283, -182 },
  { "OLD HAWAIIAN, Mean Solution",     "OHI-M", ELLIPSE_IN,  201, -228, -346 },
  { "OLD HAWAIIAN, Hawaii",            "OHI-A", ELLIPSE_IN,  229, -222, -348 },

  { "OLD HAWAIIAN, Kauai",             "OHI-B", ELLIPSE_IN,  185, -233, -337 },
  { "OLD HAWAIIAN, Maui",              "OHI-C", ELLIPSE_IN,  205, -233, -355 },
  { "OLD HAWAIIAN, Oahu",              "OHI-D", ELLIPSE_IN,  198, -226, -347 },
  { "PITCAIRN ASTRO 1967, Pitcairn Island", 
                                       "PIT",   ELLIPSE_IN,  185,  165,   42 },
  { "SANTO (DOS) 1965, Espirito Santo Island",
                                       "SAE",   ELLIPSE_IN,  170,   42,   84 },
  { "VITI LEVU 1916, Viti Levu Island (Fiji Islands)",
                                       "MVS",   ELLIPSE_CD,   51,  391,  -36 },
  { "WAKE-ENIWETOK 1960, Marshall Islands",
                                       "ENW",   ELLIPSE_HO,  102,   52,  -38 },

  { "WAKE ISLAND ASTRO 1952, Wake Atoll",
                                       "WAK",   ELLIPSE_IN,  276,  -57,  149 },

  /* APPENDIX C */
  /* LOCAL GEODETIC DATUMS RELATED TO WGS84 THROUGH NON-SATELLITE TIES */

  { "BUKIT RIMPAH, Bangka and Belitung Islands (Indonesia)",
                                       "BUR",   ELLIPSE_BR, -384,  664,  -48 },
  { "CAMP AREA ASTRO, Camp McMurdo Area, Antarctica",
                                       "CAZ",   ELLIPSE_IN, -104, -129,  239 },
  { "EUROPEAN 1950, Iraq, Israel, Jordan, Kuwait, Lebanon, Saudi Arabia and "
    "Syria",                           "EUR-S", ELLIPSE_IN, -103, -106, -141 },
  { "GUNUNG SEGARA, Kalimantan (Indonesia)",
                                       "GSE",   ELLIPSE_BR, -403,  684,   41 },
  { "HERAT NORTH, Afghanistan",        "HEN",   ELLIPSE_IN, -333, -222,  114 },
  
  { "HERMANNSKOGEL, Yugoslavia (Prior to 1990) Slovenia, Croatia, Bosnia "
    "and Herzegovina and Serbia",      "HER",   ELLIPSE_BR,  682, -203,  480 },
  { "INDIAN, Pakistan",                "IND-P", ELLIPSE_EF,  283,  682,  231 },
  { "PULKOVO 1942, Russia",            "PUK",   ELLIPSE_KA,   28, -130,  -95 },
  { "TANANARIVE OBSERVATORY 1925, Madagascar",
                                       "TAN",   ELLIPSE_IN, -189, -242,  -91 },
  { "VOIROL 1874, Tunisia and Algeria","VOI",   ELLIPSE_CD,  -73, -247,  227 },
  { "YACARE, Uruguay",                 "YAC",   ELLIPSE_IN, -155,  171,   37 },

};

const StandardEllipse* _precompiled_ellipse(EllipseID id)
{
  unsigned iid = (unsigned)id;
  if(iid>=NUMOF(reference_ellipse))iid=(unsigned)ELLIPSE_WGS84;
  return &reference_ellipse[iid];
}

Ellipse* standard_ellipse(EllipseID id)
{
  return copy_ellipse(_precompiled_ellipse(id));
}

Ellipse* copy_ellipse(const StandardEllipse* e)
{
  Ellipse* eout = (Ellipse*)malloc(sizeof(Ellipse));
  eout->name    = (char*)malloc(strlen(e->name)+1);
  eout->id_code = (char*)malloc(strlen(e->id_code)+1);
  eout->a       = e->a;
  eout->e2      = e->e2;
  strcpy(eout->name, e->name);
  strcpy(eout->id_code, e->id_code);
  return eout;
}

Ellipse* copy_ellipse(const Ellipse* e)
{
  Ellipse* eout = (Ellipse*)malloc(sizeof(Ellipse));
  eout->name    = (char*)malloc(strlen(e->name)+1);
  eout->id_code = (char*)malloc(strlen(e->id_code)+1);
  eout->a       = e->a;
  eout->e2      = e->e2;
  strcpy(eout->name, e->name);
  strcpy(eout->id_code, e->id_code);
  return eout;
}

void free_ellipse(Ellipse* e)
{
  free(e->name);
  free(e->id_code);
  free(e);
}

#ifdef ENUM_VERIFY_MAIN
#include<iostream>
#include<set>

int main()
{
  std::set<std::string> datums;
  for(unsigned idatum=0;idatum<NUMOF(datum_transformation);idatum++)
    {
      std::string id(datum_transformation[idatum].id_code);
      for(std::string::iterator ichar = id.begin(); ichar!=id.end(); ichar++)
	if(*ichar == '-')*ichar='_';
      id = std::string("DATUM_")+id;
      if(datums.find(id) != datums.end())
	std::cout << "    /* DUPLICATE " << id << " */" << std::endl;
      else
	std::cout << "    " << id << ',' << std::endl;
      datums.insert(id);
    }
}

#endif
