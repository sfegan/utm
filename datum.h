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
 * $Id: datum.h,v 1.4 2006/12/01 18:16:08 sfegan Exp $
 *
 *****************************************************************************/

#ifndef DATUM_H
#define DATUM_H

enum EllipseID
  {
    ELLIPSE_AA = 0,
    ELLIPSE_AN,    ELLIPSE_BR,    ELLIPSE_BN,    ELLIPSE_CC,    ELLIPSE_CD,
    ELLIPSE_EB,    ELLIPSE_EA,    ELLIPSE_EC,    ELLIPSE_EF,    ELLIPSE_EE,
    ELLIPSE_ED,    ELLIPSE_RF,    ELLIPSE_HE,    ELLIPSE_HO,    ELLIPSE_ID,
    ELLIPSE_IN,    ELLIPSE_KA,    ELLIPSE_AM,    ELLIPSE_FA,    ELLIPSE_SA,
    ELLIPSE_WD,    ELLIPSE_WE
  };


#define ELLIPSE_AUSTRALIAN    ELLIPSE_AN
#define ELLIPSE_BESSEL        ELLIPSE_BR
#define ELLIPSE_CLARKE_1866   ELLIPSE_CC
#define ELLIPSE_CLARKE_1880   ELLIPSE_CD
#define ELLIPSE_GRS80         ELLIPSE_EF
#define ELLIPSE_INT24         ELLIPSE_IN
#define ELLIPSE_WGS72         ELLIPSE_WD
#define ELLIPSE_WGS84         ELLIPSE_WE

struct StandardEllipse
{
  const char* name;
  const char* id_code;
  double a;
  double e2;
};

struct Ellipse
{
  char* name;
  char* id_code;
  double a;
  double e2;
};

#if !defined(__cplusplus)
typedef enum GridZone GridZone;
typedef enum Hemisphere Hemisphere;
typedef enum EllipseName EllipseName;
typedef enum EllipseID EllipseID;
typedef struct Ellipse Ellipse;
#endif

enum DatumID
  {
    DATUM_ADI_M = 0,
    DATUM_ADI_E,  DATUM_ADI_F,  DATUM_ADI_A,  DATUM_ADI_C,  DATUM_ADI_D,
    DATUM_ADI_B,  DATUM_AFG,    DATUM_ARF_M,  DATUM_ARF_A,  DATUM_ARF_H,
    DATUM_ARF_B,  DATUM_ARF_C,  DATUM_ARF_D,  DATUM_ARF_E,  DATUM_ARF_F,
    DATUM_ARF_G,  DATUM_ARS_M,  DATUM_ARS_A,  DATUM_ARS_B,  DATUM_PHA,
    DATUM_BID,    DATUM_CAP,    DATUM_CGE,    DATUM_DAL,    DATUM_EUR_F,
    DATUM_EUR_T,  DATUM_LEH,    DATUM_LIB,    DATUM_MAS,    DATUM_MER,
    DATUM_MIN_A,  DATUM_MIN_B,  DATUM_MPO,    DATUM_NSD,    DATUM_OEG,
    DATUM_PTB,    DATUM_PTN,    DATUM_SCK,    DATUM_SRL,    DATUM_VOR,
    DATUM_AIN_A,  DATUM_AIN_B,  DATUM_BAT,    DATUM_EUR_H,  DATUM_HKD,
    DATUM_HTN,    DATUM_IND_B,  DATUM_IND_I,  DATUM_INF_A,  DATUM_ING_A,
    DATUM_ING_B,  DATUM_INH_A,  DATUM_INH_A1, DATUM_IDN,    DATUM_KAN,
    DATUM_KEA,    DATUM_KGS,    DATUM_NAH_A,  DATUM_NAH_B,  DATUM_NAH_C,
    DATUM_FAH,    DATUM_QAT,    DATUM_SOA,    DATUM_TIL,    DATUM_TOY_M,
    DATUM_TOY_A,  DATUM_TOY_C,  DATUM_TOY_B,  DATUM_TOY_B1, DATUM_AUA,
    DATUM_AUG,    DATUM_EST,    DATUM_EUR_M,  DATUM_EUR_A,  DATUM_EUR_E,
    DATUM_EUR_G,  DATUM_EUR_K,  DATUM_EUR_B,  DATUM_EUR_I,  DATUM_EUR_J,
    DATUM_EUR_L,  DATUM_EUR_C,  DATUM_EUR_D,  DATUM_EUS,    DATUM_HJO,
    DATUM_IRL,    DATUM_OGB_M,  DATUM_OGB_A,  DATUM_OGB_B,  DATUM_OGB_C,
    DATUM_OGB_D,  DATUM_MOD,    DATUM_SPK_A,  DATUM_SPK_B,  DATUM_SPK_C,
    DATUM_SPK_D,  DATUM_SPK_E,  DATUM_SPK_F,  DATUM_SPK_G,  DATUM_CCD,
    DATUM_CAC,    DATUM_NAS_C,  DATUM_NAS_B,  DATUM_NAS_A,  DATUM_NAS_D,
    DATUM_NAS_V,  DATUM_NAS_W,  DATUM_NAS_Q,  DATUM_NAS_R,  DATUM_NAS_E,
    DATUM_NAS_F,  DATUM_NAS_G,  DATUM_NAS_H,  DATUM_NAS_I,  DATUM_NAS_J,
    DATUM_NAS_O,  DATUM_NAS_P,  DATUM_NAS_N,  DATUM_NAS_T,  DATUM_NAS_U,
    DATUM_NAS_L,  DATUM_NAR_A,  DATUM_NAR_E,  DATUM_NAR_B,  DATUM_NAR_C,
    DATUM_NAR_H,  DATUM_NAR_D,  DATUM_BOO,    DATUM_CAI,    DATUM_CHU,
    DATUM_COA,    DATUM_PRP_M,  DATUM_PRP_A,  DATUM_PRP_B,  DATUM_PRP_C,
    DATUM_PRP_D,  DATUM_PRP_E,  DATUM_PRP_F,  DATUM_PRP_G,  DATUM_PRP_H,
    DATUM_HIT,    DATUM_SAN_M,  DATUM_SAN_A,  DATUM_SAN_B,  DATUM_SAN_C,
    DATUM_SAN_D,  DATUM_SAN_E,  DATUM_SAN_F,  DATUM_SAN_J,  DATUM_SAN_G,
    DATUM_SAN_H,  DATUM_SAN_I,  DATUM_SAN_K,  DATUM_SAN_L,  DATUM_SIR,
    DATUM_ZAN,    DATUM_AIA,    DATUM_ASC,    DATUM_SHB,    DATUM_BER,
    DATUM_DID,    DATUM_FOT,    DATUM_GRA,    DATUM_ISG,    DATUM_LCF,
    DATUM_ASM,    DATUM_NAP,    DATUM_FLO,    DATUM_PLN,    DATUM_POS,
    DATUM_PUR,    DATUM_QUO,    DATUM_SAO,    DATUM_SAP,    DATUM_SGM,
    DATUM_TDC,    DATUM_ANO,    DATUM_GAA,    DATUM_IST,    DATUM_KEG,
    DATUM_MIK,    DATUM_REU,    DATUM_AMA,    DATUM_ATF,    DATUM_TRN,
    DATUM_ASQ,    DATUM_IBE,    DATUM_CAO,    DATUM_CHI,    DATUM_GIZ,
    DATUM_EAS,    DATUM_GEO,    DATUM_GUA,    DATUM_DOB,    DATUM_JOH,
    DATUM_KUS,    DATUM_LUZ_A,  DATUM_LUZ_B,  DATUM_MID,    DATUM_MID_87,
    DATUM_OHA_M,  DATUM_OHA_A,  DATUM_OHA_B,  DATUM_OHA_C,  DATUM_OHA_D,
    DATUM_OHI_M,  DATUM_OHI_A,  DATUM_OHI_B,  DATUM_OHI_C,  DATUM_OHI_D,
    DATUM_PIT,    DATUM_SAE,    DATUM_MVS,    DATUM_ENW,    DATUM_WAK,
    DATUM_BUR,    DATUM_CAZ,    DATUM_EUR_S,  DATUM_GSE,    DATUM_HEN,
    DATUM_HER,    DATUM_IND_P,  DATUM_PUK,    DATUM_TAN,    DATUM_VOI,
    DATUM_YAC,

    DATUM_WGS84 = 1000,
    DATUM_WGS72,
  };

#define DATUM_NAD27   DATUM_NAS_C

struct DatumTransformationParameters
{
  const char* name;
  const char* id_code;
  EllipseID ellipse;
  double delta_x;
  double delta_y;
  double delta_z;
};

#if !defined(__cplusplus)
typedef struct DatumTransformationParameters DatumTransformationParameters;
#endif

const StandardEllipse* _precompiled_ellipse(EllipseID id);
Ellipse* standard_ellipse(EllipseID id);
Ellipse* copy_ellipse(const Ellipse* e);
Ellipse* copy_ellipse(const StandardEllipse* e);
void free_ellipse(Ellipse* e);

#endif // DATUM_H
