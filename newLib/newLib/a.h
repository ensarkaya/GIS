#pragma once
#include <iostream>
#include "gdal.h"
#include "gdal_priv.h"
#include "gdal_utils.h"
#include "alg\gdal_alg.h"
#include "commonutils.h"
#include "gdal_utils_priv.h"
#include "ogrsf_frmts\ogrsf_frmts.h"
#include "cpl_string.h"
#include <Math.h>

using namespace std;


extern "C" double interpolation(double q11, double q12, double q21, double q22, int C, int R, double c, double r);
extern "C" GDALDataset* gdalCreateDataset(const char* path);
extern "C" void getGEoTransform(GDALDataset *dataset, unsigned int &width, unsigned int & height, double & leftUpperX, double &leftUpperY, double &pixelSizeX, double &pixelSizeY);
extern "C" double readElevationVal(const char input[], double X, double Y);
extern "C" void readElevationValInner(GDALDataset * dataset, double X, double Y, double &elevation);
extern "C" void measureDistance(double sLat, double sLong, double eLat, double eLong, double & distance);
extern "C" void Sphtan(double a, double b, double c, double & angab, double & angac);
extern "C" void Sepang(double ph1, double th1, double dist, double ph2, double th2, double & an1, double & an2);
extern "C" void FindAzimuth(double lat1, double lon1, double lat2, double lon2, double & an1, double &  an2);
extern "C" void FindNewCoord(double longitude, double latitude, double azimuth, double distance, double &  rLongitude, double &  rLatitude);
extern "C" void pathProfile(const char input[], double lat1, double lon1, double lat2, double lon2, double arr[]);
extern "C" int size(double sLat, double sLong, double eLat, double eLong);
extern "C" void los(double arr[],int size, double height1, double heigth2,double ret[]);
extern "C" void lineCheck(const char input[], double lat1, double lon1, double lat2, double lon2, double height1, double height2, double arr2[], double coorLat[], double coorLon[]);
extern "C" int shapeFile(const char input[], double lat1, double lon1, double height1, double height2, double radius, double** arr2d, double** coor2d);
extern "C" void createSHP(double ** arr2d, double ** coor, double radius, string output);
extern "C" static void CreateElevAttrib(const char* pszElevAttrib, OGRLayerH hLayer);
extern "C" int contour(const char * pszSrcFilename1, const char *pszDstFilename1, double dfInterval1);

struct GDALRasterizeOptions
{
	/*! output format. Use the short format name. */
	char *pszFormat;

	/*! the progress function to use */
	GDALProgressFunc pfnProgress;

	/*! pointer to the progress data variable */
	void *pProgressData;

	bool bCreateOutput;
	bool b3D;
	bool bInverse;
	char **papszLayers;
	char *pszSQL;
	char *pszDialect;
	char *pszBurnAttribute;
	char *pszWHERE;
	std::vector<int> anBandList;
	std::vector<double> adfBurnValues;
	char **papszRasterizeOptions;
	char **papszTO;
	double dfXRes;
	double dfYRes;
	char **papszCreationOptions;
	GDALDataType eOutputType;
	std::vector<double> adfInitVals;
	int bNoDataSet;
	double dfNoData;
	OGREnvelope sEnvelop;
	bool bGotBounds;
	int nXSize, nYSize;
	OGRSpatialReferenceH hSRS;
	bool bTargetAlignedPixels;
};
struct GDALGridOptions
{
	/*! output format. Use the short format name. */
	char *pszFormat;

	/*! allow or suppress progress monitor and other non-error output */
	bool bQuiet;

	/*! the progress function to use */
	GDALProgressFunc pfnProgress;

	/*! pointer to the progress data variable */
	void *pProgressData;

	char            **papszLayers;
	char            *pszBurnAttribute;
	double          dfIncreaseBurnValue;
	double          dfMultiplyBurnValue;
	char            *pszWHERE;
	char            *pszSQL;
	GDALDataType    eOutputType;
	char            **papszCreateOptions;
	int             nXSize;
	int             nYSize;
	double          dfXMin;
	double          dfXMax;
	double          dfYMin;
	double          dfYMax;
	bool            bIsXExtentSet;
	bool            bIsYExtentSet;
	GDALGridAlgorithm eAlgorithm;
	void            *pOptions;
	char            *pszOutputSRS;
	OGRGeometry     *poSpatialFilter;
	bool            bClipSrc;
	OGRGeometry     *poClipSrc;
	char            *pszClipSrcDS;
	char            *pszClipSrcSQL;
	char            *pszClipSrcLayer;
	char            *pszClipSrcWhere;
	bool             bNoDataSet;
	double           dfNoDataValue;
};