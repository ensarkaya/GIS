#include "a.h"
/// <summary>
/// interpolation function uses bilinear interpolation formula from 
///  "https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.1144-3-200111-S!!PDF-E.pdf" 
/// </summary>
/// <param name="q11"> height of q11 </param>
/// <param name="q12"> height of q12 </param>
/// <param name="q21">height of q21 </param>
/// <param name="q22">height of q22 </param>
/// <param name="C">constant- column number </param>
/// <param name="R">constant- row number </param>
/// <param name="c"> column of the point that we want to interpolate </param>
/// <param name="r"> row of the point that we want to interpolate</param>
double interpolation(double q11, double q12, double q21, double q22, int C, int R, double c, double r) {
	double y =
		q11 * ((R + 1 - r)*(C + 1 - c)) +
		q12 * ((R + 1 - r)*(c - C)) +
		q21 * ((r - R)*(C + 1 - c)) +
		q22 * ((r - R)*(c - C));
	return y;
}
/// <summary>
///  Creates a GDALDataset by using the path of the map
///  
/// </summary>
/// <param name="path"> Path of the input map </param>
GDALDataset* gdalCreateDataset(const char* path) {
	GDALAllRegister();

	GDALDataset  *dataset = (GDALDataset *)GDALOpen(path, GA_ReadOnly);

	if (dataset == NULL) {
		cout << "Error1 : The map cannot be found. Press any key to exit " <<endl;
		string s;
		cin >> s;
		exit(1);
	}
	return dataset;
}
/// <summary>
///  This funtion gets an GDALDataset and returns its necessary geotansform values 
/// </summary>
void getGEoTransform(GDALDataset *dataset, unsigned int &width, unsigned int & height, double & leftUpperX, double &leftUpperY, double &pixelSizeX, double &pixelSizeY)
{
	width = dataset->GetRasterXSize();
	height = dataset->GetRasterYSize();

	// Reading image coordinates
	double geoTransform[6];
	if (dataset->GetGeoTransform(geoTransform) == CE_None) {
		leftUpperX = geoTransform[0];
		leftUpperY = geoTransform[3];
		pixelSizeX = geoTransform[1];
		pixelSizeY = geoTransform[5];
	}
	else {
		cout << "Error2 : Failed read geotransform. Press any button to exit " << endl;
		string s;
		cin >> s;
		exit(2);
	}
}

/// <summary>
///  This function reads interpolated elevation value of a point given lat and lon values
/// </summary>
/// <param name="input">Path of the input map </param>
/// <param name="X"> Latitude(GEO-WGS84) </param>
/// <param name="Y"> Longitude(GEO-WGS84) </param>
double readElevationVal(const char input[], double X, double Y)
{
	double a=X, b=Y;
	double elevation = 0;
	GDALDataset * dataset = gdalCreateDataset(input);
	// Get image width and height values in tems of pixels
	unsigned int width, height;// This value cannot be lower than 0 so it's a unsigned int

	double leftUpperX, leftUpperY, pixelSizeX, pixelSizeY;
	getGEoTransform(dataset, width, height, leftUpperX, leftUpperY, pixelSizeX, pixelSizeY);

	// Getting the elevation raster band 
	GDALRasterBand  *elevationBand = dataset->GetRasterBand(1);

	double pixelX = 0.0, pixelY = 0.0, pixelX1 = 0.0, pixelY1 = 0.0;
	int intPixelX = 0, intPixelY = 0, intPixelX1 = 0, intPixelY1 = 0;

	pixelX = (a - (pixelSizeX / 2.0) - leftUpperX) / pixelSizeX;
	pixelY = (b + (pixelSizeX / 2.0) - leftUpperY) / pixelSizeY;
	intPixelX = (int)pixelX;
	intPixelY = (int)pixelY;

	//bilinear interpolation 
	int arr3[4] = { 0,0,0,0 };
	if (width >= intPixelX && height >= intPixelY)
	{
		elevationBand->RasterIO(GF_Read, intPixelX, intPixelY, 2, 2, &arr3[0], 2, 2, GDT_Int32, 0, 0);
		elevation = interpolation(arr3[0], arr3[1], arr3[2], arr3[3], intPixelX, intPixelY, pixelX, pixelY);
	}
	else {
		cout << "Error3 : Pixels out of range. Press any button to exit " << endl;
		string s;
		cin >> s;
		exit(3);
	}
	dataset->~GDALDataset();
	return elevation;
}
void readElevationValInner(GDALDataset * dataset, double X, double Y, double &elevation)
{
	// Get image width and height values in tems of pixels
	unsigned int width, height;// This value cannot be lower than 0 so it's a unsigned int

	double leftUpperX, leftUpperY, pixelSizeX, pixelSizeY;
	getGEoTransform(dataset, width, height, leftUpperX, leftUpperY, pixelSizeX, pixelSizeY);

	// Getting the elevation raster band 
	GDALRasterBand  *elevationBand = dataset->GetRasterBand(1);

	double pixelX = 0.0, pixelY = 0.0, pixelX1 = 0.0, pixelY1 = 0.0;
	int intPixelX = 0, intPixelY = 0, intPixelX1 = 0, intPixelY1 = 0;

	pixelX = (X - (pixelSizeX / 2.0) - leftUpperX) / pixelSizeX;
	pixelY = (Y + (pixelSizeX / 2.0) - leftUpperY) / pixelSizeY;
	intPixelX = (int)pixelX;
	intPixelY = (int)pixelY;

	//bilinear interpolation 
	int arr3[4] = { 0,0,0,0 };
	if (width >= intPixelX && height >= intPixelY)
	{
		elevationBand->RasterIO(GF_Read, intPixelX, intPixelY, 2, 2, &arr3[0], 2, 2, GDT_Int32, 0, 0);
		elevation = interpolation(arr3[0], arr3[1], arr3[2], arr3[3], intPixelX, intPixelY, pixelX, pixelY);
	}
	else {
		cout << "Error4 : Pixels out of range. Press any button to exit " << endl;
		string s;
		cin >> s;
		exit(4);
	}
}
/// <summary>
///  This funtion first takes 2 points and measures the distance between them. Then returns an 
///  integer value which indicates how many pieces will we have if we get values from every 100 meters.
/// </summary>
/// <param name="sLat"> start point Latitude </param>
/// <param name="sLong"> start point Longitude </param>
/// <param name="eLat"> end point Latitude </param>
/// <param name="eLong">end point Longitude </param>
int size(double sLat, double sLong, double eLat, double eLong)
{
	double const Pideg = 1.74532925199432957692369076848e-2;
	double const earthRadius = 6371.0;
	double distance;
	double dlon = Pideg * (eLong - sLong);
	double dlat = Pideg * (eLat - sLat);
	double a = (sin(dlat / 2) * sin(dlat / 2)) +
		cos(Pideg * (sLat)) * cos(Pideg * (eLat)) * (sin(dlon / 2) * sin(dlon / 2));
	double angle = 2 * asin(sqrt(a));
	distance = angle * earthRadius;
	double temp2 = distance * 10; // distance in 100 meters
	int size = (int)temp2; // size of the array
	return size;
}
void measureDistance(double sLat, double sLong, double eLat, double eLong, double & distance)
{
	double const Pideg = 1.74532925199432957692369076848e-2;
	double const earthRadius = 6371.0;
	double dlon = Pideg * (eLong - sLong);
	double dlat = Pideg * (eLat - sLat);
	double a = (sin(dlat / 2) * sin(dlat / 2)) +
		cos(Pideg * (sLat)) * cos(Pideg * (eLat)) * (sin(dlon / 2) * sin(dlon / 2));
	double angle = 2 * asin(sqrt(a));
	distance = angle * earthRadius;
}
/// <summary>
///  This subroutine calculates angles of a spherical triangle.
/// </summary>
/// <param name="a"> a,b,c:   (sides (radians) of a spherical triangle)</param>
/// <param name="b"> a,b,c:   (sides (radians) of a spherical triangle)</param>
/// <param name="c"> a,b,c:   (sides (radians) of a spherical triangle)</param>
/// <result>
///  angab :  (spherical angle (radians) between sides a and b)
///  angac :  (spherical angle (radians) between sides a and c)
/// </result>
static void Sphtan(double a, double b, double c, double & angab, double & angac)
{
	//calculate the semi-perimeter of the spherical triangle.
	double p = (a + b + c) / 2.0;
	double spa = abs(sin(p - a) / sin(p));
	double spb = abs(sin(p - b));
	double spc = abs(sin(p - c));
	if (spb < 1E-15) spb = 1E-15;
	if (spc < 1E-15) spc = 1E-15;

	//calculate 2 angles of the spherical triangle.
	double t1 = sqrt(spa * spb / spc);
	double t2 = sqrt(spa * spc / spb);

	angab = 2.0 * atan(t1);
	angac = 2.0 * atan(t2);
}

/// <summary>
///   This subroutine calculates the angle an1 between north and point 2
///   as seen from point 1, and the angle an2 between north and point 1
///   as seen from point 2.   an1 and an2 are measured (degs), clockwise
///   from north.   the polar coordinates of point 1 (ph1,th1) and of
///   point 2 (ph2,th2) are given in degs;   the geodesic distance, dist,
///   between points 1 and 2 is given in km.
/// </summary>
/// <param name="ph1"> (longitude and latitude of point 1)</param>
/// <param name="th1"> (longitude and latitude of point 1)</param>
/// <param name="dist">((km) distance between points 1 and 2)</param>
/// <param name="ph2"> (longitude and latitude of point 2)</param>
/// <param name="th2"> (longitude and latitude of point 2)</param>
/// <param name="an1">(angle between north and point 2 as seen from point 1,
///              rotating clockwise)</param>
/// <param name="an2">(angle between north and point 1 as seen from point 2,
///              rotating clockwise)</param>
/// <remarks>constants:
///  pideg=1.74532925199432957692369076848e-2: (conversion factor from degs to rads)
///  pirad=57.295779513082320876798154814105 : (conversion factor from rads to degs) </remarks>
static void Sepang(double ph1, double th1, double dist, double ph2, double th2, double & an1, double & an2)
{
	double const Pideg = 1.74532925199432957692369076848e-2;
	double const Pirad = 57.295779513082320876798154814105;
	double const earthRadius = 6371.0;
	//determine whether linear approximation or spherical calc. is needed
	double thdif = abs(th2 - th1);
	if (thdif < 1E-15)
		thdif = 1E-15;
	double th = th1 > th2 ? th1 : th2;
	double phdif = abs(ph1 - ph2) * cos(th * Pideg);
	double facang = atan(phdif / thdif);
	if (dist * sin(facang) >= 10.0)
	{
		if (abs(ph1 - ph2) >= 1E-10)
		{
			//calculate the sides of the spherical triangle (radians).
			double dot = dist / earthRadius;
			//convert degrees to radians.
			double rlat1 = (90.0 - th1) * Pideg;
			double rlat2 = (90.0 - th2) * Pideg;
			// calculate the angles one, two of the spherical triangle.
			double one, two;
			Sphtan(dot, rlat1, rlat2, one, two);

			//convert degrees to radians.
			one = one * Pirad;
			two = two * Pirad;

			//calc.the wanted angles;  bearings calculated clockwise rel. to north.
			if (ph1 < ph2)
			{
				an1 = one;
				an2 = 360.0 - two;
				return;
			}
			an1 = 360.0 - one;
			an2 = two;
			return;
		}

		// points lie on the same meridian.
		an1 = 0.0;
		an2 = 0.0;
		if (th1 >= th2) an1 = 180.0;
		if (th1 < th2) an2 = 180.0;
		return;
	}

	//small distances; use plane trig.
	double fac = abs((ph2 - ph1) / thdif);
	if (ph2 >= ph1)
	{
		if (th2 >= th1)
		{
			an1 = Pirad * atan(fac * cos(th2 * Pideg));
			an2 = 180 + Pirad * atan(fac * cos(th1 * Pideg));
			return;
		}
		an1 = 180.0 - Pirad * atan(fac * cos(th2 * Pideg));
		an2 = 360.0 - Pirad * atan(fac * cos(th1 * Pideg));
		return;
	}
	if (th2 >= th1)
	{
		an1 = 360.0 - Pirad * atan(fac * cos(th2 * Pideg));
		an2 = 180.0 - Pirad * atan(fac * cos(th1 * Pideg));
		return;
	}
	an1 = 180.0 + Pirad * atan(fac * cos(th2 * Pideg));
	an2 = Pirad * atan(fac * cos(th1 * Pideg));
}

/// <summary>
///  This subroutine calculates the distance and relative orientations
///  between two pointsi PT1 and PT2, on the earth:
///    AN1 is the angle between north and PT2 as seen from PT1;
///    AN2 is the angle between north and PT1 as seen from pt2;
///    AN1 and AN2 are measured (degs), clockwise from north;
///    the geodesic distance, dist, between PT1 and PT2 is given in KM.
///
/// The polar coordinates of PT1 (long1,lati1) and of PT2 (long2,lati2)
/// are input in degress and minutes and seconds.
/// </summary>
/// <param name="lon1"> (decimal) longitude in transmitter</param>
/// <param name="lat1"> (decimal) latitude in transmitter</param>
/// <param name="lon2"> (decimal) longitude in receiver</param>
/// <param name="lat2"> (decimal) latitude in receiver</param>
/// <param name="an1"> output, the angle between north and pt2, viewed from pt1 </param>
/// <param name="an2"> output,  the angle between north and pt1, viewed from pt2</param>

void FindAzimuth(double lat1, double lon1, double lat2, double lon2, double & an1, double &  an2)
{
	double dist;
	measureDistance(lat1, lon1, lat2, lon2, dist);
	Sepang(lon1, lat1, dist, lon2, lat2, an1, an2);
}

/// <summary>
///  This function finds the new coordinate information of a point which has "distance" distance away and "azimuth" azimuth angle away from given point.
/// </summary>
/// <param name="longitude">(decimal) longitude</param>
/// <param name="latitude">(decimal) latitude</param>
/// <param name="azimuth">(degree) azimuth</param>
/// <param name="distance">(km) distance</param>
/// <param name="rLongitude">output, (decimal) longitude</param>
/// <param name="rLatitude">output, (decimal) latitude</param>
void FindNewCoord(double longitude, double latitude, double azimuth, double distance, double &  rLongitude, double &  rLatitude)
{
	double const Pideg = 1.74532925199432957692369076848e-2;
	double const Pirad = 57.295779513082320876798154814105;
	double const earthRadius = 6371.0;
	const double tolerance = 1e-6;
	if (abs(distance) < tolerance)
	{
		rLatitude = latitude;
		rLongitude = longitude;
		return;
	}
	while (azimuth <= 0) azimuth = azimuth + 360;
	while (azimuth >= 360.0) azimuth = azimuth - 360;
	double theta = distance / earthRadius;
	double enlrad = latitude * Pideg;
	double north2 = cos(azimuth * Pideg) * sin(theta) * cos(enlrad) + cos(theta) * sin(enlrad);
	double enl = asin(north2);
	if (abs(enl - 90.0) < 1E-8) enl = 90 - 1E-8;

	if (abs(azimuth) > tolerance)
	{
		double p1 = cos(theta) - sin(enlrad) * sin(enl);
		double p2 = cos(enl) * cos(enlrad);
		double p3;
		if ((p1 / p2) > 1.0) p3 = 1.0;
		else if ((p1 / p2) < -1.0) p3 = -1.0;
		else p3 = p1 / p2;
		double east2 = acos(p3) * Pirad;
		if (azimuth > 180.0) rLongitude = (longitude - east2);
		else rLongitude = (longitude + east2);
	}
	else rLongitude = longitude;
	rLatitude = enl * Pirad;
}
/// <summary>
///  This funtions takes the path of map and 2 points and returns an array 
///  that shows elevation values from each 100 meters
/// </summary>
/// <param name="input"> Path of map</param>
/// <param name="lat1">1st latitude</param>
/// <param name="lon1">1st longitude </param>
/// <param name="lat2"> 2nd lat </param>
/// <param name="lon2">2nd lon </param>
/// <param name="arr">Elevation values from each 100 meters between 2 points </param>
void pathProfile(const char input[], double lat1, double lon1, double lat2, double lon2, double arr[])
{
	GDALDataset * dataset = gdalCreateDataset(input);
	double dist, an1, an2, newLon2, newLat2;
	measureDistance(lat1, lon1, lat2, lon2, dist);
	FindAzimuth(lat1, lon1, lat2, lon2, an1, an2);
	double temp2 = dist * 10; // distance in 100 meters
	int size = (int)temp2; // size of the array

	for (int j = 1; j < size + 1; j++)
	{
		arr[j - 1] = 0;
		double elevation2;
		FindNewCoord(lon1, lat1, an1, 0.1*j, newLon2, newLat2);
		readElevationValInner(dataset, newLon2, newLat2, elevation2);
		arr[j - 1] = elevation2;
	}
	dataset->~GDALDataset();
}

/// <summary>
///  This function finds minimum line of sight between 2 points
/// </summary>
/// <param name="arr">The output array of path profile </param>
/// <param name="size"> size of the input array </param>
/// <param name="height1"> height of the transmitter from ground </param>
/// <param name="height2">height of the receiver from ground </param>
/// <param name="ret[0]">distance (km) </param>
/// <param name="ret[1]">minumun Clearance (m) </param>
/// <param name="ret[2]">control value if = 1 there are smt between 2 points, if = 0 2 points see each other</param>
void los(double arr[],int size, double height1, double height2, double ret[]) 
{
	if (arr == NULL)
		exit(9);
	double h1 = height1, h2 = height2;
	bool ctrl = false;
	h1 += arr[0];
	h2 += arr[size - 1];
	double min = 99999;
	int index = -1;
	double temp;
	if (h1 >= h2)
	{
		for (int j = 0; j < size; j++) {
			temp = h1 - ((h1 - h2) / size)*j;
			temp = temp - arr[j];
			if (temp < min)
			{
				min = temp;
				index = j;
			}
			if (temp < 0)
				ctrl = true;
		}
	}
	else
	{
		for (int j = 0; j < size; j++) {
			temp = h1 + ((h2 - h1) / size)*j;
			temp = temp - arr[j];
			if (temp < min)
			{
				min = temp;
				index = j;
			}
			if (temp < 0)
				ctrl = true;
		}
	}
	if (ctrl)
	{
		ret[0] = index * 0.1;
		ret[1] = min;
		ret[2] = 1;
	}
	else {
		ret[0] = -2;
		ret[1] = -2;
		ret[2]=0;
	}	
}


/// <summary>
///  This function returns an array which has binary values demonstrates visibility of the points on a line from a point (0=visible, 1=not visible)
/// </summary>
/// <param name="input">The output array of path profile </param>
/// <param name="lat1">1st latitude</param>
/// <param name="lon1">1st longitude </param>
/// <param name="lat2"> 2nd lat </param>
/// <param name="lon2">2nd lon </param>
/// <param name="height1"> height of the transmitter from ground </param>
/// <param name="height2">height of the receiver from ground </param>
/// <param name="arr2">control values if = 1 there are smt between 2 points, if = 0 2 points see each other </param>
/// <param name="coorLat"> latitudes of the points in order </param>
/// <param name="coorLon"> longitudes of the points in order</param>
void lineCheck(const char input[], double lat1, double lon1, double lat2, double lon2, double height1, double height2, double arr2[], double coorLat[], double coorLon[])
{
	int s = size(lat1, lon1, lat2, lon2);
	s *= 2;
	double azimuth1 = 0, azimuth2 = 0;
	FindAzimuth(lat1, lon1, lat2, lon2, azimuth1, azimuth2);

	double newLat = 0, newLon = 0, tHeight = 0, rHeight = 0;//transmitter height, receiver height
	double* temp = new double[s];
	pathProfile(input, lat1, lon1, lat2, lon2, temp);
	for (int i = 0; i < s; i++)
	{
		//cout <<"LN-" <<i <<"-"<< s << endl;
		newLat = 0, newLon = 0, tHeight = 0, rHeight = 0;//transmitter height, receiver height
		FindNewCoord(lon1, lat1, azimuth1, i*0.05, newLon, newLat);
		int s2 = size(lat1, lon1, newLat, newLon);
		s2 *= 2;
		coorLat[i] = newLat;
		coorLon[i] = newLon;

		double * my = new double[3];

		los(temp, s2, height1, height2, my);
		arr2[i] = my[2];
		delete[] my;
	}
	delete[] temp;
}

/// <summary>
///  This function repeats lineCheck funtion 360 times(according to angleInterval) and returns the number of repeats
/// </summary>
/// <param name="input">The output array of path profile </param>
/// <param name="lat1">1st latitude</param>
/// <param name="lon1">1st longitude </param>
/// <param name="height1"> height of the transmitter from ground </param>
/// <param name="height2">height of the receiver from ground </param>
/// <param name="radius">radius of the circle (km) </param>
/// <param name="arr2d"> double array of lines and their binary values come from lineCheck method </param>
/// <param name="coor2d"> First half is latitude values and 2nd half is longitude values of the points in arr2d</param>
int shapeFile(const char input[], double lat1, double lon1, double height1, double height2, double radius, double** arr2d, double** coor2d)
{
	GDALDataset * dataset = gdalCreateDataset(input);
	double angleInterval = 1.0;
	int	size = (int)(360.0 / angleInterval); // size of the return array

	double spElev;
	readElevationValInner(dataset, lon1, lat1, spElev);//Start point elevation
	spElev += height1; // Elevation of transmitter

	unsigned int width, height;
	double lux, luy, px, py, an1, an2, nLon, nLat;
	getGEoTransform(dataset, width, height, lux, luy, px, py);
	FindAzimuth(lat1, lon1, luy, lux, an1, an2);
	FindNewCoord(lon1, lat1, an1, radius, nLon, nLat);
	lineCheck(input, lat1, lon1, nLat, nLon, height1, height2, arr2d[0], coor2d[0], coor2d[size]);
	for (int i = 1; i < size; i++)
	{
		double temp = 0;
		if (an1 + i >= 360)
			temp = an1 + (i* angleInterval - 360);
		else
			temp = an1 + (i* angleInterval);
		FindNewCoord(lon1, lat1, temp, radius, nLon, nLat);
		lineCheck(input, lat1, lon1, nLat, nLon, height1, height2, arr2d[i], coor2d[i], coor2d[size + i]);
	}
	return size;
}

/// <summary>
///  This function gets the values frrom shapefile functions and a radius and a string which shows the name of the output file
///  then it creates a .shp file in "C:\\Users\\Ensar\\Documents\\map\\output.shp". You have to change the string according to your pc
/// </summary>
/// <param name="arr2d">double array of lines and their binary values come from shapefile method </param>
/// <param name="coor">First half is latitude values and 2nd half is longitude values of the points in arr2d </param>
/// <param name="radius"> km </param>
void createSHP(double ** arr2d, double ** coor, double radius, string output) {
	const char *pszDriverName = "ESRI Shapefile";
	double angle = 1.0;
	GDALDriver *poDriver;
	GDALAllRegister();
	poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
	if (poDriver == NULL)
	{
		printf("%s driver not available.\n", pszDriverName);
		exit(1);
	}
	GDALDataset *poDS;
	string str = "C:\\Users\\Ensar\\Documents\\map\\";
	str += output;
	str += ".shp";
	cout << str << endl;
	const char * tmp = str.c_str();
	poDS = poDriver->Create(tmp, 0, 0, 0, GDT_Float64, NULL);
	if (poDS == NULL)
	{
		printf("Creation of output file failed.\n");
		exit(1);
	}
	OGRLayer *poLayer;
	const char * tmp2 = output.c_str();
	poLayer = poDS->CreateLayer(tmp2, NULL, wkbPoint, NULL);
	if (poLayer == NULL)
	{
		printf("Layer creation failed.\n");
		exit(1);
	}
	OGRFieldDefn oField("LOS", OFTInteger);
	oField.SetWidth(32);
	if (poLayer->CreateField(&oField) != OGRERR_NONE)
	{
		printf("Creating LOS field failed.\n");
		exit(1);
	}
	int cnt = 0;
	int interval = (int)(360.0 / angle);
	double dist = radius * 20;
	for (int i = 0; i < interval; i++)
	{
		for (int j = 0; j < dist; j++)
		{
			int p = (int)arr2d[i][j];
			if (p == 0)
			{
				OGRFeature *poFeature;
				poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
				poFeature->SetField("LOS", p);
				OGRPoint pt;
				pt.setX(coor[interval + i][j]);
				pt.setY(coor[i][j]);
				poFeature->SetGeometry(&pt);
				if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
				{
					printf("Failed to create feature in shapefile.\n");
					exit(1);
				}
				OGRFeature::DestroyFeature(poFeature);
			}
		}
	}
	GDALClose(poDS);
}
///Taken from apps / gdal library 
static void CreateElevAttrib(const char* pszElevAttrib, OGRLayerH hLayer)
{
	OGRFieldDefnH hFld = OGR_Fld_Create(pszElevAttrib, OFTReal);
	OGRErr eErr = OGR_L_CreateField(hLayer, hFld, FALSE);
	OGR_Fld_Destroy(hFld);
	if (eErr == OGRERR_FAILURE)
	{
		exit(1);
	}
}
/// <summary>
///  This functions creates contour lines of a map
///  
/// </summary>
/// <param name="pszSrcFilename1">input map path </param>
/// <param name="pszDstFilename1"> output shapefile path </param>
/// <param name="dfInterval1"> interval between consecutive lines </param>
int contour(const char * pszSrcFilename1, const char *pszDstFilename1, double dfInterval1)
{
	bool b3D = false;
	int bNoDataSet = FALSE;
	bool bIgnoreNoData = false;
	int nBandIn = 1;
	double dfInterval = dfInterval1;
	double dfNoData = 0.0;
	double dfOffset = 0.0;
	double dfExpBase = 0.0;
	const char *pszSrcFilename = pszSrcFilename1;
	const char *pszDstFilename = pszDstFilename1;
	const char *pszElevAttrib = nullptr;
	const char *pszElevAttribMin = nullptr;
	const char *pszElevAttribMax = nullptr;
	const char *pszFormat = nullptr;
	char **papszDSCO = nullptr;
	char **papszLCO = nullptr;
	//	double adfFixedLevels[1000];
	int nFixedLevelCount = 0;
	const char *pszNewLayerName = "contour";
	bool bQuiet = false;
	GDALProgressFunc pfnProgress = nullptr;
	bool bPolygonize = false;

	GDALAllRegister();
	OGRRegisterAll();

	/* -------------------------------------------------------------------- */
	/*						Open source raster file.                        */
	/* -------------------------------------------------------------------- */

	pszElevAttrib = "elev";
	GDALDatasetH hSrcDS = GDALOpen(pszSrcFilename, GA_ReadOnly);
	if (hSrcDS == nullptr)
		exit(2);

	GDALRasterBandH hBand = GDALGetRasterBand(hSrcDS, 1);
	if (hBand == nullptr)
	{
		CPLError(CE_Failure, CPLE_AppDefined,
			"Band %d does not exist on dataset.",
			1);
		exit(2);
	}

	if (!bNoDataSet && !bIgnoreNoData)
		dfNoData = GDALGetRasterNoDataValue(hBand, &bNoDataSet);


	/* -------------------------------------------------------------------- */
	/*      Try to get a coordinate system from the raster.                 */
	/* -------------------------------------------------------------------- */
	OGRSpatialReferenceH hSRS = nullptr;

	const char *pszWKT = GDALGetProjectionRef(hSrcDS);

	if (pszWKT != nullptr && strlen(pszWKT) != 0)
		hSRS = OSRNewSpatialReference(pszWKT);
	/* -------------------------------------------------------------------- */
	/*      Create the output file.                                         */
	/* -------------------------------------------------------------------- */
	CPLString osFormat;
	if (pszFormat == nullptr)
	{
		std::vector<CPLString> aoDrivers =
			GetOutputDriversFor(pszDstFilename, GDAL_OF_VECTOR);
		if (aoDrivers.empty())
		{
			CPLError(CE_Failure, CPLE_AppDefined,
				"Cannot guess driver for %s", pszDstFilename);
			exit(10);
		}
		else
		{
			if (aoDrivers.size() > 1)
			{
				CPLError(CE_Warning, CPLE_AppDefined,
					"Several drivers matching %s extension. Using %s",
					CPLGetExtension(pszDstFilename), aoDrivers[0].c_str());
			}
			osFormat = aoDrivers[0];
		}
	}
	else
	{
		osFormat = pszFormat;
	}

	OGRSFDriverH hDriver = OGRGetDriverByName(osFormat.c_str());

	if (hDriver == nullptr)
	{
		fprintf(stderr, "Unable to find format driver named %s.\n",
			osFormat.c_str());
		exit(10);
	}

	OGRDataSourceH hDS =
		OGR_Dr_CreateDataSource(hDriver, pszDstFilename, papszDSCO);
	if (hDS == nullptr)
		exit(1);

	OGRLayerH hLayer =
		OGR_DS_CreateLayer(hDS, pszNewLayerName, hSRS,
			bPolygonize ? (b3D ? wkbMultiPolygon25D : wkbMultiPolygon)
			: (b3D ? wkbLineString25D : wkbLineString),
			papszLCO);
	if (hLayer == nullptr)
		exit(1);

	OGRFieldDefnH hFld = OGR_Fld_Create("ID", OFTInteger);
	OGR_Fld_SetWidth(hFld, 8);
	OGR_L_CreateField(hLayer, hFld, FALSE);
	OGR_Fld_Destroy(hFld);

	if (pszElevAttrib)
	{
		CreateElevAttrib(pszElevAttrib, hLayer);
	}

	if (pszElevAttribMin)
	{
		CreateElevAttrib(pszElevAttribMin, hLayer);
	}

	if (pszElevAttribMax)
	{
		CreateElevAttrib(pszElevAttribMax, hLayer);
	}
	/**/
	char** options = nullptr;
	int iIDField = OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(hLayer), "ID");
	int iElevField = (pszElevAttrib == nullptr) ? -1 :
		OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(hLayer),
			pszElevAttrib);

	int iElevFieldMin = (pszElevAttribMin == nullptr) ? -1 :
		OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(hLayer),
			pszElevAttribMin);

	int iElevFieldMax = (pszElevAttribMax == nullptr) ? -1 :
		OGR_FD_GetFieldIndex(OGR_L_GetLayerDefn(hLayer),
			pszElevAttribMax);
	options = CSLAppendPrintf(options, "LEVEL_INTERVAL=%f", dfInterval);
	if (bNoDataSet) {
		options = CSLAppendPrintf(options, "NODATA=%.19g", dfNoData);
	}
	if (iIDField != -1) {
		options = CSLAppendPrintf(options, "ID_FIELD=%d", iIDField);
	}
	if (iElevField != -1) {
		options = CSLAppendPrintf(options, "ELEV_FIELD=%d", iElevField);
	}
	if (iElevFieldMin != -1) {
		options = CSLAppendPrintf(options, "ELEV_FIELD_MIN=%d", iElevFieldMin);
	}
	if (iElevFieldMax != -1) {
		options = CSLAppendPrintf(options, "ELEV_FIELD_MAX=%d", iElevFieldMax);
	}

	CPLErr eErr = GDALContourGenerateEx(hBand, hLayer, options, pfnProgress, nullptr);

	CSLDestroy(options);
	OGR_DS_Destroy(hDS);
	GDALClose(hSrcDS);

	if (hSRS)
		OSRDestroySpatialReference(hSRS);

	CSLDestroy(papszDSCO);
	CSLDestroy(papszLCO);
	GDALDestroyDriverManager();
	OGRCleanupAll();

	return (eErr == CE_None) ? 0 : 1;
}
