import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class Transformer
{
	public static void main(String[] args)
	{
		if (args.length != 2)
		{
			System.out.println("Usage: java Transformer <latitude> <longitude>");
			System.exit(0);
		}
		else
		{
			double latitude = Double.parseDouble(args[0]);
			double longitude = Double.parseDouble(args[1]);
			Transformer transformer = new Transformer(latitude, longitude);
		}
	}
    
	public Transformer(double latitude, double longitude)
	{
		int zone = getZoneNumber(longitude);
		char zone_letter = getZoneLetter(latitude);
	    ArrayList origin = getZoneOrigin(latitude, longitude);
	    System.out.println("Latitude " + latitude + ", longitude " + longitude + " is in zone " + zone + " " + zone_letter);
	    
	    double x_origin = (double)origin.get(1);
	    double y_origin = (double)origin.get(0);
	    
	    System.out.println("The longitude of the zone origin is " + x_origin);
	    System.out.println("The latiude of the zone origin is " + y_origin);
	    
	    double longitude_degree_length = getLongitudeDegreeLength(y_origin);
	    longitude_degree_length       += getLongitudeDegreeLength(latitude);
	    longitude_degree_length       /= 2;
	    
	    double latitude_degree_length = getLatitudeDegreeLength(y_origin);
	    latitude_degree_length       += getLatitudeDegreeLength(latitude);
	    latitude_degree_length       /= 2;
	    
	    double longitude_delta = longitude - x_origin;
	    
    	double east = longitude_delta * longitude_degree_length;
    	east       += 500000;
    	
    	double north = latitude * latitude_degree_length;
    	if(latitude < 0)
    		north += 10000000;
	    
	    System.out.println("Latitude " + latitude + ", longitude " + longitude + " is in zone " + zone + " " + zone_letter);
	    System.out.println("The location is " + String.format("%.4f", east) + " meters east of the false zone origin.");
	    System.out.println("The location is " + String.format("%.4f", north) + " meters north of the false zone origin.");
	    
	}
	
	public double getHaversinDistance(double lat1, double long1, double lat2, double long2)
	{
	    double lat_delta  = lat2 - lat1;
	    double long_delta = long2 - long1;
	    
	    lat1 = Math.toRadians(lat1);
	    lat2 = Math.toRadians(lat2);
	    lat_delta = Math.toRadians(lat_delta);
	    long_delta = Math.toRadians(long_delta);
	    
	    double a = Math.sin(lat_delta / 2) * Math.sin(lat_delta / 2) + Math.sin(long_delta / 2) * Math.sin(long_delta / 2) * Math.cos(lat1) * Math.cos(lat2);
	    double b = Math.sqrt(a);
	    double c = Math.sqrt(1 - a);
	    double d = 2 * Math.atan2(b, c);
	    
	    double radius = 6367000;
		return(d * radius);
	}
	
	public double getVincentyDistance(double lat1, double lon1, double lat2, double lon2) {

		  double a = 6378137, b = 6356752.314245, f = 1 / 298.257223563;

		  double L = Math.toRadians(lon2 - lon1);

		  double U1 = Math.atan((1 - f) * Math.tan(Math.toRadians(lat1)));

		  double U2 = Math.atan((1 - f) * Math.tan(Math.toRadians(lat2)));

		  double sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);

		  double sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);

		  double cosSqAlpha, sinSigma, cos2SigmaM, cosSigma, sigma;

		  double lambda = L, lambdaP, iterLimit = 100;

		  do {

		    double sinLambda = Math.sin(lambda), cosLambda = Math.cos(lambda);

		    sinSigma = Math.sqrt( (cosU2 * sinLambda)

		      * (cosU2 * sinLambda)

		      + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda)

		      * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda)

		    );

		    if (sinSigma == 0) return 0;

		    cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

		    sigma = Math.atan2(sinSigma, cosSigma);

		    double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;

		    cosSqAlpha = 1 - sinAlpha * sinAlpha;

		    cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

		    double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

		    lambdaP = lambda;

		    lambda =  L + (1 - C) * f * sinAlpha  
		      * (sigma + C * sinSigma 
		        * (cos2SigmaM + C * cosSigma
		          *(-1 + 2 * cos2SigmaM * cos2SigmaM)
		        )
		      );

		  } while (Math.abs(lambda - lambdaP) > 1e-12 && --iterLimit > 0);

		  if (iterLimit == 0) return 0;

		  double uSq = cosSqAlpha * (a * a - b * b) / (b * b);

		  double A = 1 + uSq / 16384

		      * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));

		  double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

		  double deltaSigma = 

		    B * sinSigma

		      * (cos2SigmaM + B / 4

		        * (cosSigma 

		          * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM

		            * (-3 + 4 * sinSigma * sinSigma)

		              * (-3 + 4 * cos2SigmaM * cos2SigmaM)));

		  double s = b * A * (sigma - deltaSigma);
		  return s;
		}
	
	    double getLatitudeDistance(double latitude_degrees1, double latitude_degrees2, int number_of_segments)
	    {
	    	double delta     = latitude_degrees2 - latitude_degrees1;
	    	double increment = delta / number_of_segments;
	    	
	    	double distance = 0;
	    	double current = latitude_degrees1;
	    	for(int i = 0; i < number_of_segments; i++)
	    	{
	    	    double degree_length1 = getLatitudeDegreeLength(current);
	    	    double degree_length2 = getLatitudeDegreeLength(current + increment);
	    	    double average        = (degree_length1 + degree_length2) / 2;
	    	    distance             += average * increment;
	    	    
	    	    current += increment;
	    	}
	    	
	    	return distance;
	    }
	 
	    int getZoneNumber(double longitude_degrees)
	    {
	    	int zone = (int)Math.floor(31 + longitude_degrees / 6);
	    	return(zone);
	    }
	    
	    char getZoneLetter(double longitude_degrees)
	    {
	    	longitude_degrees += 80;
	    	longitude_degrees /= 8;
	    	longitude_degrees = Math.floor(longitude_degrees);
	    	
	    	byte char_value = (byte)longitude_degrees;
	    	
	    	if(char_value < 6)
	    	{
	    	    char_value += 67;	
	    	}
	    	else if(char_value < 11)
	    	{
	    	    char_value += 68;	
	    	}
	    	else
	    	{
	    	    char_value += 69;	
	    	}
	    	
	    	char zone_letter = (char)char_value;
	    	return zone_letter;
	    }
	    
	    ArrayList getZoneOrigin(double latitude, double longitude)
	    {
	    	ArrayList origin = new ArrayList();
	    	
	    	
	    	double latitude_origin = 0;
	    	int zone = getZoneNumber(longitude);
	    	double longitude_origin = (zone - 31) * 6 + 3;
	    	
	    	origin.add(latitude_origin);
	    	origin.add(longitude_origin);
	    	
	    	return origin;
	    }
	    
	    double getLatitudeDegreeLength(double latitude_degrees)
	    {
	        double latitude_radians = Math.toRadians(latitude_degrees);	
	        
	        double length = 111132.92 - 559.82 * Math.cos(2 * latitude_radians) + 1.175 * Math.cos(4 * latitude_radians) 
	                        - 0.0023 * Math.cos(6 * latitude_radians);
	        return length;
	    }
	 
	    double getLatitudeDegreeLength2(double latitude_degrees)
	    {
	        double latitude_radians = Math.toRadians(latitude_degrees);	
	        double radius           = 6371000;
	        double circumference    = 2 * Math.PI * radius * Math.cos(latitude_radians);
	        double length           = circumference / 360;
	        
	        return 2 * length;
	    }
	    
	    double getLongitudeDegreeLength(double latitude_degrees)
	    {
	        double latitude_radians = Math.toRadians(latitude_degrees);	
	        
	        double length = 111412.84 * Math.cos(latitude_radians) - 93.5 * Math.cos(3 * latitude_radians) 
	                        + 0.118 * Math.cos(5 * latitude_radians);
	        return length;
	    }
	    
	    double getLongitudeDegreeLength2(double latitude_degrees)
	    {
	    	double latitude_radians = Math.toRadians(latitude_degrees);	
	        double radius           = 6371000;
	        double circumference    = 2 * Math.PI * radius * Math.cos(latitude_radians);
	        double length           = circumference *  Math.cos(latitude_radians) / 360;
	        
	        return 2 * length;
	    }
	    
	    double getCurvature(double latitude_degrees)
	    {
	    	double radius           = 6378137.;
	    	double esquare          = 0.081819191 * 0.081819191;
	    	double latitude_radians = Math.toRadians(latitude_degrees);
	    	double latitude_sin     = Math.sin(latitude_radians);
	    	double lsquare          = latitude_sin * latitude_sin;
	    	double root             = Math.sqrt(1. - esquare * lsquare);
	    	double cubic            = root * root * root;
	    	double curvature        = radius * ( 1. - esquare / cubic);
	    	
	    	return curvature;
	    }
	    
	    double getScaleFactor(double easting, double latitude_degrees)
	    {
	    	double scale_factor     = 0.9996;
	    	double curvature        = getCurvature(latitude_degrees);
	    	double latitude_radians = Math.toRadians(latitude_degrees);
	    	
	    	
	    	scale_factor *= easting - 500000.;
	    	scale_factor *= Math.sin(latitude_radians);
	    	scale_factor *= Math.tan(latitude_radians);
	    	scale_factor /= 2 * curvature * curvature;
	    	
	    	return scale_factor;
	    }
}