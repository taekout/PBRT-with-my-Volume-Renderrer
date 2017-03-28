#include "../core/stdafx.h"
#include "FireColor.h"

FireColor::FireColor(void)
{
}


FireColor::~FireColor(void)
{
}

double FireColor::ComputePlanck(
	const double &T,        // temperature (Kelvin)
	const double &lambda)   // wavelength (meter)
{
	static const double h = 6.62606896e-34;   // Plank constant
	static const double c = 2.99792458e+8;    // Speed of light
	static const double k = 1.38064880e-23;   // Boltzmann constant
	static const double arg1 = 2 * M_PI * h * c * c;
	static const double arg2 = ( h * c ) / k;
	return (arg1 * pow(lambda, -5.0)) / (exp(arg2 / (lambda * T)) - 1.0);
}

// convert XYZ tristimulus values to RGB
Color3d FireColor::XYZ2sRGB(const double &X, const double &Y, const double &Z)
{
	// convert XYZ to sRGB
	return Color3d(
		X *  3.2406 + Y * -1.5372 + Z * -0.4986,
		X * -0.9689 + Y *  1.8758 + Z *  0.0415,
		X *  0.0557 + Y * -0.2040 + Z *  1.0570 );
}

Color3d FireColor::XYZ2RGBStandardByCIE(const double &X, const double &Y, const double &Z)
{
	// convert XYZ to sRGB
	return Color3d(
		X *  3.1956 + Y * 2.4478 + Z * -0.1434,
		X * -2.5455 + Y * 7.0492 + Z * 0.9963,
		Z );
}

Color3d FireColor::ComputeBlackbody(const double &temperature, double &power)
{
	power = 0;
	double X = 0, Y = 0, Z = 0;
	for (unsigned k = 0; k < numSteps; k++ ) {
		// convert to nanometre
		double lambda = (wavelengthMin + k * wavelengthStep) * 1e-9;
		double I = ComputePlanck(temperature, lambda);
		power += I * wavelengthStep;
		X += I * colorMatchingFunc[k][0];
		Y += I * colorMatchingFunc[k][1];
		Z += I * colorMatchingFunc[k][2];
	}
	power /= double( numSteps );

	// normalise the result
	double nor = 1 / std::max( X, std::max( Y, Z ) );

	X *= nor, Y *= nor, Z *= nor;

	return XYZ2RGBStandardByCIE(X, Y, Z);
}

void FireColor::ComputeBlackbodyRamp()
{
	unsigned rampWidth = 512, rampHeight = 50;
	unsigned tempMin = 1000;
	unsigned tempMax = 10000;
	Color3d *result = new Color3d[ rampWidth ];
	double power, temperature;
	for (unsigned i = 0; i < rampWidth; ++i) {
		temperature = tempMin + (tempMax - tempMin) * i / double(rampWidth - 1);
		result[ i ] = ComputeBlackbody(temperature, power );
	}
	// convert double rgb to unsigned bytes
	unsigned numBytes = 3 * rampWidth;
	unsigned char *lines = new unsigned char[ numBytes ];
	unsigned r, g, b;
	double scale = 1 / 2.6; // arbitrary division to keep values < 1
	for (unsigned i = 0, k = 0; i < rampWidth; ++i, k += 3) {
		r = (unsigned char)(std::max( 0., std::min( 1., result[i].r * scale ) ) * 255 + 0.5);
		g = (unsigned char)(std::max( 0., std::min( 1., result[i].g * scale ) ) * 255 + 0.5);
		b = (unsigned char)(std::max( 0., std::min( 1., result[i].b * scale ) ) * 255 + 0.5);
		lines[ k ] = r, lines[ k + 1 ] = g, lines[ k + 2 ] = b;
	}
	// save the result to PPM file
	std::ofstream ofs("./blackbodyramp.ppm");
	ofs << "P6\n" << rampWidth << " " << rampHeight << "\n255\n";
	for (unsigned j = 0; j < rampHeight; ++j) {
		ofs.write((char*)lines, numBytes);
	}
	ofs.close();
	if(result)
		delete [] result;
	if(lines)
		delete [] lines;
}

float FireColor::ComputeStefanBoltzmann(const float &temperature)
{
	return pow(temperature, 4.0f) * 0.000000056703f;
}

Color3d FireColor::ComputePlanckLocus(const double & T)
{
	double xc, yc, arg1 = 1e9 / (T * T * T), arg2 = 1e6 / (T * T), arg3 = 1e3 / T;
	if ( T <= 1000.0 ) {
		return Color3d(0.,0.,0.);
	}
	else if ( T > 1000. && T < 1667.0 ) {
		// Use Wien's Displacement Law to approximate the wavelength.
		//WienDisplacement(T);
		xc = -0.2661239 * arg1 - 0.2343580 * arg2 + 0.8776956 * arg3 + 0.179910;

		double u = (0.860117757 + 1.54118254 * 1e-4 * T + 1.28641212 * 1e-7 * T * T) /
			(1.0 + 8.42420235 * 1e-4 * T + 7.08145163 * 1e-7 * T * T);
		double v = (0.317398726 + 4.22806245 * 1e-5 * T + 4.20481691 * 1e-8 * T * T) /
			(1.0 - 2.89741816 * 1e-5 * T + 1.61456053 * 1e-7 * T * T);
		double x = (-6.0 * u * v) / (4*v*(4.0*v - 2.0 - u*v));
		double y = (4.0*x / u + 2.0*x - 3.0) / 12.0;
		/*double y = ( (1.882 * u + 5.5932) * (2.9088 / 1.882) - 2.9088 * u) / (2.0 * u - 1.9116) / (1.0 - (1.882 * u + 5.5932) * (12.0 * v - 7.8972) / (2.0 * u - 1.9116));
		double x = (2.0 * y * v - 7.8972 * y + 2.9088 * v) / (1.882 * v);*/
		double z = 1.0 - x - y;
		// Problem! : I know x, y. Then I can get X, Y, Z if I know Y(luminance). How to determine Y?
		// Once I know X, Y, Z.  I can just use some standard conversion 3x3 matrix from XYZ to RGB.

		Color3d resultcolor;
		return XYZ2RGBStandardByCIE(x, y, z).ClampNegative();
		
		//return resultcolor;
	}
	else if ( T >= 1667 && T <= 4000 ) {
		xc = -0.2661239 * arg1 - 0.2343580 * arg2 + 0.8776956 * arg3 + 0.179910;
	}
	else if ( T > 4000 && T <= 25000 ) {
		xc = -3.0258469 * arg1 + 2.1070379 * arg2 + 0.2226347 * arg3 + 0.240390;
	}
	double xc3=xc*xc*xc, xc2=xc*xc;
	if ( T > 1000.0 && T < 1667.0 ) {
		yc = -1.1063814 * xc3 - 1.34811020 * xc2 + 2.18555832 * xc - 0.20219683;
	}
	else if ( T < 1667.0 ) {
		yc = -1.1063814 * xc3 - 1.34811020 * xc2 + 2.18555832 * xc - 0.20219683;
	}
	if ( T >= 1667 && T <= 2222 ) {
		yc = -1.1063814 * xc3 - 1.34811020 * xc2 + 2.18555832 * xc - 0.20219683;
	}
	else if ( T > 2222 && T <= 4000 ) {
		yc = -0.9549476 * xc3 - 1.37418593 * xc2 + 2.09137015 * xc - 0.16748867;
	}
	else if ( T > 4000 && T <= 25000 ) {
		yc = +3.0817580 * xc3 - 5.87338670 * xc2 + 3.75112997 * xc - 0.37001483;
	}
	double X = xc / yc;
	double Y = 1;
	double Z = (1 - xc - yc) / yc; // Y = 1.0, that is why I don't multiply it by Y.
	return XYZ2RGBStandardByCIE(X, Y, Z).ClampNegative();
}

wavelength FireColor::WienDisplacement(const double & T)
{
	return 0.;
}

