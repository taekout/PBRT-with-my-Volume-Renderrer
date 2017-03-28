#pragma once

// CIE color matching functions
static double colorMatchingFunc[][3]={
	{0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
	{0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
	{0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
	{0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
	{0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
	{0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
	{0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
	{0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
	{0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
	{0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
	{0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
	{0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
	{0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
	{0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
	{1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
	{1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
	{0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
	{0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
	{0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
	{0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
	{0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
	{0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
	{0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
	{0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
	{0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
	{0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
	{0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
};

const double wavelengthMin = 380;
const double wavelengthMax = 780;
const double wavelengthStep = 5;
const unsigned numSteps = (wavelengthMax - wavelengthMin) / wavelengthStep + 1;

template<typename T>
class Color3
{
public:    
	Color3() : r(0), g(0), b(0) {}
	Color3( T rr, T gg, T bb) : r(rr), g(gg), b(bb) {}
	T r, g, b;
	Color3 operator / (double w) {
		return Color3(r/w, g/w, b/w);
	}
	Color3 operator * (double w) {
		return Color3(r*w, g*w, b*w);
	}
	Color3 operator *= (float w) {
		r *= w; g *= w; b *= w;
		return Color3(r,g,b);
	}
	bool operator == (Color3 c) {
		return (c.r == r) && (c.g == g) && (c.b == b);
	}
	bool operator != (Color3 c) {
		return (  (c.r != r) || (c.g != g) || (c.b != b) );
	}
	inline Color3 & ClampNegative() {
		r = (r < 0.) ? 0. : r;
		g = (g < 0.) ? 0. : g;
		b = (b < 0.) ? 0. : b;
		return *this;
	}
};

typedef Color3<double>	Color3d;
typedef double			wavelength;

class FireColor
{
public:
	FireColor(void);
	virtual ~FireColor(void);
	inline double ComputePlanck(const double &T, const double &lambda);	// temperature (Kelvin), wavelength (meter)
	Color3d XYZ2sRGB(const double &X, const double &Y, const double &Z);
	Color3d XYZ2RGBStandardByCIE(const double &X, const double &Y, const double &Z);
	Color3d ComputeBlackbody(const double &temperature, double &power);
	void ComputeBlackbodyRamp();
	float ComputeStefanBoltzmann(const float &temperature);
	Color3d ComputePlanckLocus(const double & T);
	wavelength WienDisplacement(const double & T);
};

