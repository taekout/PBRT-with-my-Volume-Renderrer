#pragma once
#include "spectrum.h"
#include "memory.h"

class PointLight;

namespace mipmapdata {
typedef enum {
	TEXTURE_REPEAT,
	TEXTURE_BLACK,
	TEXTURE_CLAMP
} Wrapping;
}

struct RGBColor
{
	float c[3];
	RGBColor() {
		c[0] = c[1] = c[2] = 0.f;
	}
	RGBColor operator+ (const RGBColor &s2) const {
		RGBColor ret = *this;
		for(int i = 0 ; i < 3 ; i++) {
			ret.c[i] += s2.c[i];
		}
		return ret;
	}
	RGBColor operator* (const RGBColor &s2) const {
		RGBColor ret = *this;
		for(int i = 0 ; i < 3 ; i++) {
			ret.c[i] *= s2.c[i];
		}
		return ret;
	}
	RGBColor operator* (const float s2) const {
		RGBColor ret = *this;
		for(int i = 0 ; i < 3 ; i++) {
			ret.c[i] *= s2;
		}
		return ret;
	}

	operator Spectrum() const;
};

class TMipmap
{
	
protected:
	static const int initValue = -1000;

	int level;
	int nx, ny, nz;
	MemoryArena arena;

public:
	std::vector<RGBColor ***> levelMaps;

	TMipmap(int _nx, int _ny, int _nz, int _level);
	virtual ~TMipmap(void);
	void InitializeMap(RGBColor *** volume);

	void Blur2X2X2(RGBColor ***srcMap, RGBColor ***dstMap, int _z, int _y, int _x);

	int X() const { return nx; }
	int Y() const { return ny; }
	int Z() const { return nz; }
	int Levels() const { return level; }
	Spectrum LookUp(int level, float _nx, float _ny, float _nz);
	Spectrum LAttenuatedBasedOnDistance(Point p, PointLight * light);

private:
	void Halve(int &x, int &y, int &z) {
		if(z % 2 == 0) z /= 2;
		else	z = z / 2 + 1;
		if(y % 2 == 0) y /= 2;
		else	y = y / 2 + 1;
		if(x % 2 == 0) x /= 2;
		else	x = x / 2 + 1;
	}
};

