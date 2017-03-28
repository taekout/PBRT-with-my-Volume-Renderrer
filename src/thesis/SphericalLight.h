#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_DIFFUSE_H
#define PBRT_LIGHTS_DIFFUSE_H

// lights/diffuse.h*
#include "pbrt.h"
#include "light.h"
#include "primitive.h"
#include "../shapes/sphere.h"

class HenyeyGreensteinSampleGenerator;

// SphericalLight Declarations
class SphericalLight {
public:
	// SphericalLight Public Methods
	SphericalLight(Point center, float radius, Spectrum le, HenyeyGreensteinSampleGenerator *pHGMultipleScatteringRayGenerator);
	virtual ~SphericalLight();

	Spectrum L(const Point &p, const Normal &n, const Vector &w) const {
		return Dot(n, w) > 0.f ? Lemit : 0.f;
	}
	Spectrum Power();
	bool IsDeltaLight() const { return false; }
	float Pdf(const Point &, const Vector &) const;
	Spectrum Sample_L(const Point &P, float pEpsilon, const LightSample &ls, float time,
		Vector *wo, float *pdf, VisibilityTester *visibility) const;
	Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
		float time, Ray *ray, Normal *Ns, float *pdf) const;


	//TAekyu Shin
	inline Point GetLightPos() {
		return lightPos;
	}
	void SetLightPos(Point pos) {
		lightPos = pos;
	}

protected:
	// SphericalLight Protected Data
	Spectrum Lemit;	// Lemit = L(Light Intensity) * sc(scale)
	Sphere *sphere;
	float area;

	Point lightPos;
	float radius;
	RNG *rng;

	HenyeyGreensteinSampleGenerator * hgSampleGenerator;
	int nSampleHG;
};

#endif // PBRT_LIGHTS_DIFFUSE_H
