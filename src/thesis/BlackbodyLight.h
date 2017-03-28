#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_BLACKBODY_H
#define PBRT_LIGHTS_BLACKBODY_H

// lights/point.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
#include "../thesis/FireColor.h"

// BlackbodyLight Declarations
class BlackbodyLight : public Light {
public:
    // BlackbodyLight Public Methods
    BlackbodyLight(const Transform &light2world, const Spectrum &intensity);
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *vis) const;
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1,
                      float u2, float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVisibility, float time, RNG &rng, Spectrum *coeffs) const;

	// Thesis by Taekyu Shin
	inline Point GetPoint(int lightIndex) { return *lightPositions.at(lightIndex); }
	inline Spectrum GetIntensity(int lightIndex) { return *lightIntensity.at(lightIndex); }
	void AddLight(Point p, Spectrum intensity, Spectrum s_s, Spectrum s_a, Spectrum den);
private:
    // BlackbodyLight Private Data
	int nLights;
	std::vector<Point *> lightPositions; // Later use KD-tree to optimize. For ex, I could put together those lights. To do that, I need use std::map.
	// See Wiki KD-tree. See section (Volumetric Objects). It talks about orthogonal range search. I could try that.
	std::vector<Spectrum *> lightIntensity;
	std::vector<Spectrum *> sigma_s;
	std::vector<Spectrum *> sigma_a;
	std::vector<Spectrum *> sootDensity;
};


/* BlackbodyLight *CreateBlackbodyLight(const Transform &light2world,
        const ParamSet &paramSet); */
 BlackbodyLight *CreateBlackbodyLight(const Transform &light2world,
        const ParamSet &paramSet, vector<Point> &lightPos, vector<Light*> *multipleBlackbodyLights);
BlackbodyLight *CreateFireLights(const Transform &light2world, const ParamSet &paramSet,
	vector<Point> &lightPos, vector<Color3d> &lightColor, vector<Light*> *multipleBlackbodyLights);

#endif // PBRT_LIGHTS_BLACKBODY_H
