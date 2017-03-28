#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_POINT_H
#define PBRT_LIGHTS_POINT_H

// lights/point.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
#include "../thesis/FireColor.h"

// PointLight Declarations
class PointLight : public Light {
public:
    // PointLight Public Methods
    PointLight(const Transform &light2world, const Spectrum &intensity, float backGroundI, int bVoxelLight);
	virtual Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *vis) const;
    virtual Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    virtual Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1,
                      float u2, float time, Ray *ray, Normal *Ns, float *pdf) const;
    virtual float Pdf(const Point &, const Vector &) const;
    virtual void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVisibility, float time, RNG &rng, Spectrum *coeffs) const;

	// Thesis by Taekyu Shin
	inline Point GetPoint(void) { return lightPos; }
	inline Spectrum GetIntensity(void) { return Intensity; }
	inline void SetSampleID(int nID) { nStratifiedSampleID = nID; }

	virtual Spectrum Le(const RayDifferential &r) const;
private:
    // PointLight Private Data
    Point lightPos;
    Spectrum Intensity;
	int bVoxelLight; // 0 : Pointlight, 1: VoxelLight, 2: Spherical Light
	int nStratifiedSampleID;
	RNG m_rng;

	float backgroundI;
};


/* PointLight *CreatePointLight(const Transform &light2world,
        const ParamSet &paramSet); */
PointLight *CreatePointLight(const Transform &light2world,
        const ParamSet &paramSet, vector<Point> &lightPos, vector<Light*> *multiplePointLights);
PointLight *CreateFirePointLights(VolumeRegion * vgd, const Transform &light2world, const ParamSet &paramSet,
	vector<Point> &lightPos, vector<Color3d> &lightColor, vector<Light*> *multiplePointLights);

#endif // PBRT_LIGHTS_POINT_H
