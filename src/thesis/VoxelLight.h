#if defined(_MSC_VER)
#pragma once
#endif

#ifndef THESIS_VOXELLIGHT
#define THESIS_VOXELLIGHT

// lights/point.h*
#include "pbrt.h"
#include "light.h"
#include "shape.h"
#include "../thesis/FireColor.h"

// VoxelLight Declarations
class VoxelLight : public Light {
public:
    // VoxelLight Public Methods
    VoxelLight(const Transform &light2world, const Spectrum &intensity, const BBox & voxelExtent, RNG * rng);
	virtual ~VoxelLight();
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
	Point GetPoint(void) const;
	Point GetCenterPoint(void) const;
	inline Spectrum GetIntensity(void) { return Intensity; }
private:
    // VoxelLight Private Data
	BBox extent;
	Point extentSize;
	Point midPoint;
    Spectrum Intensity;
	RNG *rng;
};


/* VoxelLight *CreateVoxelLight(const Transform &light2world,
        const ParamSet &paramSet); */
VoxelLight *CreateVoxelLight(const Transform &light2world,
        const ParamSet &paramSet, vector<BBox> &lightPos, vector<Light*> *multipleVoxelLights);
VoxelLight *CreateFireVoxelLights(VolumeRegion * vgd, const Transform &light2world, const ParamSet &paramSet,
	vector<BBox> & lightVoxels, vector<Color3d> &lightColor, vector<Light*> *multipleVoxelLights);

#endif // THESIS_VOXELLIGHT
