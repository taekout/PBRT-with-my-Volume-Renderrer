#if defined(_MSC_VER)
#pragma once
#endif
  
#ifndef PBRT_INTEGRATORS_SINGLE_H
#define PBRT_INTEGRATORS_SINGLE_H

// integrators/single.h*
#include "volume.h"
#include "integrator.h"

class HenyeyGreensteinSampleGenerator;
class MediaBidir;

// SingleScatteringIntegrator Declarations
class SingleScatteringIntegrator : public VolumeIntegrator {
public:
    // SingleScatteringIntegrator Public Methods
    SingleScatteringIntegrator(float ss, int nDirectionsForEyeScatter, int nDirectionsForLightScatter, float rayMarchingStepSize);
	virtual ~SingleScatteringIntegrator(void);
    Spectrum Transmittance(const Scene *, const Renderer *,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene);
    Spectrum Li(const Scene *, const Renderer *, const RayDifferential &ray,
         const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena, MediaBidir *pMediaBidir = NULL) const;
	virtual Spectrum LiForNormalSingleScatterRay(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *transmittance, MemoryArena &arena, MediaBidir *pMediaBidir = NULL) const;
private:
    // SingleScatteringIntegrator Private Data
    float stepSize;
    int tauSampleOffset, scatterSampleOffset;

	// Thesis by Taekyu Shin
	int nDirectionsForEyeScatter;
	int nDirectionsForLightScatter;
	HenyeyGreensteinSampleGenerator *pHGScatterRayGenerator;
	HenyeyGreensteinSampleGenerator *pHGLightRayGenerator;
	HenyeyGreensteinSampleGenerator *pHGSliceSampler;
	float rayMarchingStepSize;
	//Renderer *pMetroRenderer;
};


SingleScatteringIntegrator *CreateSingleScatteringIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_SINGLE_H
