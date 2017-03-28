
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

// core/integrator.h*
#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "light.h"
#include "reflection.h"
#include "sampler.h"
#include "material.h"
#include "probes.h"
#include "renderer.h"

// Integrator Declarations
class Integrator {
public:
    // Integrator Interface
    virtual ~Integrator();
    virtual void Preprocess(const Scene *scene, const Camera *camera,
                            const Renderer *renderer) {
    }
    virtual void RequestSamples(Sampler *sampler, Sample *sample,
                                const Scene *scene) {
    }
};


class SurfaceIntegrator : public Integrator {
public:
    // SurfaceIntegrator Interface
    virtual Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const = 0;
};


Spectrum UniformSampleAllLights(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Point &p, const Normal &n, const Vector &wo,
    float rayEpsilon, float time, BSDF *bsdf, const Sample *sample, RNG &rng,
    const LightSampleOffsets *lightSampleOffsets,
    const BSDFSampleOffsets *bsdfSampleOffsets);
Spectrum UniformSampleOneLight(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Point &p, const Normal &n, const Vector &wo,
    float rayEpsilon, float time, BSDF *bsdf,
    const Sample *sample, RNG &rng, int lightNumOffset = -1,
    const LightSampleOffsets *lightSampleOffset = NULL,
    const BSDFSampleOffsets *bsdfSampleOffset = NULL);
Spectrum EstimateDirect(const Scene *scene, const Renderer *renderer,
    MemoryArena &arena, const Light *light, const Point &p,
    const Normal &n, const Vector &wo, float rayEpsilon, float time, const BSDF *bsdf,
    RNG &rng, const LightSample &lightSample, const BSDFSample &bsdfSample,
    BxDFType flags);
Spectrum SpecularReflect(const RayDifferential &ray, BSDF *bsdf, RNG &rng,
    const Intersection &isect, const Renderer *renderer, const Scene *scene,
    const Sample *sample, MemoryArena &arena);
Spectrum SpecularTransmit(const RayDifferential &ray, BSDF *bsdf, RNG &rng,
    const Intersection &isect, const Renderer *renderer, const Scene *scene,
    const Sample *sample, MemoryArena &arena);
Distribution1D *ComputeLightSamplingCDF(const Scene *scene);

#endif // PBRT_CORE_INTEGRATOR_H
