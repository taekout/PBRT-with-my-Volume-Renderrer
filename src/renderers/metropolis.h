#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_METROPOLIS_H
#define PBRT_RENDERERS_METROPOLIS_H

// renderers/metropolis.h*
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
struct MLTSample;
class DirectLightingIntegrator;
class MediaBidir;
struct LightingSample;

// Metropolis Declarations
struct PathVertex;

class SamplerRenderer; //////////////////////////////////// Give this a thought!!!!!!!!!!!!!!!! You might wanna remove this and change all 'SamplerRenderer' to just 'Renderer'.

class MetropolisRenderer : public Renderer {
public:
    // MetropolisRenderer Public Methods
    MetropolisRenderer(int perPixelSamples, int nBootstrap,
        int directPixelSamples, float largeStepProbability,
        bool doDirectSeparately, int maxConsecutiveRejects, int maxDepth,
        Camera *camera, bool doBidirectional, SamplerRenderer * p);
    ~MetropolisRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // MetropolisRenderer Private Methods
    Spectrum PathL(const MLTSample &sample, const Scene *scene,
        MemoryArena &arena, const Camera *camera,
        const Distribution1D *lightDistribution, PathVertex *cameraPath,
        PathVertex *lightPath, RNG &rng) const;
    Spectrum Lpath(const Scene *scene, const PathVertex *path, int pathLength,
        MemoryArena &arena, const vector<LightingSample> &samples,
        RNG &rng, float time, const Distribution1D *lightDistribution,
        const RayDifferential &escapedRay, const Spectrum &escapedAlpha) const;
    Spectrum Lbidir(const Scene *scene,
        const PathVertex *cameraPath, int cameraPathLength,
        const PathVertex *lightPath, int lightPathLength,
        MemoryArena &arena, const vector<LightingSample> &samples,
        RNG &rng, float time, const Distribution1D *lightDistribution,
        const RayDifferential &escapedRay, const Spectrum &escapedAlpha) const;

    // MetropolisRenderer Private Data
    Camera *camera;
    bool bidirectional;
    uint32_t nDirectPixelSamples, nPixelSamples, maxDepth;
    uint32_t largeStepsPerPixel, nBootstrap, maxConsecutiveRejects;
    DirectLightingIntegrator *directLighting;
    AtomicInt32 nTasksFinished;
    friend class MLTTask;

public:
	Renderer * samplerRenderer;
};


MetropolisRenderer *CreateMetropolisRenderer(const ParamSet &params,
    Camera *camera, SamplerRenderer *pSamplerRenderer);
MediaBidir *CreateMetropolisRendererForMedia(const ParamSet &params, MetropolisRenderer *pRenderer);
#endif // PBRT_RENDERERS_METROPOLIS_H
