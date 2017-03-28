#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_SAMPLERRENDERER_H
#define PBRT_RENDERERS_SAMPLERRENDERER_H

// renderers/samplerrenderer.h*
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
#include <time.h>

class MediaBidir;

// SamplerRenderer Declarations
class SamplerRenderer : public Renderer {
public:
    // SamplerRenderer Public Methods
    SamplerRenderer(Sampler *s, Camera *c, SurfaceIntegrator *si,
                    VolumeIntegrator *vi, bool visIds, MediaBidir * mediaRenderer);
    ~SamplerRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
private:
    // SamplerRenderer Private Data
    bool visualizeObjectIds;
    Sampler *sampler;
    Camera *camera;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
	MediaBidir *mediaRenderer;
};



// SamplerRendererTask Declarations
class SamplerRendererTask : public Task {
public:
    // SamplerRendererTask Public Methods
    SamplerRendererTask(const Scene *sc, Renderer *ren, Camera *c,
                        ProgressReporter &pr, Sampler *ms, Sample *sam, 
                        bool visIds, int tn, int tc)
      : reporter(pr)
    {
        scene = sc; renderer = ren; camera = c; mainSampler = ms;
        origSample = sam; visualizeObjectIds = visIds; taskNum = tn; taskCount = tc;
		srand(time(NULL));
    }
    void Run();
private:
    // SamplerRendererTask Private Data
    const Scene *scene;
    const Renderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    bool visualizeObjectIds;
    int taskNum, taskCount;
};



#endif // PBRT_RENDERERS_SAMPLERRENDERER_H
