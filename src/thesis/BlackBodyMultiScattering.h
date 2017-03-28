#pragma once
#include "pbrt.h"
#include "spectrum.h"
#include "TMipmap.h"
#include "../src/volumes/volumegrid.h"
#include "integrator.h"
#include "montecarlo.h"

class BlackBodyMultiScattering : public VolumeIntegrator
{
public:
	Mutex * mutex;
	Preprocessor * mPreprocessor;
	Scene * mScene;
	VolumeGridDensity * mVolume;
	float stepSize;

private:
	TMipmap *mipmap;
	MemoryArena arena;

	int tauSampleOffset, scatterSampleOffset;

	const float thetaSquare;

	Spectrum unitSpectrum;
public:
	BlackBodyMultiScattering(float stepSize, int nxMipmap, int nyMipmap, int nzMipmap);
	virtual ~BlackBodyMultiScattering(void);

	virtual Spectrum Li(const Scene *scene, const Renderer *renderer,
		const RayDifferential &ray, const Sample *sample, RNG &rng,
		Spectrum *transmittance, MemoryArena &arena, MediaBidir *pMediaBidir = NULL) const;
	virtual Spectrum LiForNormalSingleScatterRay(const Scene *scene, const Renderer *renderer,
		const RayDifferential &ray, const Sample *sample, RNG &rng,
		Spectrum *transmittance, MemoryArena &arena, MediaBidir *pMediaBidir = NULL) const;
	virtual Spectrum Transmittance(const Scene *scene,
		const Renderer *renderer, const RayDifferential &ray,
		const Sample *sample, RNG &rng, MemoryArena &arena) const;
	void RequestSamples(Sampler *sampler, Sample *sample,
		const Scene *scene);

	inline TMipmap * GetMipmap(void) { return mipmap; }
	bool Debug_Intersect(Point &p, RayDifferential & _ray) const;

protected:
	Spectrum P_MS(float theta, const Spectrum & l) const;
	float P(float theta) const;
	Spectrum ComputeBlurWidth(float thetaSquare, const Spectrum l, float S, const Spectrum a, const Spectrum b, const Spectrum l2) const;
};

BlackBodyMultiScattering *CreateBlackBodyMultiScattering(const ParamSet &params);
