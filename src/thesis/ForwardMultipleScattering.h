#if defined(_MSC_VER)
#pragma once
#endif
#include "volume.h"
#include "integrator.h"

class HenyeyGreensteinSampleGenerator;
class MediaBidir;

// ForwardMultipleScattering Declarations
class ForwardMultipleScattering : public VolumeIntegrator {
public:
	// ForwardMultipleScattering Public Methods
	ForwardMultipleScattering(float ss, int nMultipleScatterSamples, bool bVariableMultScatter, float HRRPUVUnit, int nMaxScatter, float fHowForward = 0.64f);
	virtual ~ForwardMultipleScattering(void);
	Spectrum Transmittance(const Scene *, const Renderer *,
		const RayDifferential &ray, const Sample *sample, RNG &rng,
		MemoryArena &arena) const;
	void RequestSamples(Sampler *sampler, Sample *sample,
		const Scene *scene);
	Spectrum Li(const Scene *, const Renderer *, const RayDifferential &ray,
		const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena, MediaBidir *pMediaBidir = NULL) const;

	inline Spectrum * LightTransport(const VolumeRegion * vr, const Light * light, int lightIndex, const Point & p, const LightSample & ls, const RayDifferential &ray, Vector & wo, float & pdf, VisibilityTester &vis, const Scene * scene, const Renderer * renderer, RNG & rng, MemoryArena & arena ) const;

	void GetClosestLight( int nLights, const Scene * scene, Point p, float &minDistance, int &closestLightIndex ) const;

	Spectrum LiForScattering(const Scene *, const Renderer *, const RayDifferential &ray,
		const Sample *sample, RNG &rng) const;
	Spectrum LiForNormalSingleScatterRay(const Scene *scene, const Renderer *renderer,
		const RayDifferential &ray, const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena, MediaBidir *pMediaBidir) const;

private:
	// ForwardMultipleScattering Private Data
	float stepSize;
	int tauSampleOffset, scatterSampleOffset;

	// Thesis by Taekyu Shin
	int nMultipleScatterSamples;
	float fHowForward;
	HenyeyGreensteinSampleGenerator *pHGMultipleScatteringRayGenerator;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator;

	Mutex * mutex;
	vector<Vector *> gIsoDir;

	// var mult scattering
	bool bVariableMultScatter;
	float HRRPUVUnit;
	int nMaxScatter;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator2;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator4;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator6;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator8;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator10;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator12;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator14;
	HenyeyGreensteinSampleGenerator *pHGIsotropicScatteringRayGenerator16;
};


ForwardMultipleScattering *CreateForwardMultipleScattering(const ParamSet &params);
