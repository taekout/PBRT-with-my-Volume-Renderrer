// integrators/single.cpp*
#include "stdafx.h"
#include "integrators/single.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"
#include "../lights/point.h"
#include "../thesis/HenyeyGreensteinSampleGenerator.h"
#include "../thesis/MediaBidir.h"
#include "../volumes/volumegrid.h"

/*
1. create no-sobolpoints(1024) in singlescattering class constructor.
2. figure out the pattern of sobol points.
   (Like : 0 - 100 : some similarly oriented vectors. 900 - 1024: the opposite direction vectors.
*/

SingleScatteringIntegrator::SingleScatteringIntegrator(float ss, int nDirectionsForEyeScatter, int nDirectionsForLightScatter, float _rayMarchingStepSize)
{
	this ->nDirectionsForEyeScatter = nDirectionsForEyeScatter;
	this ->nDirectionsForLightScatter = nDirectionsForLightScatter;
	this ->pHGScatterRayGenerator = new HenyeyGreensteinSampleGenerator(this ->nDirectionsForEyeScatter);
	this ->pHGLightRayGenerator = new HenyeyGreensteinSampleGenerator(this ->nDirectionsForLightScatter, 0.64f);
	stepSize = ss;
	rayMarchingStepSize = _rayMarchingStepSize;
	/*
	pMetroRenderer = NULL; pSobelGenerator = new SobelGenerator();
	pSobelGenerator ->DistributedVectorGenerator(nSobelSequenceOverSphere, "joe-kuo-other-4.5600", 2.f);
	RNG rng(3);
	this ->largeRayDirections = new vector<Vector *>;
	pSobelGenerator ->PickRandomPoints(rng, nSamplesToPickFromSobel, rayDirIndices, *largeRayDirections);

	// Generate Light Large Paths.
	pLightPathGenerator = new SobelGenerator();
	pLightPathGenerator ->DistributedVectorGenerator(1024, "joe-kuo-other-4.5600", 1.f);
	RNG rng2(4);
	this ->largeLightDirections  = new vector<Vector *>;
	pLightPathGenerator ->PickRandomPoints(rng2, 128, lightDirIndices, *largeLightDirections);
	// Now we have the light large paths. Both ray directions and light directions will be used in bidirectional path tracing function().
	*/
}

SingleScatteringIntegrator::~SingleScatteringIntegrator(void)
{
	delete pHGScatterRayGenerator;
	delete pHGLightRayGenerator;
}

// SingleScatteringIntegrator Method Definitions
void SingleScatteringIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}


Spectrum SingleScatteringIntegrator::Transmittance(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}


Spectrum SingleScatteringIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *T, MemoryArena &arena, MediaBidir *pMediaBidir) const {
    VolumeRegion *vr = scene->volumeRegion;
    float t0, t1;
	if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
        *T = 1.f;
        return 0.f;
    }
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);

    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize); // stepSize 
    float step = (t1 - t0) / nSamples; // one step distance (for real). Do not use stepSize, but use step. That is the real one.
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step; // Two oneD(=one dimensional) patterns are used for selecting which light to sample.
	// A two-dimensional pattern is used for selecting points on area light sources.
	float *lightNum = arena.Alloc<float>(nSamples);
	LDShuffleScrambled1D(1, nSamples, lightNum, rng);
	float *lightComp = arena.Alloc<float>(nSamples);
	LDShuffleScrambled1D(1, nSamples, lightComp, rng);
	float *lightPos = arena.Alloc<float>(2*nSamples);
	LDShuffleScrambled2D(1, nSamples, lightPos, rng);
	int sampOffset = 0;
    // Compute sample patterns for single scattering samples
    for (int i = 0; i < nSamples ; ++i, t0 += step) { // In the single scattering equ, it is the most out integral part. Ex> Integral of [0,t]
        // Advance to sample at _t0_ and update _T_
        pPrev = p; // save previous position before stepping forward(ray marching)
        p = ray(t0); // step forwarded position == p
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth); // What is tauRay? -> To get the transmittance, it is necessary to get the previous Ray and then compute the prev tau.
        Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau); // p.880  Tr(pi->p) = Tr(pi -> p_i-1) * Tr(p_i-1 -> p).

        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }

        // Compute single-scattering source term at _p_
		Lv += Tr * vr->Lve(p, w, ray.time); // This part is to compute emission-term at p. ==> In fire smoke, Emission = 0.
        Spectrum ss = vr->sigma_s(p, w, ray.time); // Compute sigma_s(sigma_s is out scattering coefficient. sigmal = sigma_a + sigma_s.)
        if (!ss.IsBlack() && scene->lights.size() > 0) {
            // Add contribution of _light_ due to scattering at _p_
                        // This part computes the equ. on p. 884
			// pick a closest light.
			float minDistance = 10000000.0f;
			int closestLightIndex = -1;
			int nLights = scene ->lights.size();
			for(int lightNum = 0 ; lightNum < nLights ; lightNum++) {
				PointLight * pLight = (PointLight *)scene ->lights.at(lightNum);
				Point lightPt = pLight ->GetPoint();
				float dist = Distance(p, lightPt);
				if(minDistance > dist) {
					minDistance = dist;
					closestLightIndex = lightNum;
				}
			}

			PointLight *light = (PointLight *)scene->lights[closestLightIndex]; // select only one light.(it should be okay as long as it takes enough number of samples along the ray
			Spectrum closestIntensity = light ->GetIntensity();

			float pdf = 0.f;
            VisibilityTester vis;
            Vector wo;
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis); // I don't know why. In case of point lights, pdf is just set to 1. --> I got it: It simply attenuates the power based on pdf.(like L /= pdf) so For PointLight, pdf = 1. For area light, pdf = some value based on outgoing angle and sold angle.
			Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);

			PointLight *sunLight = (PointLight *)scene->lights[0]; // light outside the box is placed at index 0.
			Spectrum sunL = sunLight->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
			Spectrum sunLd = sunL * vis.Transmittance(scene, renderer, NULL, rng, arena);

            Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * (Ld + sunLd) * float(nLights) /
                    pdf;
			
        }
		sampOffset++;
    }
    *T = Tr;
    return Lv * step;
}

Spectrum SingleScatteringIntegrator::LiForNormalSingleScatterRay(const Scene *scene, const Renderer *renderer,
	const RayDifferential &ray, const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena, MediaBidir *pMediaBidir) const
{
    VolumeRegion *vr = scene->volumeRegion;
    float t0 = -1.f, t1 = -1.f;
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);

	if(t0 == -1.f) { printf("You don't use LiForNormalSingleScatterRay.() Deprecated."); exit(-1); }
    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;

    // Compute sample patterns for single scattering samples
    float *lightNum = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
    float *lightComp = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
    float *lightPos = arena.Alloc<float>(2*nSamples);
    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
    uint32_t sampOffset = 0;
    for (int i = 0; i < nSamples ; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
        Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau);

        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }

        // Compute single-scattering source term at _p_
        Lv += Tr * vr->Lve(p, w, ray.time);
        Spectrum ss = vr->sigma_s(p, w, ray.time);
        if (!ss.IsBlack() && scene->lights.size() > 0) {
//            int nLights = scene->lights.size();
			int ln = scene ->nNormalLightIndex;
            Light *light = scene->lights[ln];
            // Add contribution of _light_ due to scattering at _p_
            float pdf;
            VisibilityTester vis;
            Vector wo;
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            
            if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {
                Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);
                Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * Ld /
                        pdf;
            }
        }
        ++sampOffset;
    }
    *T = Tr;
    return Lv * step;
}


SingleScatteringIntegrator *CreateSingleScatteringIntegrator(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
	//int maxRayMarching = params.FindOneInt("maxRayMarching", 7);
	int nDirectionsForLightScatter = params.FindOneInt("nLargeLightScatterRays", 1024);
	int nDirectionsForEyeScatter = params.FindOneInt("nLargeScatterRays", 16);
	float rayMarchingStepSize = params.FindOneFloat("rayMarchingStepSize", 0.25);
	return new SingleScatteringIntegrator(stepSize, nDirectionsForEyeScatter, nDirectionsForLightScatter, rayMarchingStepSize);
}


