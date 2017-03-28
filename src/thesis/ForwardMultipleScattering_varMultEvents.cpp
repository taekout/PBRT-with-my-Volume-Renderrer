// integrators/single.cpp*
#include "stdafx.h"
#include "ForwardMultipleScattering.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"
#include "../lights/point.h"
#include "../thesis/VoxelLight.h"
#include "../thesis/HenyeyGreensteinSampleGenerator.h"
#include "../thesis/MediaBidir.h"
#include "../volumes/volumegrid.h"
#include "../thesis/BlackBodyVolume.h"
#include "../lights/diffuse.h"

#define SUNLIGHT	0
#define TAEKYUDEBUG	0

#if SUNLIGHT
extern Light * gSunLight;
#endif

#if TAEKYUDEBUG
void SanityCheck(const Vector & vec) {
	
	if( !IsCloseWithEpisilon(vec.LengthSquared(), 1.0f) ) {
		printf("OMG. SanityCheck was not passed. The vector length is not 1.0.\n");
		exit(-15);
	}
	else if( vec.HasNaNs() ) {
		printf("OMG. SanityCheck was not passed. NaN value for the vector here.\n");
		exit(-17);
	}
}
#endif

ForwardMultipleScattering::ForwardMultipleScattering(float ss, int _nMultipleScatterSamples, float _fHowForward)
{

	nMultipleScatterSamples = _nMultipleScatterSamples;
	if( nMultipleScatterSamples == 0 ) {
		this->pHGMultipleScatteringRayGenerator2 = NULL;
		this->pHGMultipleScatteringRayGenerator4 = NULL;
		this->pHGMultipleScatteringRayGenerator6 = NULL;
		this->pHGMultipleScatteringRayGenerator8 = NULL;
		this->pHGMultipleScatteringRayGenerator10 = NULL;
		this->pHGIsotropicScatteringRayGenerator = NULL;
	}
	else {
		this->pHGMultipleScatteringRayGenerator2 = new HenyeyGreensteinSampleGenerator(2, _fHowForward);
		this->pHGMultipleScatteringRayGenerator4 = new HenyeyGreensteinSampleGenerator(4, _fHowForward);
		this->pHGMultipleScatteringRayGenerator6 = new HenyeyGreensteinSampleGenerator(6, _fHowForward);
		this->pHGMultipleScatteringRayGenerator8 = new HenyeyGreensteinSampleGenerator(8, _fHowForward);
		this->pHGMultipleScatteringRayGenerator10 = new HenyeyGreensteinSampleGenerator(10, _fHowForward);
		this->pHGIsotropicScatteringRayGenerator = new HenyeyGreensteinSampleGenerator(this->nMultipleScatterSamples, 0.01f);	// Always Isotropic scattering.
	}
	this->fHowForward = _fHowForward;
	stepSize = ss;
}

ForwardMultipleScattering::~ForwardMultipleScattering(void)
{
	delete pHGMultipleScatteringRayGenerator2;
	delete pHGMultipleScatteringRayGenerator4;
	delete pHGMultipleScatteringRayGenerator6;
	delete pHGMultipleScatteringRayGenerator8;
	delete pHGMultipleScatteringRayGenerator10;
	delete pHGIsotropicScatteringRayGenerator;
	pHGMultipleScatteringRayGenerator2 = pHGMultipleScatteringRayGenerator4 = pHGMultipleScatteringRayGenerator6 = pHGMultipleScatteringRayGenerator8 = pHGMultipleScatteringRayGenerator10 = NULL;
	pHGIsotropicScatteringRayGenerator = NULL;
}

// ForwardMultipleScattering Method Definitions
void ForwardMultipleScattering::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}


Spectrum ForwardMultipleScattering::Transmittance(const Scene *scene,
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

Spectrum * ForwardMultipleScattering::LightTransport(const VolumeRegion * vr, const Light * light, int lightIndex, const Point & p, const LightSample & ls, const RayDifferential &ray, Vector & wo, float & pdf, VisibilityTester &vis, const Scene * scene, const Renderer * renderer, RNG & rng, MemoryArena & arena ) const
{
	// Check Voxel b/f sampling lights.
	BlackBodyVolume *bbv = ((BlackBodyVolume *)vr);
	Spectrum *Ld = arena.Alloc<Spectrum> ();
	bool bVoxelizationSuccess = false;
	if(bbv->bStoreDataInVoxel && lightIndex != -1) { // If Sun light, lightIndex == -1.
		int _x, _y, _z;
		bbv->GetVoxelIndex(p, _x, _y, _z);
		pdf = 1.f;
		bVoxelizationSuccess = bbv->GetPowerStoredInVoxel(lightIndex, _x, _y, _z, *Ld);
		if(bVoxelizationSuccess == false) {
			//Store
			Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
			*Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);
			bbv->StorePowerInVoxel(lightIndex, _x, _y, _z, *Ld);
			return Ld;
		}
		else
			return Ld;
	}

	// sampling lights.
	Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
	*Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);

	return Ld;
}

bool g_bClosestLight = true;

Spectrum ForwardMultipleScattering::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *T, MemoryArena &arena, MediaBidir *pMediaBidir) const {
    VolumeRegion *vr = scene->volumeRegion;
	BlackBodyVolume *bbv = NULL;
	if( !(bbv = dynamic_cast<BlackBodyVolume *> (vr)) ) {
		printf("Only blackbody volume is accepted.\n");
		exit(-91);
	}

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
		if(!bbv->IsInsideExtent(p)) continue;
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
			
			int lightIndex = -1;
			int nLights = scene ->lights.size();
			float minDistance = 10000000.0f;

			if(g_bClosestLight)
				GetClosestLight(nLights, scene, p, minDistance, lightIndex);
			else {
				lightIndex = min(Floor2Int(lightNum[sampOffset] * nLights), nLights-1);
			}

			Light *light;
			light = (PointLight *) scene->lights[lightIndex]; // select only one light.(it should be okay as long as it takes enough number of samples along the ray
			//Spectrum closestIntensity = light ->GetIntensity();

			float pdf = 0.f;
            VisibilityTester vis;
            Vector wo;
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
// I don't know why. In case of point lights, pdf is just set to 1. --> I got it: It simply attenuates the power based on pdf.(like L /= pdf) so For PointLight, pdf = 1. For area light, pdf = some value based on outgoing angle and sold angle.
			Spectrum *Ld = LightTransport(vr, light, lightIndex, p, ls, ray, wo, pdf, vis, scene, renderer, rng, arena);


			// <Multiple Scattering Steps>
			// 1. HG Creates N samples. 
			Spectrum scatterSpectrum(0.f);
			int nScatters = 0;
			float HRRPUVgradient = 0.0f;
			Vector & pivotVectorForMult = ((BlackBodyVolume *)vr)->GetGradientVector(p, HRRPUVgradient);
			int nQuotient = HRRPUVgradient / 100;
			nQuotient = (nQuotient > 10) ? 10 : nQuotient;
			nScatters = nQuotient; //multScatterSamples.size();
			nScatters /= 2; nScatters *= 2;
			if(pHGMultipleScatteringRayGenerator2 != NULL && nScatters != 0) {
				std::vector<Vector *> multScatterSamples;
				float lengthSquared = pivotVectorForMult.LengthSquared();
				
				if(lengthSquared == 0.0f)
					pHGIsotropicScatteringRayGenerator->TransformDirectionIntoCurrentRaySpace( ray.d, multScatterSamples );
				else {
#if TAEKYUDEBUG
					SanityCheck(pivotVectorForMult);
#endif
					
					switch(nScatters) {
					case 2:
						pHGMultipleScatteringRayGenerator2->TransformDirectionIntoCurrentRaySpace( ray.d, multScatterSamples );
						break;
					case 4:
						pHGMultipleScatteringRayGenerator4->TransformDirectionIntoCurrentRaySpace( ray.d, multScatterSamples );
						break;
					case 6:
						pHGMultipleScatteringRayGenerator6->TransformDirectionIntoCurrentRaySpace( ray.d, multScatterSamples );
						break;
					case 8:
						pHGMultipleScatteringRayGenerator8->TransformDirectionIntoCurrentRaySpace( ray.d, multScatterSamples );
						break;
					case 10:
						pHGMultipleScatteringRayGenerator10->TransformDirectionIntoCurrentRaySpace( ray.d, multScatterSamples );
						break;
					default:
						printf("Wrong mult ray numbers. HRRPUV gradient\n");
						exit(-93);
					}
					
					//pHGMultipleScatteringRayGenerator->TransformDirectionIntoCurrentRaySpace( pivotVectorForMult, multScatterSamples );
				}
				// 2. Call LiMultipleScattering(ray);
				

				if(nScatters != 0) {

					float pdfSum = 0.0f;
					for(int k = 0 ; k < nScatters ; k++) {
						//Ray multray((const Point &) p, *multScatterSamples.at(k), 0.0f);
						RayDifferential multray((const Point &) p, *multScatterSamples.at(k), 0.0f);
						multray.d.Normalize();
						//Ray tempEyeRay = ray;
						Spectrum inScatterL;
#if TAEKYUDEBUG
						if(!IsCloseWithEpisilon(multray.d.Length(), 1.f) ) {
							printf("Interestingly, multray.d is not normalized.:%f", multray.d.Length());
							exit(-11);
						}
						if(!IsCloseWithEpisilon(ray.d.Length(), 1.f)) {
							printf("Danger! Not Normalized. : %f", ray.d.Length());
							exit(-75);
						}
#endif
						float dotp = //max(Dot( multray.d, ray.d ), 0.0f);
							PhaseHG(multray.d, ray.d, fHowForward);
						if(dotp != 0)
							inScatterL = LiForScattering( scene, renderer, multray, sample, rng, arena ) * dotp;
#if TAEKYUDEBUG
						SanityCheck(inScatterL);
#endif
						if(lengthSquared == 0.0f) {
							float pdfForOneMult = pHGIsotropicScatteringRayGenerator->GetPDF(k);
							scatterSpectrum += inScatterL;// * pdfForOneMult;
							pdfSum += pdfForOneMult;
						}
						else {
							float pdfForOneMult;
							switch(nScatters) {
							case 2:
								pdfForOneMult = pHGMultipleScatteringRayGenerator2->GetPDF(k);
								break;
							case 4:
								pdfForOneMult = pHGMultipleScatteringRayGenerator4->GetPDF(k);
								break;
							case 6:
								pdfForOneMult = pHGMultipleScatteringRayGenerator6->GetPDF(k);
								break;
							case 8:
								pdfForOneMult = pHGMultipleScatteringRayGenerator8->GetPDF(k);
								break;
							case 10:
								pdfForOneMult = pHGMultipleScatteringRayGenerator10->GetPDF(k);
								break;
							default:
								printf("Wrong mult ray numbers. HRRPUV gradient\n");
								exit(-93);
							}
							//float pdfForOneMult = pHGMultipleScatteringRayGenerator->GetPDF(k);
							scatterSpectrum += inScatterL * pdfForOneMult;
							pdfSum += pdfForOneMult;
						}

					}
					scatterSpectrum = scatterSpectrum * (4 * 3.141592 /*/ pdfSum*/ / nScatters);
				}
				
				for(int k = 0 ; k < nScatters ; k++) {
					delete multScatterSamples.at(k);
				}
				multScatterSamples.clear();	// Does it crash here?
			}

			
#if SUNLIGHT
			PointLight *sunLight = (PointLight *)gSunLight;
			Spectrum *sunLd = LightTransport(vr, sunLight, -1, p, ls, ray, wo, pdf, vis, scene, renderer, rng, arena);
			//Spectrum *Ld = LightTransport(light, p, ls, ray, wo, pdf, vis, scene, renderer, rng, arena);
			//Spectrum sunL = sunLight->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
			//Spectrum sunLd = sunL * vis.Transmittance(scene, renderer, NULL, rng, arena);
#endif
#ifdef TAEKYUDEBUG
			Spectrum tempLv = Lv;
#endif
            Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * (*Ld +
#if SUNLIGHT
				*sunLd +
#endif
				scatterSpectrum) * float(nLights) /
                    pdf;
#ifdef TAEKYUDEBUG
			if( Lv.At(0) < tempLv.At(0) || Lv.At(1) < tempLv.At(1) || Lv.At(2) < tempLv.At(2) ) {
				printf("Know this thing happens. After summation, the value decreased.\n");
				exit(-88);
			}
#endif
			
        }
		sampOffset++;
    }
    *T = Tr;

    return Lv * step;
}

/*
Ignore T (Transmission). T is used to attenuate the surface material shading. But it does not count for scatter rays since it does not affect material shading.
*/
Spectrum ForwardMultipleScattering::LiForScattering(const Scene *scene, const Renderer *renderer,
	const RayDifferential &ray, const Sample * sample, RNG &rng, MemoryArena &arena) const
{
	VolumeRegion *vr = scene->volumeRegion;
	float t0, t1;
	if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
		return 0.f;
	}
	t0 = (t0 < 0.f) ? 0.f : t0;
	if(t1 <= t0)
		return 0.f;
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
		if( ! ((BlackBodyVolume *) vr)->IsInsideExtent(p) ) continue;
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
			
			int lightIndex = -1;
			int nLights = scene ->lights.size();
			float minDistance = 10000000.0f;
			if(g_bClosestLight)
				GetClosestLight(nLights, scene, p, minDistance, lightIndex);
			else {
				lightIndex = min(Floor2Int(lightNum[sampOffset] * nLights), nLights-1);
			}

			Light *light;
			light = (PointLight *) scene->lights[lightIndex]; // select only one light.(it should be okay as long as it takes enough number of samples along the ray

			//Spectrum closestIntensity = light ->GetIntensity();

			float pdf = 0.f;
			VisibilityTester vis;
			Vector wo;
			LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
				lightPos[2*sampOffset+1]);
// I don't know why. In case of point lights, pdf is just set to 1. --> I got it: It simply attenuates the power based on pdf.(like L /= pdf) so For PointLight, pdf = 1. For area light, pdf = some value based on outgoing angle and sold angle.
			Spectrum *Ld = LightTransport(vr, light, lightIndex, p, ls, ray, wo, pdf, vis, scene, renderer, rng, arena);
			/*Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis); 
			Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);*/
#if SUNLIGHT
			PointLight *sunLight = (PointLight *)gSunLight;
			Spectrum *sunLd = LightTransport(vr, sunLight, -1, p, ls, ray, wo, pdf, vis, scene, renderer, rng, arena);
			/*Spectrum sunL = sunLight->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
			Spectrum sunLd = sunL * vis.Transmittance(scene, renderer, NULL, rng, arena);*/
#endif

			Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * (*Ld
#if SUNLIGHT
				+ *sunLd
#endif
				) * float(nLights) /
				pdf;
#if TAEKYUDEBUG
			SanityCheck(Lv);
#endif
		}
		sampOffset++;
	}
#if TAEKYUDEBUG
	SanityCheck(Lv);
	Spectrum test = Lv * step;
	SanityCheck(test);
#endif
	
	return Lv * step;
}

Spectrum ForwardMultipleScattering::LiForNormalSingleScatterRay(const Scene *scene, const Renderer *renderer,
	const RayDifferential &ray, const Sample *sample, RNG &rng, Spectrum *T, MemoryArena &arena, MediaBidir *pMediaBidir) const
{
	return Spectrum(.0f);
}

ForwardMultipleScattering *CreateForwardMultipleScattering(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
	//int maxRayMarching = params.FindOneInt("maxRayMarching", 7);
	int nMultipleScatterSamples = params.FindOneInt("multiplescatterSamples", 5);
	if(nMultipleScatterSamples > 0) {
		printf("Multiple Scattering option detected.\n");
	}
	float fForwardScattering = params.FindOneFloat("howforward", 0.01f);
	g_bClosestLight = params.FindOneBool("closestLight", true);
	return new ForwardMultipleScattering(stepSize, nMultipleScatterSamples, fForwardScattering);
}


void ForwardMultipleScattering::GetClosestLight( int nLights, const Scene * scene, Point p, float &minDistance, int &closestLightIndex ) const
{
	minDistance = 10000000.0f;
	for(int lightNum = 0 ; lightNum < nLights ; lightNum++) {
		PointLight * pLight = (PointLight *)scene ->lights.at(lightNum);
		Point lightPt = pLight ->GetPoint();
		float dist = Distance(p, lightPt);
		if(minDistance > dist) {
			minDistance = dist;
			closestLightIndex = lightNum;
		}
	}
}
