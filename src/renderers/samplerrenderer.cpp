// renderers/samplerrenderer.cpp*
#include "stdafx.h"
#include "renderers/samplerrenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "metropolis.h"
#include "../thesis/MediaBidir.h"
#include "../thesis/MultipleScatteringIntegrator.h"
#include "../thesis/BlackBodyMultiScattering.h"
#include "../thesis/BlackBodyVolume.h"
#include "../thesis/Preprocessor.h"

static uint32_t hash_myown(char *key, uint32_t len)
{
    uint32_t   hash, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
} 

// SamplerRendererTask Definitions
void SamplerRendererTask::Run() { // In Direct Rendering, this will be run. In case of metropolis? --> MLTTask::Run()
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            if (visualizeObjectIds) { // I guess this is the intersection case. For example, volume rendering does not need intersection? so visualizeObjectIds == false
                if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
                    // random shading based on shape id...
                    uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                    uint32_t h = hash_myown((char *)ids, sizeof(ids));
                    float rgb[3] = { (h & 0xff), (h >> 8) & 0xff, (h >> 16) & 0xff };
                    Ls[i] = Spectrum::FromRGB(rgb);
                    Ls[i] /= 255.f;
                }
                else
                    Ls[i] = 0.f;
            }
            else { // volume Rendering
            if (rayWeight > 0.f)
                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i]);
            else {
                Ls[i] = 0.f;
                Ts[i] = 1.f;
            }

            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                printf("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
				exit(-80);
            }
            else if (Ls[i].y() < -1e-5) {
                printf("Negative luminance value, %f, returned"
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
				exit(-81);
            }
            else if (isinf(Ls[i].y())) {
                printf("Infinite luminance value returned"
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
				exit(-82);
            }
			/*else if(Ls[i].At(0) == 0.0f || Ls[i].At(1 ) == 0.0f || Ls[i].At(2) == 0.0f) {
				printf("Something to look into here.\n");
			}*/
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// SamplerRenderer Method Definitions
SamplerRenderer::SamplerRenderer(Sampler *s, Camera *c,
                                 SurfaceIntegrator *si, VolumeIntegrator *vi,
                                 bool visIds, MediaBidir * mediaRenderer) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
	this ->mediaRenderer = mediaRenderer;
}


SamplerRenderer::~SamplerRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}

Point gStratifiedSamples[SAMPLE_CASES][STRATIFIEDSAMPLENUMBER];

void StraitifiedSample3D(void) {
	const int _nx = STRATIFIEDSAMPLENUMBER_X, _ny = STRATIFIEDSAMPLENUMBER_Y, _nz = STRATIFIEDSAMPLENUMBER_Z; // STRATIFIEDSAMPLENUMBER = 10 x 10 x 10.
	static RNG rng;
	float dx = 1.f / _nx, dy = 1.f / _ny, dz = 1.f / _nz;
	for(int sampleID = 0 ; sampleID < SAMPLE_CASES ; sampleID++) {
		for(int z = 0 ; z < _nz ; z++) {
			for (int y = 0; y < _ny; ++y) {
				for (int x = 0; x < _nx; ++x) {
					float jx = rng.RandomFloat(); // make sure RandomFloat() returns 0 <= x <= 1.0.
					float jy = rng.RandomFloat();
					float jz = rng.RandomFloat();
					gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].x = min((x + jx) * dx, OneMinusEpsilon) - 0.5f;
					gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].y = min((y + jy) * dy, OneMinusEpsilon) - 0.5f;
					gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].z = min((z + jz) * dz, OneMinusEpsilon) - 0.5f;
					if( gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].HasNaNs() ) {
						printf("NaN Point.\n");
						exit(-76);
					}
					else if( gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].x >= 0.5f ||
						gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].y >= 0.5f ||
						gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].z >= 0.5f)
					{
						printf("Oh no...");
						exit(-54);
					}
					else if( gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].x <= -0.5f ||
						gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].y <= -0.5f ||
						gStratifiedSamples[sampleID][z * _ny * _nx + y * _nx + x].z <= -0.5f)
					{
						printf("Oh no...2");
						exit(-54);
					}
				}
			}
		}
	}

}


void SamplerRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
	if(dynamic_cast<PremozeMultipleScatteringIntegrator *>(volumeIntegrator) != NULL) {
		//Preprocess.... 
		PremozeMultipleScatteringIntegrator * pIntegrator = dynamic_cast<PremozeMultipleScatteringIntegrator *>(volumeIntegrator);
		pIntegrator->mScene = (Scene *)scene;
		if(pIntegrator != NULL) {
			if(scene->lights.size() > 0) {
				if(dynamic_cast<VolumeGridDensity *>(scene->volumeRegion) != NULL) {
					pIntegrator->mPreprocessor = new Preprocessor((VolumeGridDensity *)scene->volumeRegion, pIntegrator->stepSize, true);
					pIntegrator->mPreprocessor->PreprocessVolume(pIntegrator->GetMipmap(), (vector<Light *> &) scene->lights);
				}
				else {
					printf("Volume Data is not the type of VolumeGridDensity.\n");
					exit(-50);
				}
			}
		}
	}
	else if(dynamic_cast<BlackBodyMultiScattering *>(volumeIntegrator) != NULL) {
		//Preprocess.... 
		BlackBodyMultiScattering * pIntegrator = dynamic_cast<BlackBodyMultiScattering *>(volumeIntegrator);
		pIntegrator->mScene = (Scene *)scene;
		if(pIntegrator != NULL) {
			if(scene->lights.size() > 0) {
				if(dynamic_cast<BlackBodyVolume *>(scene->volumeRegion) != NULL) {
					pIntegrator->mPreprocessor = new Preprocessor((VolumeGridDensity *)scene->volumeRegion, pIntegrator->stepSize, true);
					pIntegrator->mPreprocessor->PreprocessVolume(pIntegrator->GetMipmap(), (vector<Light *> &) scene->lights);
				}
				else {
					printf("Volume Data is not the type of VolumeGridDensity.\n");
					exit(-50);
				}
			}
		}
	}

	// Here, I initialize point light 
	{

		RNG rng;
		for(int i = 0 ; i < scene->lights.size() ; i++) {
			Light * pLight = scene->lights.at(i);
			PointLight * pPointLight = dynamic_cast<PointLight *> (pLight);
			if(pPointLight == NULL) {
				printf("PointLights are allowed.\n");
				exit(-98);
			}
			int ID = (int)rng.RandomFloat() * 5.0;
			if(ID < 0) ID = 0;
			if(ID >= 5) ID = 4;
			pPointLight->SetSampleID(ID);
			if(ID == -1) { printf("Wrong ID"); exit(-33); }
		}
		StraitifiedSample3D();
	}
    
	PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);

    // Create and launch _SamplerRendererTask_s for rendering image

    // Compute number of _SamplerRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16)); // ?
    nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Rendering");
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i) // I guess here they divide and assign each task to threads.?
        renderTasks.push_back(new SamplerRendererTask(scene, this, camera,
                                                      reporter, sampler, sample, 
                                                      visualizeObjectIds, 
                                                      nTasks-1-i, nTasks));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();

    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();
    PBRT_FINISHED_RENDERING();
	PremozeMultipleScatteringIntegrator * pMultIntegrator = dynamic_cast<PremozeMultipleScatteringIntegrator *>(volumeIntegrator);
	if(pMultIntegrator != NULL) {
		printf("\n\n%d times : %f = %f / %d\n", pMultIntegrator->mPreprocessor->nPhaseEvents,
			pMultIntegrator->mPreprocessor->sumN / pMultIntegrator->mPreprocessor->nPhaseEvents,
			pMultIntegrator->mPreprocessor->sumN, pMultIntegrator->mPreprocessor->nPhaseEvents);
	}
	BlackBodyMultiScattering * pBlackIntegrator = dynamic_cast<BlackBodyMultiScattering *>(volumeIntegrator);
	if(pBlackIntegrator != NULL) {
		printf("\n\n%d times : %f = %f / %d\n", pBlackIntegrator->mPreprocessor->nPhaseEvents,
			pBlackIntegrator->mPreprocessor->sumN / pBlackIntegrator->mPreprocessor->nPhaseEvents,
			pBlackIntegrator->mPreprocessor->sumN, pBlackIntegrator->mPreprocessor->nPhaseEvents);
	}
    // Clean up after rendering and store final image
    delete sample;
    camera->film->WriteImage();

	//_CrtDumpMemoryLeaks();

}


Spectrum SamplerRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect))
		// TEsted surfaceIntegrator => surfaceIntegrator->Li() results in the material color + direct light shading on the surface.
		Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
									rng, arena);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Li += scene->lights[i]->Le(ray); // 0.f for PointLight
    }
        // If volumeIntegrator is EmissionIntegrator, then it computes 2nd term in equation (16.2)
	// Tested volumeIntegrator => volumeIntegrator->Li() results in the volumetric lighting in the air.
	Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng, T, arena, this ->mediaRenderer);
//	if(scene ->nNormalLightIndex != -1)
//		Lvi += volumeIntegrator ->LiForNormalSingleScatterRay(scene, this, ray, sample, rng, T, arena, this ->mediaRenderer);
    return *T * Li + Lvi;
}


Spectrum SamplerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


