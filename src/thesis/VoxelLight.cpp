#include "stdafx.h"
#include "pbrt.h"
#include "VoxelLight.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"
#include "BlackBodyVolume.h"
#include <time.h>

#define TAEKYUDEBUG_DEPRECATED_VOXELLIGHT_NOT_USED		0

// VoxelLight Method Definitions
VoxelLight::VoxelLight(const Transform &light2world,
	const Spectrum &intensity, const BBox & voxelExtent, RNG * _rng)
	: Light(light2world) {
		Point lightPos = LightToWorld(Point(0,0,0));
		extent = voxelExtent;
		Vector sizeOfExtent = voxelExtent.pMax - voxelExtent.pMin;
		midPoint = (voxelExtent.pMax + voxelExtent.pMin) * 0.5f;
#if TAEKYUDEBUG_DEPRECATED_VOXELLIGHT_NOT_USED
		if(midPoint.x != lightPos.x || midPoint.y != lightPos.y || midPoint.z != lightPos.z) {
			printf("OMG light points are diffferent.\n");
			exit(-1);
		}
#endif
		extentSize.x = sizeOfExtent.x;
		extentSize.y = sizeOfExtent.y;
		extentSize.z = sizeOfExtent.z;
		Intensity = intensity;

		this->rng = _rng;
#if TAEKYUDEBUG_DEPRECATED_VOXELLIGHT_NOT_USED
		if(Intensity.At(0) > 600.f) {
			printf("");
		}
#endif
}

VoxelLight::~VoxelLight() {
}

Spectrum VoxelLight::Sample_L(const Point &p, float pEpsilon,
	const LightSample &ls, float time, Vector *wi, float *pdf,
	VisibilityTester *visibility) const {
		Point tentativeLightPosition = GetCenterPoint();
		*wi = Normalize(tentativeLightPosition - p);
		*pdf = 1.f;
		visibility->SetSegment(p, pEpsilon, tentativeLightPosition, 0., time);
		return Intensity / DistanceSquared(tentativeLightPosition, p);
}


Spectrum VoxelLight::Power(const Scene *) const {
	return 4.f * M_PI * Intensity;
}

VoxelLight *CreateFireVoxelLights(VolumeRegion * vgd, const Transform &light2world, const ParamSet &paramSet,
	vector<BBox> & lightVoxels, vector<Color3d> &lightColor, vector<Light*> *multipleVoxelLights) {

	RNG *rng = new RNG(); // Memory Leak.
	Spectrum I;// = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.
	Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0)); // scale
	Point P = paramSet.FindOnePoint("from", Point(0,0,0)); // from --> Light location --> Only this needs to be multiply updated.
	Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world; // This has to be computed for every position.
	if(vgd == NULL) { printf("Volume Pointer Null. Get out.\n"); exit(-12); }
	multipleVoxelLights ->reserve(MAX_LIT_VOXELS);
	BlackBodyVolume * pvgd = dynamic_cast<BlackBodyVolume *>(vgd);
	if(pvgd == NULL) {
		printf("It must be VolumeGridDensity type of a volume.\n"); exit(-11);
	}
	//add fire light ==> I.
	int nx = -1, ny = -1, nz = -1;
	float sizex = -1.f, sizey = -1.f, sizez= -1.f;
	pvgd->GetDimension(nx, ny, nz);
	pvgd->GetExtentSize(sizex, sizey, sizez);
	sizex /= nx; sizey /= ny; sizez /= nz;
	vector<Color3d>::iterator colorIt = lightColor.begin();
	for(vector<BBox>::iterator it = lightVoxels.begin() ; it < lightVoxels.end() ; it++, colorIt++) {
		BBox box = (*it);
		Color3d c = *colorIt;
		float fc[3] = {c.r, c.g, c.b};
		I = Spectrum::FromRGB(fc);
		Point midP = (box.pMin + box.pMax) / 2.0f;
		Transform l2w = Translate(Vector(midP.x, midP.y, midP.z)) * light2world;
		if( I.At(0) < 0.0000001 && I.At(1) < 0.0000001 && I.At(2) < 0.0000001 )
			;
		else
			multipleVoxelLights ->push_back(new VoxelLight(l2w, I, box, rng)); //dynamically created lights have also the same power as indicated.

	}

	I = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.

	Point sunLightMax(P.x + sizex / (float)nx / 2.f, P.y + sizey / (float)ny / 2.f, P.z + sizez / (float)nz / 2.f);
	Point sunLightMin(P.x - sizex / nx / 2.f, P.y - sizey / ny / 2.f, P.z - sizez / nz / 2.f);

	BBox sunLightBox(sunLightMin, sunLightMax);
	return new VoxelLight(l2w, I * sc, sunLightBox, rng);
}

VoxelLight *CreateVoxelLight(const Transform &light2world,
	const ParamSet &paramSet, vector<BBox> &lightVoxels, vector<Light*> *multipleVoxelLights) {
		Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.
		Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0)); // scale
		Point P = paramSet.FindOnePoint("from", Point(0,0,0)); // from --> Light location --> Only this needs to be multiply updated.
		Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world; // This has to be computed for every position.
		multipleVoxelLights ->reserve(MAX_LIT_VOXELS);
		RNG *rng = new RNG(); // Memory Leak
		for(vector<BBox>::iterator it = lightVoxels.begin() ; it < lightVoxels.end() ; it++) {
			BBox box = (*it);
			Point p = (box.pMin + box.pMax) / 2.0f;
			Transform l2w = Translate(Vector(p.x, p.y, p.z)) * light2world;
			multipleVoxelLights ->push_back(new VoxelLight(l2w, I * sc, box, rng)); //dynamically created lights have also the same power as indicated.
		}

		Point sunLightMax(P.x + 0.05f, P.y + 0.05f, P.z + 0.05f);
		Point sunLightMin(P.x - 0.05f, P.y - 0.05f, P.z - 0.05f);

		BBox sunLightBox(sunLightMin, sunLightMax);
		return new VoxelLight(l2w, I * sc, sunLightBox, rng);
}


float VoxelLight::Pdf(const Point &, const Vector &) const {
	return 0.;
}


Spectrum VoxelLight::Sample_L(const Scene *scene, const LightSample &ls,
	float u1, float u2, float time, Ray *ray, Normal *Ns,
	float *pdf) const {
		printf("Deprecated : Sample_L() function.\n");
		exit(-1);
		return Spectrum(-1.f);
		/*
		*ray = Ray(GetPoint(), UniformSampleSphere(ls.uPos[0], ls.uPos[1]),
			0.f, INFINITY, time);
		*Ns = (Normal)ray->d;
		*pdf = UniformSpherePdf();
		return Intensity;
		*/
}


void VoxelLight::SHProject(const Point &p, float pEpsilon, int lmax,
	const Scene *scene, bool computeLightVisibility, float time,
	RNG &rng, Spectrum *coeffs) const {
	printf("VoxelLight SHPRoject has been deprecated.\n");
	exit(-30);
}

Point VoxelLight::GetPoint(void) const
{
#if TAEKYUDEBUG_DEPRECATED_VOXELLIGHT_NOT_USED
	return midPoint;
#endif
#if !TAEKYUDEBUG_DEPRECATED_VOXELLIGHT_NOT_USED
	// pick a random position in the voxel.
	Point p = extentSize * rng->RandomFloat();
	p += extent.pMin;

	return p;
#endif
}

Point VoxelLight::GetCenterPoint(void) const
{
	return midPoint;
}