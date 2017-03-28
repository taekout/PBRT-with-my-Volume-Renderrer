#include "stdafx.h"
#include "lights/point.h"
#include "diffuse.h"
#include "./shapes/sphere.h"
#include "sh.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"
#include "./thesis/BlackBodyVolume.h"

#ifndef NDEBUG
#define TAEKYUDEBUG		1
#endif

// PointLight Method Definitions
PointLight::PointLight(const Transform &light2world,
                       const Spectrum &intensity, float _backGroundI, int _bVoxelLight)
    : Light(light2world) {
    lightPos = LightToWorld(Point(0,0,0));
    Intensity = intensity;
	bVoxelLight = _bVoxelLight;
	nStratifiedSampleID = -1;
	backgroundI = _backGroundI;
#if TAEKYUDEBUG
	if(Intensity.At(0) > 600.f) {
		printf("");
	}
#endif
}

float fVoxelSize = 0.024444444f;
void NaiveRandomNewPoint(const Point & p, Point & outNewPoint)
{
	Point newLightPos = p;
	float tmp = ((float)rand()) / (float)RAND_MAX;
	tmp -= 0.5;
	tmp *= fVoxelSize;
	newLightPos.x += tmp;

	tmp = ((float)rand()) / (float)RAND_MAX;
	tmp -= 0.5;
	tmp *= fVoxelSize;
	newLightPos.y += tmp;

	tmp = ((float)rand()) / (float)RAND_MAX;
	tmp -= 0.5;
	tmp *= fVoxelSize;
	newLightPos.z += tmp;
	outNewPoint = newLightPos;
}

extern Point gStratifiedSamples[SAMPLE_CASES][STRATIFIEDSAMPLENUMBER];

Spectrum PointLight::Sample_L(const Point &p, float pEpsilon,
         const LightSample &ls, float time, Vector *wi, float *pdf,
         VisibilityTester *visibility) const {
#if TAEKYUDEBUG
	if(nStratifiedSampleID == -1) {
		printf("Wrong ID in Sample_L function");
		exit(-11);
	}
#endif

	if(bVoxelLight == 0) { // Point Light
		*wi = Normalize(lightPos - p);
		*pdf = 1.f;
		visibility->SetSegment(p, pEpsilon, lightPos, 0., time);
		return Intensity / DistanceSquared(lightPos, p);
	}
	else if(bVoxelLight == 1) { // Voxel Light
		int index = ls.uPos[0] * (float)(STRATIFIEDSAMPLENUMBER - 1);
		Point newLightPos = gStratifiedSamples[nStratifiedSampleID][index] * fVoxelSize;
		newLightPos += lightPos;
		//newLightPos = Point(0.,0.,0.);
		//NaiveRandomNewPoint(lightPos, newLightPos);	//uncomment this if you want just naive random sampling.

		*wi = Normalize(newLightPos - p);
		*pdf = 1.f;
		visibility->SetSegment(p, pEpsilon, newLightPos, 0., time);
		return Intensity / DistanceSquared(newLightPos, p);
	}
	else if(bVoxelLight == 2) { // Spherical Light
		Point newLightPos;
		do {
			newLightPos = lightPos;
			float tmp = ((float)rand()) / (float)RAND_MAX;
			tmp -= 0.5;
			tmp *= fVoxelSize;
			newLightPos.x += tmp;

			tmp = ((float)rand()) / (float)RAND_MAX;
			tmp -= 0.5;
			tmp *= fVoxelSize;
			newLightPos.y += tmp;

			tmp = ((float)rand()) / (float)RAND_MAX;
			tmp -= 0.5;
			tmp *= fVoxelSize;
			newLightPos.z += tmp;
		} while( (newLightPos - lightPos).Length() > fVoxelSize );

		*wi = Normalize(newLightPos - p);
		*pdf = 1.f;
		visibility->SetSegment(p, pEpsilon, newLightPos, 0., time);
		return Intensity / DistanceSquared(newLightPos, p);
	}
	else {
		printf("light kind not exist.");
		exit(-13);
	}
}

Spectrum PointLight::Le(const RayDifferential &) const {
	return Spectrum(backgroundI);
}

Spectrum PointLight::Power(const Scene *) const {
    return 4.f * M_PI * Intensity;
}

PointLight *CreateFirePointLights(VolumeRegion *vgd, const Transform &light2world, const ParamSet &paramSet,
	vector<Point> &lightPos, vector<Color3d> &lightColor, vector<Light*> *multiplePointLights) {
    Spectrum I;// = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0)); // scale
    Point P = paramSet.FindOnePoint("from", Point(0,0,0)); // from --> Light location --> Only this needs to be multiply updated.	
	int _bVoxelLight = paramSet.FindOneInt("VoxelLight", 0);
	float backgroundColor = paramSet.FindOneFloat("backgroundI", 0.001f);
    Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world; // This has to be computed for every position.
	if(vgd == NULL) { printf("Volume Pointer Null. Get out.\n"); exit(-12); }
	multiplePointLights ->reserve(MAX_LIT_VOXELS);
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
	for(vector<Point>::iterator it = lightPos.begin() ; it < lightPos.end() ; it++, colorIt++) {
		Point p = (*it);
		Color3d c = *colorIt;
		float fc[3] = {c.r, c.g, c.b};
		I = Spectrum::FromRGB(fc);
		Transform l2w = Translate(Vector(p.x, p.y, p.z)) * light2world;
		Transform w2l = Translate(Vector(-p.x, -p.y, -p.z)) * light2world;
		multiplePointLights->push_back(new PointLight(l2w, I, backgroundColor, _bVoxelLight)); //dynamically created lights have also the same power as indicated.
		
	}

	I = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.
    return new PointLight(l2w, I * sc, backgroundColor, _bVoxelLight);
}

PointLight *CreatePointLight(const Transform &light2world,
        const ParamSet &paramSet, vector<Point> &lightPos, vector<Light*> *multiplePointLights) {
    Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0)); // scale
    Point P = paramSet.FindOnePoint("from", Point(0,0,0)); // from --> Light location --> Only this needs to be multiply updated.
	float backgroundColor = paramSet.FindOneFloat("backgroundI", 0.001f);
    Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world; // This has to be computed for every position.
	multiplePointLights ->reserve(MAX_LIT_VOXELS);
	for(vector<Point>::iterator it = lightPos.begin() ; it < lightPos.end() ; it++) {
		Point p = (*it);
		Transform l2w = Translate(Vector(p.x, p.y, p.z)) * light2world;
		multiplePointLights ->push_back(new PointLight(l2w, I * sc, backgroundColor, false)); //dynamically created lights have also the same power as indicated.
	}

    return new PointLight(l2w, I * sc, backgroundColor, false);
}


float PointLight::Pdf(const Point &, const Vector &) const {
    return 0.;
}

Spectrum PointLight::Sample_L(const Scene *scene, const LightSample &ls,
        float u1, float u2, float time, Ray *ray, Normal *Ns,
        float *pdf) const {
    *ray = Ray(lightPos, UniformSampleSphere(ls.uPos[0], ls.uPos[1]),
               0.f, INFINITY, time);
    *Ns = (Normal)ray->d;
    *pdf = UniformSpherePdf();
    return Intensity;
}


void PointLight::SHProject(const Point &p, float pEpsilon, int lmax,
        const Scene *scene, bool computeLightVisibility, float time,
        RNG &rng, Spectrum *coeffs) const {
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    if (computeLightVisibility &&
        scene->IntersectP(Ray(p, Normalize(lightPos - p), pEpsilon,
                              Distance(lightPos, p), time)))
        return;
    // Project point light source to SH
    float *Ylm = ALLOCA(float, SHTerms(lmax));
    Vector wi = Normalize(lightPos - p);
    SHEvaluate(wi, lmax, Ylm);
    Spectrum Li = Intensity / DistanceSquared(lightPos, p);
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = Li * Ylm[i];
}


