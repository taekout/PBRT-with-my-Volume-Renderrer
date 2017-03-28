#include "stdafx.h"
#include "blackbodylight.h"
#include "sh.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"

// BlackbodyLight Method Definitions
BlackbodyLight::BlackbodyLight(const Transform &light2world,
	const Spectrum &intensity)
	: Light(light2world), nLights(0) {
}

Spectrum BlackbodyLight::Sample_L(const Point &p, float pEpsilon,
	const LightSample &ls, float time, Vector *wi, float *pdf,
	VisibilityTester *visibility) const {
	return Spectrum(0.f);
}


Spectrum BlackbodyLight::Power(const Scene *) const {
	return 0.f;
}


BlackbodyLight *CreateFireLights(const Transform &light2world, const ParamSet &paramSet,
	vector<Point> &lightPos, vector<Color3d> &lightColor, vector<Light*> *multipleBlackbodyLights) {
		Spectrum I;// = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.
		Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0)); // scale
		Point P = paramSet.FindOnePoint("from", Point(0,0,0)); // from --> Light location --> Only this needs to be multiply updated.
		Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world; // This has to be computed for every position.
		multipleBlackbodyLights ->reserve(MAX_LIT_VOXELS);
		//add fire light ==> I.
		vector<Color3d>::iterator colorIt = lightColor.begin();
		for(vector<Point>::iterator it = lightPos.begin() ; it < lightPos.end() ; it++, colorIt++) {
			Point p = (*it);
			Color3d c = *colorIt;
			float fc[3] = {c.r, c.g, c.b};
			I = Spectrum::FromRGB(fc);
			Transform l2w = Translate(Vector(p.x, p.y, p.z)) * light2world;
			multipleBlackbodyLights ->push_back(new BlackbodyLight(l2w, I)); //dynamically created lights have also the same power as indicated.
		}

		return new BlackbodyLight(l2w, I * sc);
}

BlackbodyLight *CreateBlackbodyLight(const Transform &light2world,
	const ParamSet &paramSet, vector<Point> &lightPos, vector<Light*> *multipleBlackbodyLights) {
		Spectrum I = paramSet.FindOneSpectrum("I", Spectrum(1.0)); // Light Color and Radiance itself.
		Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0)); // scale
		Point P = paramSet.FindOnePoint("from", Point(0,0,0)); // from --> Light location --> Only this needs to be multiply updated.
		Transform l2w = Translate(Vector(P.x, P.y, P.z)) * light2world; // This has to be computed for every position.
		multipleBlackbodyLights ->reserve(MAX_LIT_VOXELS);
		for(vector<Point>::iterator it = lightPos.begin() ; it < lightPos.end() ; it++) {
			Point p = (*it);
			Transform l2w = Translate(Vector(p.x, p.y, p.z)) * light2world;
			multipleBlackbodyLights ->push_back(new BlackbodyLight(l2w, I * sc)); //dynamically created lights have also the same power as indicated.
		}

		return new BlackbodyLight(l2w, I * sc);
}


float BlackbodyLight::Pdf(const Point &, const Vector &) const {
	return 0.;
}


Spectrum BlackbodyLight::Sample_L(const Scene *scene, const LightSample &ls,
	float u1, float u2, float time, Ray *ray, Normal *Ns,
	float *pdf) const {
	return Spectrum(0.f);
}


void BlackbodyLight::SHProject(const Point &p, float pEpsilon, int lmax,
	const Scene *scene, bool computeLightVisibility, float time,
	RNG &rng, Spectrum *coeffs) const {

}

void BlackbodyLight::AddLight(Point p, Spectrum intensity, Spectrum s_s, Spectrum s_a, Spectrum den)
{
	nLights ++;
	this->lightIntensity.push_back(new Spectrum(intensity));
	this->lightPositions.push_back(new Point(p));
	this->sigma_a.push_back(new Spectrum(s_a));
	this->sigma_s.push_back(new Spectrum(s_s));
	this->sootDensity.push_back(new Spectrum(den));
}
