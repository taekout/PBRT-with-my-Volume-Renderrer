#include "stdafx.h"
#include "SphericalLight.h"
#include "paramset.h"
#include "montecarlo.h"
#include "HenyeyGreensteinSampleGenerator.h"
#include <time.h>

// SphericalLight Method Definitions

SphericalLight::SphericalLight(Point _center, float _radius, Spectrum _le, HenyeyGreensteinSampleGenerator *_pHGMultipleScatteringRayGenerator) {
		Lemit = _le;
		lightPos = _center;
		this->radius = _radius;
		sphere = NULL;
		rng = new RNG(clock());
		sphere->Area();

		hgSampleGenerator = _pHGMultipleScatteringRayGenerator;
		nSampleHG = hgSampleGenerator->GetSampleNumber();
}

SphericalLight::~SphericalLight() {
	if(rng)
		delete rng;
}

Spectrum SphericalLight::Power(void) {
	return Lemit * area * M_PI;
}

float SphericalLight::Pdf(const Point &p, const Vector &wi) const {
	return sphere->Pdf(p, wi);
}

Spectrum SphericalLight::Sample_L(const Point &p, float pEpsilon,
	const LightSample &ls, float time, Vector *wi, float *pdf,
	VisibilityTester *visibility) const {
		Normal ns;
		Vector sampleVec = hgSampleGenerator->GetPoint3D(rng->RandomUInt() % nSampleHG);
		Point ps = lightPos+ sampleVec * radius * rng->RandomFloat();
		*wi = Normalize(ps - p);
		*pdf = sphere->Pdf(p, *wi);
		visibility->SetSegment(p, pEpsilon, ps, 1e-3f, time);
		Spectrum Ls = L(ps, ns, -*wi);
		
		return Ls;
}


Spectrum SphericalLight::Sample_L(const Scene *scene,
	const LightSample &ls, float u1, float u2, float time,
	Ray *ray, Normal *Ns, float *pdf) const {

		//Vector dir = UniformSampleSphere(u1, u2);
		Vector dir = hgSampleGenerator->GetPoint3D(rng->RandomUInt() % nSampleHG );
		dir = dir * rng->RandomFloat() * radius;
		if (Dot(dir, *Ns) < 0.) dir *= -1.f;
		*ray = Ray(lightPos, dir, 1e-3f, INFINITY, time);
		*pdf = sphere->Pdf(lightPos, dir) * INV_TWOPI;
		Spectrum Ls = L(lightPos, *Ns, dir);
		
		return Ls;
}


