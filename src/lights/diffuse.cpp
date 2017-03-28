// lights/diffuse.cpp*
#include "stdafx.h"
#include "lights/diffuse.h"
#include "paramset.h"
#include "montecarlo.h"

// DiffuseAreaLight Method Definitions
DiffuseAreaLight::~DiffuseAreaLight() {
    delete shapeSet;
}


DiffuseAreaLight::DiffuseAreaLight(const Transform &light2world,
        const Spectrum &le, int ns, const Reference<Shape> &s)
    : AreaLight(light2world, ns) {
    Lemit = le;
    shapeSet = new ShapeSet(s);
    area = shapeSet->Area();

	lightPos = Point(-10000.f, -10000.f, -10000.f);
}


Spectrum DiffuseAreaLight::Power(const Scene *) const {
    return Lemit * area * M_PI;
}


AreaLight *CreateDiffuseAreaLight(const Transform &light2world, const ParamSet &paramSet,
        const Reference<Shape> &shape) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
	else if(PbrtOptions.doubleQuickRender) nSamples = max(1, nSamples / 4);
    return new DiffuseAreaLight(light2world, L * sc, nSamples, shape);
}


Spectrum DiffuseAreaLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_AREA_LIGHT_STARTED_SAMPLE();
    Normal ns;
    Point ps = shapeSet->Sample(p, ls, &ns);
    *wi = Normalize(ps - p);
    *pdf = shapeSet->Pdf(p, *wi);
    visibility->SetSegment(p, pEpsilon, ps, 1e-3f, time);
    Spectrum Ls = L(ps, ns, -*wi);
    PBRT_AREA_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float DiffuseAreaLight::Pdf(const Point &p, const Vector &wi) const {
    return shapeSet->Pdf(p, wi);
}


Spectrum DiffuseAreaLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_AREA_LIGHT_STARTED_SAMPLE();
    Point org = shapeSet->Sample(ls, Ns);
    Vector dir = UniformSampleSphere(u1, u2);
    if (Dot(dir, *Ns) < 0.) dir *= -1.f;
    *ray = Ray(org, dir, 1e-3f, INFINITY, time);
    *pdf = shapeSet->Pdf(org) * INV_TWOPI;
    Spectrum Ls = L(org, *Ns, dir);
    PBRT_AREA_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


