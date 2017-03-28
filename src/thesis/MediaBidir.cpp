#include "MediaBidir.h"
#include "../lights/point.h"
#include "../volumes/volumegrid.h"
#include "../thesis/HenyeyGreensteinSampleGenerator.h"
#include "../core/transform.h"

MediaBidir::MediaBidir(MetropolisRenderer *pRenderer, int nLargeRays, int nMaxRayMarchings, 
					float _raymarchingEyeDistanceScale, float _raymarchingLightDistanceScale, bool _bVoxelization)
{
	pMetroRenderer = pRenderer;
	this ->nLargeRays = nLargeRays;
	this ->nMaxRayMarchings = nMaxRayMarchings;
	this ->raymarchingEyeDistanceScale = _raymarchingEyeDistanceScale;
	this ->raymarchingLightDistanceScale = _raymarchingLightDistanceScale;
	this ->bVoxelization = _bVoxelization;

	// Slice Sampler Initialization.
	this ->pHGGeneralSamples = new HenyeyGreensteinSampleGenerator(NO_SOBOLPOINTS, 0.01f);
	pHGGeneralSamples ->SortByDotProd();
	vector<SobolPoints3D *> *pVec = pHGGeneralSamples ->GetAll3DVectors();
	unsigned int Size = pVec ->size();
	float CDF[NO_SOBOLPOINTS];
	if(Size < 1) exit(-1);
	CDF[0] = pVec ->at(0) ->pdf;
	for(unsigned int i = 1 ; i < Size ; i++) {
		SobolPoints3D *pPoint = pVec ->at(i);
		CDF[i] = CDF[i - 1] + pPoint ->pdf;
	}
}


MediaBidir::~MediaBidir(void)
{
}

/* Density func
delta = + b * e ^ ( + log(density * b / a) ); // b = 1/128, a = 1/1024.
originally it was : x +- b * e ^ ( - log(b/a) * gamma );
*/
static void mutateForDensity(RNG &rng, float *v, float density, float min = -1.f,
                          float max = 1.f)
{
    const float a = 1.f / 1024.f, b = 1.f / 128;
	const float logTerm =  density * b / a;
    const float logRatio = logf(1.f + logTerm);
    float delta = (max - min) * b * /*expf*/(logRatio);
    if (rng.RandomFloat() < 0.5f) {
        *v += delta;
//        if (*v >= max) *v = min + (*v - max);
    }
    else {
        *v -= delta;
//        if (*v < min) *v = max - (min - *v);
    }
	if (*v < min || *v >= max) *v = min;
}

void SmallStepWhileScattering(RNG &rng, Vector & v, float density, bool bidirectional) {
        float length = v.Length();
        v.Normalize();
//      mutateForDensity(rng, &r.maxt, density, 0.f, 1.f);
//      mutateForDensity(rng, &r.mint, density, 0.f, 1.f);
        mutateForDensity(rng, &v.x, density, -1.f, 1.f);
        mutateForDensity(rng, &v.y, density, -1.f, 1.f);
        mutateForDensity(rng, &v.z, density, -1.f, 1.f);
        v *= length;
        // We don't mess up with the origin. It stays the same.
}

void MediaBidir::FreeVector(vector<vector<MediaPath *>> &paths)
{
	int nPaths = paths.size();
	for(int i = 0 ; i < nPaths ; i++) {
		vector<MediaPath *> &pathVec = paths.at(i);
		int nVec = pathVec.size();
		for( int j = 0 ; j < nVec ; j++) {
			MediaPath * mp  = pathVec.at(j);
			if(mp != NULL) delete mp;
		}
	}
}

//vector<SobolPoints3D *>
void MediaBidir::FreeVector(vector<SobolPoints3D *> &paths)
{
	int nPaths = paths.size();
	for(int i = 0 ; i < nPaths ; i++) {
		SobolPoints3D * sp = paths.at(i);
		if(sp != NULL) delete sp;
	}
}

Spectrum MediaBidir::Transmittance(const Scene *scene, const Point &p0, const Point &p1, RNG &rng, float stepSize)
{
    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
	step = 4.f * stepSize;
	offset = rng.RandomFloat();
	Spectrum ta = tau(scene, p0, p1, rng, step, offset);
    return Exp(-ta);
}

Spectrum MediaBidir::tau(const Scene *scene, const Point &p0, const Point &p1, const RNG &rng, float stepSize, float u) const { // u == randomFloat
	float t0, t1; //t0 = curPoint, t1= next point.
	//rn == normalized ray. from o with direction r.d
	Vector dir(p1-p0);
	float length = dir.Length();
	if (length == 0.f) return 0.f;
	Ray rn(p0, dir, 0.f); // check rn has the right length. check the beginPoint and endPoint.
	//Ray rn(p0, r.d / length, r.mint * length, r.maxt * length, r.time);
	//if (!IntersectP(rn, &t0, &t1)) return 0.;
	t0 = 0.f; t1 = dir.Length();
	Spectrum tau(0.);
	t0 += u * stepSize;//u * stepSize;
	while (t0 < t1) {
		tau += scene ->volumeRegion ->sigma_t(rn(t0), -rn.d, rn.time);
		//tau += sigma_t(rn(t0), -rn.d, r.time);
		t0 += stepSize;
	}
	return tau * stepSize;
}

Spectrum MediaBidir::tauForUnconnected(const Scene *scene, const Point &p0, const Point &p1, const RNG &rng, float stepSize, float u) const
{
	float t0, t1; //t0 = curPoint, t1= next point.
	//rn == normalized ray. from o with direction r.d
	Vector dir(p1-p0);
	float length = dir.Length();
	if (length == 0.f) return 0.f;
	Ray rn(p0, dir, 0.f); // check rn has the right length. check the beginPoint and endPoint.
	//Ray rn(p0, r.d / length, r.mint * length, r.maxt * length, r.time);
	//if (!IntersectP(rn, &t0, &t1)) return 0.;
	t0 = 0.f; t1 = dir.Length();
	Spectrum tau(0.);
	t0 += u * stepSize;//u * stepSize;
	while (t0 < t1) {
		tau += scene ->volumeRegion ->sigma_t(rn(t0), -rn.d, rn.time); // Think of multiplying it with 10.0f to make sure unconnected paths are expensive.
		//tau += sigma_t(rn(t0), -rn.d, r.time);
		t0 += stepSize;
	}
	return tau * stepSize;
}