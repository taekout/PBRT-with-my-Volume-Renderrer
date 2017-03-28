#pragma once

#include "../core/stdafx.h"
#include "../renderers/metropolis.h"
#include "HenyeyGreensteinSampleGenerator.h"

#define THRESHOLD_CDF	0.4f
#define NO_SOBOLPOINTS 1024

class MediaBidir
{
private:
	MetropolisRenderer *pMetroRenderer;
	int nLargeRays;
	int nMaxRayMarchings;//nSmallRays;
	float raymarchingEyeDistanceScale;
	float raymarchingLightDistanceScale;
	bool bVoxelization;

public:
	// public data
	HenyeyGreensteinSampleGenerator *pHGGeneralSamples;

	// public method
	MediaBidir(MetropolisRenderer *pRenderer, int nLargeRays, int nSmallRays, 
						float _raymarchingEyeDistanceScale, float _raymarchingLightDistanceScale, bool _bVoxelization);
	virtual ~MediaBidir(void);

	// Bidir func for Media mutation and ray marching.
	static void FreeVector(vector<vector<MediaPath *>> &paths);
	static void FreeVector(vector<SobolPoints3D *> &paths);
	Spectrum BidirPathGeneration(int nClosestLightNo, const Scene *scene, vector<SobolPoints3D *> &scatterLargeDirs, vector<SobolPoints3D *> lightLargeDirs,
										vector<vector<MediaPath *> > & scatterRayPaths, vector< vector < vector<MediaPath *> > > &lightPaths, float stepSize,
										RNG &rng, Point scatterPoint, MemoryArena &arena, const Spectrum &lightIntensity);
	Spectrum BidirOnePathGeneration(const Scene * scene, const int nLight, SobolPoints3D * &scatterLargeDirection,
										SobolPoints3D *lightLargeDirection, vector<MediaPath *> & scatterRayPath,
										vector<MediaPath *> &lightPath, float stepSize, RNG &rng, Point scatterPoint,
										MemoryArena &arena, const Spectrum &lightIntensity);
	Spectrum Transmittance(const Scene *scene, const Point &p0, const Point &p1, RNG &rng, float stepSize);
	Spectrum tau(const Scene *scene, const Point &p0, const Point &p1, const RNG &rng, float stepSize, float u) const;
	Spectrum tauForUnconnected(const Scene *scene, const Point &p0, const Point &p1, const RNG &rng, float stepSize, float u) const;
	Spectrum LBidir(const Scene * scene, int nClosestLight, const HenyeyGreensteinSampleGenerator * hg, const MemoryArena &arena,
					vector<vector<MediaPath *> > &lightScatterPath, vector<vector<MediaPath *> > &eyeScatterPath,
					const Spectrum &lightIntensity, RNG & rng, float stepSize);
	Spectrum LPath(const Scene * scene, int nClosestLight, vector<MediaPath *> &oneLightPath, vector<MediaPath *> &oneEyePath,
						const Spectrum &lightIntensity, RNG &rng, float stepSize);
	Spectrum LPathOnlyGettingCloser(const Scene * scene, int nClosestLight, vector<MediaPath *> &oneLightPath, vector<MediaPath *> &oneEyePath,
						const Spectrum &lightIntensity, RNG &rng, float stepSize);
	void ConnectTheUnconnected(const Point &src, const Point &dst, const RNG & rng);
};
/*
inline static float Distance(Point a, Point b) {
	float x = a.x - b.x;
	float y = a.y - b.y;
	float z = a.z - b.z;
	return sqrt(x*x + y*y + z*z);
}*/
