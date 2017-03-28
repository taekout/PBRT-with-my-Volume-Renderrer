#pragma once
#include "pbrt.h"
#include "VoxelBox.h"
#include "../lights/point.h"
#include "TMipmap.h"

class VolumeGridDensity;

class Preprocessor
{
public:
	int nPhaseEvents;
	double sumN ;
	void AddPhaseNormalization(double n) {
		nPhaseEvents ++;
		sumN += n;
	}
	Preprocessor(VolumeGridDensity *_vgd, float _stepSize, bool _bRedoPrepressessing);
	~Preprocessor(void);

	void PreprocessVolume(TMipmap *mipmap, vector<Light *> & lights);
//	Spectrum MultipleScatteringLi(/*const Scene *scene, const Renderer *renderer,*/
//		const RayDifferential &ray, /*const Sample *sample, */RNG &rng, Spectrum *T) const;
protected:
	void ComputeLightAttenuation(VolumeGridDensity *volGrids, VoxelBoxes *vb, vector<Light *> & lights);
	void BuildVolumePyramid(TMipmap *mipmap, VoxelBoxes *vb/*LightVolume */);

private:
	VolumeGridDensity * vgd;
	MemoryArena arena;

	

public:
	float stepSize;

	bool bNotPrintOutPreprocessData;
};

