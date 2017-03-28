#include "BlackBodyMultiScattering.h"
#include "scene.h"
#include "paramset.h"
#include "parallel.h"
#include ".\lights\point.h"
#include "Preprocessor.h"
#include ".\thesis\BlackBodyVolume.h"
#include <exception>
#include <sstream>

#define WEIGHT 1

BlackBodyMultiScattering::BlackBodyMultiScattering(float _stepSize
	, int _nxMipmap, int _nyMipmap, int _nzMipmap) : thetaSquare(0.24f), stepSize(_stepSize)
{
	int minN = min(min(_nxMipmap, _nyMipmap), _nzMipmap);
	int _level = Ceil2Int(Log2((float)minN)) + 1; 
	mipmap = new TMipmap(_nxMipmap, _nyMipmap, _nzMipmap, _level);
	mutex = Mutex::Create();

	float tmp[3] = {1.0f , 1.0f , 1.0f};
	unitSpectrum = Spectrum::FromRGB(tmp);
}


BlackBodyMultiScattering::~BlackBodyMultiScattering(void)
{
}

Spectrum BlackBodyMultiScattering::Li(const Scene *scene, const Renderer *renderer,
	const RayDifferential &_ray, const Sample *sample, RNG &rng,
	Spectrum *transmittance, MemoryArena &arena, MediaBidir *pMediaBidir) const
{
	try{
	float t0, t1;
	BlackBodyVolume * vgd = (BlackBodyVolume *)scene->volumeRegion;
	if (!vgd || !vgd->IntersectP(_ray, &t0, &t1) || (t1-t0) == 0.f) {
		*transmittance = 1.f;
		return 0.f;
	}

	// Do single scattering volume integration in _vr_
	Ray ray = _ray;
	ray.d.Normalize();

	// Prepare for volume integration stepping
	int nSamples = Ceil2Int((t1-t0) / stepSize);
	float step = (t1 - t0) / nSamples;
	//float step = 0.f;
	Spectrum Tr(1.f);
	Point p = ray(t0), pPrev;
	Point startingP = p;
	Vector w = -ray.d;
	// XXX: This should be revisited. The reason is that step size is determined later. So the addition here might turn out to be too huge.
	t0 += sample->oneD[scatterSampleOffset][0] * step;
	
	// Compute sample patterns for single scattering samples
	float *lightNum = arena.Alloc<float>(nSamples);
	LDShuffleScrambled1D(1, nSamples, lightNum, rng);
	float *lightComp = arena.Alloc<float>(nSamples);
	LDShuffleScrambled1D(1, nSamples, lightComp, rng);
	float *lightPos = arena.Alloc<float>(2*nSamples);
	LDShuffleScrambled2D(1, nSamples, lightPos, rng);
	int nMaxSamples = nSamples;
	uint32_t sampOffset = 0;
	Vector emptyVec;
	Spectrum l1(0.f);
	VoxelBoxes * vb = vgd->GetVoxels();
	Spectrum LTotal(0.f);

	/*
	1. What should be the normalization factor for the distance?
	2. BlurWidth is too big. Higher level mipmapping.
	*/
	///////////// Later, You must change the step Size!!!!!!
	for (int i = 0 ; t0 < t1/*i < nSamples*/ ; ++i, t0 += step) {

		// Advance to sample at _t0_ and update _T_
		if(t0 + 1e-3f > t1)
			continue;
		if(step < 0.05)
			printf("Step size less than 0.05:%f\n",step);
		pPrev = p;
		p = ray(t0);
		//Thesis Code.
		int nLights = scene->lights.size();
		int ln = min(Floor2Int(lightNum[sampOffset] * nLights),
			nLights-1);
		PointLight *light = dynamic_cast<PointLight *>(scene->lights[ln]);
		if(light == NULL) {
			throw "Light type is not a point light.";
		}
		float s1 = Distance(startingP, p); // (1)
		float s2 = Distance(p, light->GetPoint()); // (2)
		float S = s1 + s2; // (3)
		// Sometimes, precision problems can happen, so double-check it here.
		if(!vgd->IsInsideExtent(p))
			continue;
		Spectrum b = vgd->sigma_s(p, emptyVec, 1.f); // (4)
		Spectrum a = vgd->sigma_a(p, emptyVec, 1.f);
		/*Spectrum dl = b * Distance(pPrev, p); //(4)-Approximated! Theoretically speaking, it should be added while ray-marching, but instead,it just multiplies
		l1 = l1 + dl; // (5)*/
		l1 = b * Distance(startingP, p); // (4), (5)
		// (6) has been already computed in Precomputation stage.
		int x, y, z;
		float fX, fY, fZ;
		vgd->GetVoxelIndex(p, x, y, z, fX, fY, fZ);
		Spectrum l2 = vb->GetVoxel(x, y, z)->l_2;
		Spectrum l = l1 + l2; // (7)
		// ComputeBlurWidth - returns too big. Maybe normalization of distance.
		Spectrum blurWidth = ComputeBlurWidth(thetaSquare, l.At(0), S, a.At(0), b.At(0), l2.At(0)); // (8) - Note: I think a is based on the single point, not the sum thru the path.
		if(l.At(0) != l.At(2)) {
			printf("Ipossibl! 0 and 2 are different ?");

		}
		int nLevel = Log2(blurWidth.At(0));
		nLevel = (nLevel > mipmap->Levels()) ? mipmap->Levels() : nLevel;
		nLevel = (nLevel < 1) ? 1 : nLevel;
		//nLevel = (mipmap->Levels()) - nLevel + 1;
		Spectrum L = this->mipmap->LookUp(nLevel, fX, fY, fZ); // (9) p is converted to voxel indices(floating point)
		//Spectrum L = this->mipmap->LAttenuatedBasedOnDistance(p, light);
#if MULTDEBUG
		RayDifferential debugRay(_ray);
		bool bHit = Debug_Intersect(light->GetPoint(), debugRay);
		if(bHit == true) {
			//printf("This ray goes towards the lights\n");
		}
#endif
		float theta = acos(Dot(ray.d, Normalize(light->GetPoint() - p))); // (10)
		if(IsCloseWithEpisilon(ray.d.Length(), 1.0f) == false) {
			std::string errMsg = "direction vector not normalized. Length : ";
			std::ostringstream ss;
			ss << ray.d.Length();
			errMsg += ss.str();
			throw errMsg.c_str();
		}
		//////////////////////////////////
		// Now figure out P_MS --> Notice that when the fog density is very high, it is very isotropic scattering.(g=0.0 for HG.)
		Spectrum P = P_MS(theta, l); // (11)
		// alpha = .5f * thetaSquare
		
		float C = 1.0 / l.At(0) * (.5f * thetaSquare) * theta * theta; // (12)  --> .5 * thetaSquare = alpha.
		float c = a.At(0) + b.At(0);
		float weight = exp(-(c / b.At(0)) * l.At(0)) * exp(-C) * P.At(0); // (13)
#if !WEIGHT
		weight = 1.0f;
#endif
		LTotal = LTotal + (L * weight); // Make sure here if they add it correctly.
#if MULTDEBUG
		++sampOffset;
		sampOffset = (sampOffset == nMaxSamples) ? 0 : sampOffset;
#endif
		float sintheta = sin(theta);
		sintheta = (sintheta < 1e-4) ? 1e-4 : sintheta;
		step = blurWidth.At(0) / sintheta;
	}
	*transmittance = Tr;

	return LTotal;
	}
	catch(char *errorStr) {
		printf("Exception while multi-scattering : %s\n", errorStr);
		exit(-11);
	}
}

Spectrum BlackBodyMultiScattering::Transmittance(const Scene *scene,
	const Renderer *renderer, const RayDifferential &ray,
	const Sample *sample, RNG &rng, MemoryArena &arena) const
{
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

Spectrum BlackBodyMultiScattering::LiForNormalSingleScatterRay(const Scene *scene, const Renderer *renderer,
	const RayDifferential &ray, const Sample *sample, RNG &rng,
	Spectrum *transmittance, MemoryArena &arena, MediaBidir *pMediaBidir) const
{
	return 0.f;
}

void BlackBodyMultiScattering::RequestSamples(Sampler *sampler, Sample *sample,
	const Scene *scene) {
		tauSampleOffset = sample->Add1D(1);
		scatterSampleOffset = sample->Add1D(1);
}

Spectrum BlackBodyMultiScattering::ComputeBlurWidth(float thetaSquare, const Spectrum l, float S, const Spectrum a, const Spectrum b, const Spectrum l2) const
{
	if(b.At(0) == 0.f || b.At(1) == 0.f || b.At(2) == 0.f || 
		l2.At(0) == 0.f || l2.At(1) == 0.f || l2.At(2) == 0.f) {
		std::string errorstr;
		errorstr = "divide-by-zero exception in the ComputeBlurWidth() function...";
		std::stringstream ss, ss2;
		ss << b.At(0);
		errorstr += ss.str() + "...";
		ss2 << l2.At(0);
		errorstr += ss2.str();
		throw errorstr.c_str();
	}
	float tmp[3];
	for(int i = 0 ; i < 3 ; i++) {
		tmp[i] = thetaSquare * l.At(i) * S * S / (24.f * (l.At(i) * thetaSquare * a.At(i) / b.At(i) * l.At(i) * l.At(i) / l2.At(i) ) ); //divide-by-zero error might occur.
#if MULTDEBUG
		if(tmp[i] < 0.f)
			throw "sqrt(negative number) occurs.";
#endif
		tmp[i] = sqrt(tmp[i]);
	}
	Spectrum ret(tmp[0], tmp[1], tmp[2]);
	return ret;
}

Spectrum BlackBodyMultiScattering::P_MS(float theta, const Spectrum & l) const
{
	double N = 0.000029; // gotta experiment it and set the normalization constant N manually. // When g = 0.01, max_P = 0.082. When g = 0.001, max_P = 0.798, When g = 0.64, max_P = 1.0070.
	float scatteringEvents = l.At(0);
	if(scatteringEvents < 0.0f) 
		throw "scatteringEvents Negative in P_MS() functions.\n";
	char errmsg[1200];
	if(_isnan(scatteringEvents)) 
		throw "scattering Events NaN\n";
	if(_isnan(sqrt(1.0f - exp(-scatteringEvents))))
		throw "P_MS sqrt causes NaN\n";
	if(_isnan(sqrt(scatteringEvents / sqrt(1.0f - exp(-scatteringEvents))))) {
		throw "Is NaN Spectrum occurs in P_MS\n";
	}
	sprintf(errmsg, "%f", scatteringEvents);

	Spectrum result = P( 1.0f / sqrt(scatteringEvents / (1.0f - exp(-scatteringEvents)) ) * theta ) * ( N );
	if(result.HasNaNs())
		throw "Result = NaN in P_MS()\n";
	if(result.At(0) < 0.0)
		throw "Result < 0.0 in P_MS()\n";

	{
		MutexLock mlock(*mutex);
		this->mPreprocessor->AddPhaseNormalization(result.At(0));
	}

	return result;
	// In the original paper : Impact of Multiple scattering on Stimulated Infrared Cloud Scene Images
	// 1 / N * P(theta * sqrt(scatteringEvents / (1 - exp(-scatteringEvents)) );
}

float BlackBodyMultiScattering::P(float theta) const
{
	// Henyey so that it is easier to compute N.
	const float g = 0.01f; // When g = 0.01, max_P = 0.082. When g = 0.001, max_P = 0.798, When g = 0.64, max_P = 1.0070.
	// 1. / (4 * PI) = 0.0795774715459424
	float result = 0.0795774715459424f * (1 - g*g) / pow(1.f + g*g - 2.f * g * cos(theta), 1.5f);
	if(_isnan(result)) {
		printf("%f\n", g);
		printf("%f\n", theta);
		printf("%f\n", cos(theta));
		printf("%f\n", 1.f + g*g - 2.f * g * cos(theta));
		throw "NaN result in P() \n";
	}
	if(result < 0.0f)
		throw "Negative result in P() \n";
	return 0.0795774715459424f * (1 - g*g) / pow(1.f + g*g - 2.f * g * cos(theta), 1.5f);
}

bool BlackBodyMultiScattering::Debug_Intersect(Point &p, RayDifferential & _ray) const
{
	Vector voxExtent(2.5f, 2.5f, 2.5f);
	BBox *pBox = new BBox(p - voxExtent, p + voxExtent);
	float t0, t1;
	bool bHit = pBox->IntersectP(_ray, &t0, &t1);
	return bHit;
}

extern int nxBlackbodyMipmap;
extern int nyBlackbodyMipmap;
extern int nzBlackbodyMipmap;

extern BlackBodyVolume * gBlackBodyVol;

BlackBodyMultiScattering *CreateBlackBodyMultiScattering(const ParamSet &params)
{
	float stepSize  = params.FindOneFloat("stepsize", 1.f);
	if(!gBlackBodyVol) {
		printf("Exception : gBlackBodyVol has not been initialized.");
		exit(-13);
	}
	return new BlackBodyMultiScattering(stepSize
		, nxBlackbodyMipmap, nyBlackbodyMipmap, nzBlackbodyMipmap);
}

