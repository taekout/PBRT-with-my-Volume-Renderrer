#pragma once
#include "VoxelBox.h"

struct SobolPoints2D
{
	double x;
	double y;
	SobolPoints2D(void) {
		x = y = 0.0;
	}
};

struct SobolPoints3D
{
	double x;
	double y;
	double z;

	float pdf;
	float dotProd;
	SobolPoints3D(void) {
		x = y = z = 0.0;
		pdf = 0.f; dotProd = -2.f;
	}
	SobolPoints3D(double _x, double _y, double _z, float _pdf) {
		x = _x; y = _y; z = _z; pdf = _pdf;
		// Dot Prod with the pivot Vector(-1.f, 0.f, 0.f).
		dotProd = -1.f * float(x);
	}
	operator Vector() {
		return Vector(x,y,z);
	}
};

struct MediaPath {
	VoxelPoint vp;
	int nStep;
	Vector wPrev, wNext;
	Point pPrev, pCur;
	//BSDF *bsdf;
	Spectrum Tau;
	Spectrum Phase;
	Spectrum LightIntensity;

	MediaPath() {
		wPrev.x = wPrev.y = wPrev.z = INFINITY;
		wNext.x = wNext.y = wNext.z = INFINITY;
		pPrev.x = pPrev.y = pPrev.z = INFINITY;
		pCur.x = pCur.y = pCur.z = INFINITY;
		float rgb[3] = {0};
		Tau.FromRGB(rgb);
		Phase.FromRGB(rgb);
		LightIntensity.FromRGB(rgb);
		nStep = -1;
	}
	MediaPath(float _pPrev[3], float _pCur[3], float _wPrev[3], float _wNext[3]) {
		wPrev.x = _wPrev[0]; wPrev.y = _wPrev[1]; wPrev.z = _wPrev[2];
		wNext.x = _wNext[0]; wNext.y = _wNext[1]; wNext.z = _wNext[2];
		pPrev.x = _pPrev[0]; pPrev.y = _pPrev[1]; pPrev.z = _pPrev[2];
		pCur.x = _pCur[0]; pCur.y = _pCur[1]; pCur.z = _pCur[2];
		float rgb[3] = {0};
		Tau.FromRGB(rgb);
		Phase.FromRGB(rgb);
		LightIntensity.FromRGB(rgb);
		nStep = -1;
	}
	operator Vector() {
		return Vector(wNext.x, wNext.y, wNext.z);
	}
};

class HenyeyGreensteinSampleGenerator
{
private:
	void GenerateSobolPoints2D(unsigned int N); // First generate the Sobol 2D points.
	// Second use this function to generate 3D points from the previously generated 2D sobol points.
	void GenerateWellDistributedPointsOnSphere(float g);

public:
	HenyeyGreensteinSampleGenerator(unsigned int N, float g = 0.64f);
	virtual ~HenyeyGreensteinSampleGenerator(void);
	static void TransformVector(const Vector & curDir, vector<vector<SobolPoints3D *>> &paths,
				vector<vector<SobolPoints3D *>> &transformedPaths);
	static void TransformVector(const Vector & curDir, vector<SobolPoints3D *> &paths,
				vector<SobolPoints3D *> &transformedPaths, MemoryArena &arena);
	virtual void ToString(Vector &z, vector<vector<SobolPoints3D *>> &paths, vector<vector<SobolPoints3D* >> &paths2);
	virtual void ToString(Vector &z, vector<SobolPoints3D *> &paths, vector<SobolPoints3D *> &paths2);

	inline int GetSampleNumber(void) { return nPoints; }

	/*
** z_basis is the scatter vector direction. x and y are generated with cross product.
(x_basis, y_basis, z_basis) * (0,0,1)' = (IdentityMatrix) * (SobelSequence)';
*/
	void TransformDirectionIntoCurrentRaySpace(const Vector &scatterDir, // scatter ray direction.
		vector<Vector *> & resultRayDirectionInRaySpace)	// result Rays after space conversion.
	{
		Vector arbitrary(0.f, 1.f, 0.f);
		Vector x = Cross(scatterDir, arbitrary);
		Vector y = Cross(x, scatterDir);
		float dirx, diry, dirz;
		float hgvecx, hgvecy, hgvecz;
		for(int i = 0 ; i < nPoints ; i++)
		{
			hgvecx = points3d[i].x;
			hgvecy = points3d[i].y;
			hgvecz = points3d[i].z;
			dirx = x.x * hgvecx + y.x * hgvecy + scatterDir.x * hgvecz;
			diry = x.y * hgvecx + y.y * hgvecy + scatterDir.y * hgvecz;
			dirz = x.z * hgvecx + y.z * hgvecy + scatterDir.z * hgvecz;
			resultRayDirectionInRaySpace.push_back(new Vector(dirx, diry, dirz));
		}
	}
	inline Vector GetPoint3D(int index) {
		return Vector(points3d[index].x, points3d[index].y, points3d[index].z);
	}

	inline vector<SobolPoints3D *> * GetAll3DVectors(void) {
		return	point3dVec;
	}
	
	inline float GetPDF(int index) const {
		return points3d[index].pdf;
	}
	void SortByDotProd();

private: // Data
	int nPoints;
	SobolPoints2D *points2d;
	SobolPoints3D *points3d;
	// just add vector. vector<SobolPoint3D *> point3dvec; later get the vector point and feed it to LBidir function.
	vector<SobolPoints2D *> *point2dVec;
	vector<SobolPoints3D *> *point3dVec;
};