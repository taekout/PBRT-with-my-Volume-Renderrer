#include "TMipmap.h"
#include "..\lights\point.h"

TMipmap::TMipmap(int _nx, int _ny, int _nz, int _level)
{
	nx = _nx; ny = _ny; nz = _nz;
	level = _level;
	if(level < 1) { 
		printf("\n\n Mipmap Level must be higher than 1");
		exit(-5);
	}
	for(int l = 0 ; l < level ; l++) {
		RGBColor *** head = arena.Alloc<RGBColor **> (_nz);
		levelMaps.push_back(head);
		for(int i = 0 ; i < _nz ; i++) {
			head[i] = arena.Alloc<RGBColor *> (_ny);
			for(int j =  0 ; j < _ny ; j++) {
				head[i][j] = arena.Alloc<RGBColor > (_nx);
				for(int k = 0 ; k < _nx ; k++) {
					head[i][j][k].c[0] = head[i][j][k].c[1] = head[i][j][k].c[2] = initValue;
				}
			}
		}
		Halve(_nx, _ny, _nz);
	}
}


TMipmap::~TMipmap(void)
{
	arena.FreeAll();
}

RGBColor::operator Spectrum() const
{
	Spectrum ret(c[0], c[1], c[2]);
	return ret;
}

void TMipmap::InitializeMap(RGBColor ***volume)
{
	int _x = nx, _y = ny, _z = nz;

	if(level <= 0) {
		printf("Level must be higher than 0. So exit.\n");
		exit(-11);
	}
	for(int l = 0 ; l < level - 1 ; l++) {
		Blur2X2X2(levelMaps[l], levelMaps[l+1], _z, _y, _x);
		Halve(_x, _y, _z);
	}
}

Spectrum TMipmap::LookUp(int lookupLevel, float _nx, float _ny, float _nz)
{
	// lowest level = 0
	// Trilinear Interpolation.
	if(_nx >= nx || _ny >= ny || _nz >= nz ||
		_nx < 0.f || _ny < 0.f || _nz < 0.f || lookupLevel > level
		|| lookupLevel < 1) {
		printf("Looked up somewhere out of the voxel bound.\n");
		printf("%f, %f, %f, %d\n", _nx, _ny, _nz, lookupLevel);
		exit(-4);
	}

	float multiplicant = 1.f / powf(2, lookupLevel - 1);
	_nx *= multiplicant;
	_ny *= multiplicant;
	_nz *= multiplicant;
	int nUpperX = (int)ceilf(_nx);
	int nLowerX = (int)floorf(_nx);
	int nUpperY = (int)ceilf(_ny);
	int nLowerY = (int)floorf(_ny);
	int nUpperZ = (int)ceilf(_nz);
	int nLowerZ = (int)floorf(_nz);
	int maxIndexX, maxIndexY, maxIndexZ; 
	maxIndexX = Floor2Int((nx - 1) * multiplicant); // ???
	maxIndexY = Floor2Int((ny - 1) * multiplicant);
	maxIndexZ = Floor2Int((nz - 1) * multiplicant);
	float fUpperX = nUpperX - _nx;
	float fLowerX = _nx - nLowerX;
	float fUpperY = nUpperY - _ny;
	float fLowerY = _ny - nLowerY;
	float fUpperZ = nUpperZ - _nz;
	float fLowerZ = _nz - nLowerZ;
	if(nUpperX > maxIndexX) {
		nUpperX = maxIndexX;
	}
	if(nUpperY > maxIndexY) {
		nUpperY = maxIndexY;
	}
	if(nUpperZ > maxIndexZ) {
		nUpperZ = maxIndexZ;
	}
	if(nLowerX < 0) {
		nLowerX = 0;
	}
	if(nLowerY < 0) {
		nLowerY = 0;
	}
	if(nLowerZ < 0) {
		nLowerZ = 0;
	}

	RGBColor fInterpolated;
	RGBColor ***head = levelMaps.at(lookupLevel - 1);
	//if(levelMaps.size() == lookupLevel)
	/*if(level == lookupLevel) {
		return head[0][0][0]; // When the highest level, there is only only voxel. Return the one voxel.
	}*/
	bool bValueOkay = Range(fUpperX, 0.0f, 1.0f) && Range(fUpperY, 0.0f, 1.0f) && Range(fUpperZ, 0.0f, 1.0f)
		&& Range(fLowerX, 0.0f, 1.0f) && Range(fLowerY, 0.0f, 1.0f) && Range(fLowerZ, 0.0f, 1.0f);
	if(!bValueOkay)
		throw "Mipmapping look up values are wrong.\n";

	fInterpolated =
		head[nLowerZ][nLowerY][nLowerX] * fUpperZ * fUpperY * fUpperX +
		head[nLowerZ][nLowerY][nUpperX] * fUpperZ * fUpperY * fLowerX +
		head[nLowerZ][nUpperY][nLowerX] * fUpperZ * fLowerY * fUpperX +
		head[nLowerZ][nUpperY][nUpperX] * fUpperZ * fLowerY * fLowerX +
		head[nUpperZ][nLowerY][nLowerX] * fLowerZ * fUpperY * fUpperX +
		head[nUpperZ][nLowerY][nUpperX] * fLowerZ * fUpperY * fLowerX +
		head[nUpperZ][nUpperY][nLowerX] * fLowerZ * fLowerY * fUpperX +
		head[nUpperZ][nUpperY][nUpperX] * fLowerZ * fLowerY * fLowerX;

	return fInterpolated;
}

Spectrum TMipmap::LAttenuatedBasedOnDistance(Point p, PointLight * light)
{
	Point lightPos = light->GetPoint();
	float distance = (lightPos - p).Length();
	Spectrum inten = light->GetIntensity();
	inten = inten / distance / distance; //exp(-distance);
	return inten;
}

void TMipmap::Blur2X2X2(RGBColor ***srcMap, RGBColor ***dstMap, int _z, int _y, int _x )
{
	for(int i = 0 ; i < _z ; i+=2) {
		for(int j = 0 ; j < _y ; j+=2) {
			for(int k = 0 ; k < _x ; k+=2) {
				// interpolate between 8 voxels.
				float sum[3] = {0.0f};
				int nVoxelsToInterpolate = 0;
				for(int d0 = 0 ; d0 < 2 ; d0++) {
					for(int d1 = 0 ; d1 < 2 ; d1++) {
						for(int d2 = 0 ; d2 < 2 ; d2++) {
							if(i + d0 < _z && j + d1 < _y && k + d2 < _x) {
								nVoxelsToInterpolate  ++;
								sum[0] += srcMap[i + d0][j + d1][k + d2].c[0];
								sum[1] += srcMap[i + d0][j + d1][k + d2].c[1];
								sum[2] += srcMap[i + d0][j + d1][k + d2].c[2];
							}
						}
					}
				}
				dstMap[i / 2][j / 2][k / 2].c[0] = sum[0] / nVoxelsToInterpolate;
				dstMap[i / 2][j / 2][k / 2].c[1] = sum[1] / nVoxelsToInterpolate;
				dstMap[i / 2][j / 2][k / 2].c[2] = sum[2] / nVoxelsToInterpolate;
			}
		}
	}
}
