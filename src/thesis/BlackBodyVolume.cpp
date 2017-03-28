#include "stdafx.h"
#include "BlackBodyVolume.h"
#include "paramset.h"
#include "../thesis/VoxelBox.h"
#include "../thesis/Preprocessor.h"

#include <sstream>


void SanityCheck(const Spectrum & spectrum) {
	if( spectrum.At(0) < -1e-5 || spectrum.At(1) < -1e-5 || spectrum.At(2) < -1e-5) {
		printf("OMG. SanityCheck was not passed. The luminance is negative.\n");
		exit(-15);
	}
	else if( spectrum.HasNaNs() ) {
		printf("OMG. SanityCheck was not passed. NaN value here.\n");
		exit(-17);
	}
}


void ParsePbrtFile(string & filename, float *volumedata, int _nx, int _ny, int _nz) {
	ifstream fileH;
	string token;

	fileH.open(filename);
	if(!fileH.is_open() || !fileH.good()) {
		printf("File does not exist. What???");
		exit(-31);
	}

	getline(fileH, token);
	getline(fileH, token);
	getline(fileH, token);

	Color3d zeroColor(0.,0.,0.);
	float fData = -1.f;
	for(int z = 0 ; z < _nz ; z++) {
		for(int y = 0 ; y < _ny ; y++) {
			for(int x = 0 ; x < _nx ; x++) {
				fileH >> token;
				std::stringstream ss;
				ss << token;
				ss >> fData;
				
				if(ss.fail()) {
					printf("Wrong token value : %s", token);
					exit(-10);
				}
				volumedata[z * _ny * _nx + y * _nx + x] = fData;
			}
		}
	}
	fileH >> token;
	if(token != "]") {
		printf("There is something wrong with the parsing.\n");
		exit(-11);
	}

}

int nxBlackbodyMipmap = -1;
int nyBlackbodyMipmap = -1;
int nzBlackbodyMipmap = -1;

BlackBodyVolume::BlackBodyVolume(const Spectrum &sa, const Spectrum &ss, float gg,
	const Spectrum &emit, const BBox &e, const Transform &v2w,
	int x, int y, int z, const float *d, string & filenameAbsorption, string & filenameExtinction, float coefficientScale, bool _bStoreDataInVoxel)	// d == data for the volume(bucky)
	: DensityRegion(sa, ss, gg, emit, v2w), nx(x), ny(y), nz(z), extent(e)
{
	density = new float[nx*ny*nz];
	memcpy(density, d, nx*ny*nz*sizeof(float));
	
	absorption = new float[z * y * x];
	extinction = new float[z * y * x];
	scattering = new float[z * y * x];
	voxels = new VoxelBoxes(x, y, z, (BBox *)&this->extent);

	ParsePbrtFile(filenameAbsorption, absorption, this->nx, this->ny, this->nz);
	ParsePbrtFile(filenameExtinction, this->extinction, this->nx, this->ny, this->nz);

	// scattering coefficients are computed from absorption, extinction because scattering = extinction - absorption.
	for(int i = 0 ; i < z ; i++) {
		for(int j = 0 ; j < y ; j++) {
			for(int k = 0 ; k < x ; k++) {
				
				scattering[i * y * x + j * x + k] = extinction[i * y * x + j * x + k] - absorption[i * y * x + j * x + k];
				if( scattering[i * y * x + j * x + k] < 0.0f ) {
					scattering[i * y * x + j * x + k] = 0.0f;
				}
				extinction[i * y * x + j * x + k] = extinction[i * y * x + j * x + k] * coefficientScale;
				absorption[i * y * x + j * x + k] = absorption[i * y * x + j * x + k] * coefficientScale;
			}
		}
	}

	bStoreDataInVoxel = _bStoreDataInVoxel;
	if(_bStoreDataInVoxel) {
		// initialize lPower.
		float rgb[3] = {-1.f, -1.f, -1.f};

		lPower = (Spectrum **)arena.Alloc(sizeof(Spectrum) * nz * ny * nz); // x * y * z * nLights.
		for(int i = 0 ; i < nz ; i++) {
			for(int j = 0 ; j < ny ; j++) {
				for(int k = 0 ; k < nx ; k++) {
					lPower[i * ny * nx + j * nx + k] = (Spectrum *) arena.Alloc(sizeof(Spectrum) * MAX_LIT_VOXELS);
					for(int w = 0 ; w < MAX_LIT_VOXELS ; w++) {
						lPower[i * ny * nx + j * nx + k][w] = Spectrum::FromRGB(rgb);
					}
				}
			}
		}
	}

	// Alloc HRRPUVs.
	HRRPUVs = new float [x * y * z];

	// Alloc Gradient Vectors.
	gradientVec = new Vector[x * y * z];
	gradient = new float[x * y * z];

	//sigmaa, sigmas, sigmat
	sigmaa = new float[nx*ny*nz];
	sigmat = new float[nx*ny*nz];
	sigmas = new float[nx*ny*nz];
	for(int i = 0 ; i < nz ; i++) {
		for(int j = 0 ; j < ny ; j++) {
			for(int k = 0 ; k < nx ; k++) {
				float den = density[i * ny * nx + j * nx + k];
				sigmaa[i * ny * nx + j * nx + k] = absorption[i * ny * nx + j * nx + k] * den;
				sigmas[i * ny * nx + j * nx + k] = scattering[i * ny * nx + j * nx + k] * den;
				sigmat[i * ny * nx + j * nx + k] = extinction[i * ny * nx + j * nx + k] * den;
			}
		}
	}

}

BlackBodyVolume::~BlackBodyVolume()
{
	delete[] density;
	delete voxels;
	delete [] HRRPUVs;
	delete [] gradientVec;
	delete [] gradient;

	delete [] sigmaa;
	delete [] sigmas;
	delete [] sigmat;
}

// BlackBodyVolume Method Definitions
float BlackBodyVolume::Density(const Point &Pobj) const {
	if (!extent.Inside(Pobj)) {
		return 0.f;
	}
	// Compute voxel coordinates and offsets for _Pobj_
	Vector vox = extent.Offset(Pobj);
	vox.x = vox.x * nx - 0.5;	// Originally -0.5f. But I think that it is PBRT's bug.(It's from the original code.)
	vox.y = vox.y * ny - 0.5; // Originally -0.5f.
	vox.z = vox.z * nz - 0.5; // Originally -0.5f.
	//Assert(Range(vox.x, 0.0f, (float)nx) && Range(vox.y, 0.0f, (float)ny) && Range(vox.z, 0.0f, (float)nz));
	if(vox.x < 0.f) vox.x = 0.f;
	if(vox.y < 0.f) vox.y = 0.f;
	if(vox.z < 0.f) vox.z = 0.f;

#if TAEKYUDEBUG
	if(!Range(vox.x, 0.0f, (float)nx-EPSILON_F) || !Range(vox.y, 0.0f, (float)ny-EPSILON_F) || !Range(vox.z, 0.0f, (float)nz-EPSILON_F) ) {
		char *msg = new char [100];
		printf("vox.x = %f, vox.y = %f, vox.z = %f.\nnx = %d, ny = %d, nz = %d\n", vox.x, vox.y, vox.z, nx, ny, nz);
		exit(-31);
	}
#endif
	//if(vox  is outside extent)
	//		tweak vx, vy, vz values. Like manually push vox.xyz inside the volume and set vx/vy/vz to abs(vox.x - vx). It's okay to get the boundary wrong. No big deal.

	int vx = Floor2Int(vox.x), vy = Floor2Int(vox.y), vz = Floor2Int(vox.z);
	float dx = vox.x - vx, dy = vox.y - vy, dz = vox.z - vz;

	// Trilinearly interpolate density values to compute local density
	float d00 = Lerp(dx, D(vx, vy, vz),     D(vx+1, vy, vz));
	float d10 = Lerp(dx, D(vx, vy+1, vz),   D(vx+1, vy+1, vz));
	float d01 = Lerp(dx, D(vx, vy, vz+1),   D(vx+1, vy, vz+1));
	float d11 = Lerp(dx, D(vx, vy+1, vz+1), D(vx+1, vy+1, vz+1));
	float d0 = Lerp(dy, d00, d10);
	float d1 = Lerp(dy, d01, d11);
	return Lerp(dz, d0, d1);
}

inline VoxelPoint BlackBodyVolume::GetVoxelIndex(const Point &p) const
{
	// Compute voxel coordinates and offsets for _p_
	Vector vox = extent.Offset(WorldToVolume(p));
	vox.x = vox.x * nx;
	vox.y = vox.y * ny;
	vox.z = vox.z * nz;
	int vx = Floor2Int(vox.x), vy = Floor2Int(vox.y), vz = Floor2Int(vox.z);
	vx = (vx == nx) ? vx - 1 : vx; // vx = [0, 32]. vx could be 32(==nx). so make it 31 in that case.
	vy = (vy == ny) ? vy - 1 : vy;
	vz = (vz == nz) ? vz - 1 : vz;
	return VoxelPoint(vx, vy, vz);
}

BBox BlackBodyVolume::GetVoxelInVolume(int x, int y, int z) {		// By Taekyu Shin
	if(nx <= x || x < 0 || ny <= y || y < 0 || nz <= z || z < 0) {
		printf("Ohhh... Outside the volume is looked up GetVoxelInVolume of BlackbodyVolume\n");
		exit(-31);
	}
	Point coordinates(x, y, z);
	Vector maxVec = (extent.pMax - extent.pMin);
	Point oneVoxelSize(maxVec.x / nx, maxVec.y / ny, maxVec.z / nz);
	Point minP(oneVoxelSize.x * coordinates.x, oneVoxelSize.y * coordinates.y, oneVoxelSize.z * coordinates.z);
	Point maxP = minP + oneVoxelSize;
	return BBox(minP, maxP);
}

Spectrum BlackBodyVolume::sigma_a(const Point &p, const Vector &, float) const {
	int x, y, z;
	GetVoxelIndex(p, x, y, z);
	//Spectrum result = Density(WorldToVolume(p)) * absorption[z * ny * nx + y * nx + x];
#if TAEKYUDEBUG
	SanityCheck(sigmaa[z * ny * nx + y * nx + x]);
#endif
	return sigmaa[z * ny * nx + y * nx + x];
}

Spectrum BlackBodyVolume::sigma_a(int _x, int _y, int _z, const Vector & vec) const
{
#if TAEKYUDEBUG
	SanityCheck(sigmaa[_z * ny * nx + _y * nx + _x]);
#endif
	return sigmaa[_z * ny * nx + _y * nx + _x];
}

Spectrum BlackBodyVolume::sigma_s(const Point &p, const Vector &, float) const {
	int x, y, z;
	GetVoxelIndex(p, x, y, z);
	//Spectrum result = Density(WorldToVolume(p)) * scattering[z * ny * nx + y * nx + x];
#if TAEKYUDEBUG
	SanityCheck(sigmas[z * ny * nx + y * nx + x]);
#endif
	return sigmas[z * ny * nx + y * nx + x];
}

Spectrum BlackBodyVolume::sigma_s(int _x, int _y, int _z, const Vector & vec) const {
#if TAEKYUDEBUG
	SanityCheck(sigmas[_z * ny * nx + _y * nx + _x]);
#endif
	return sigmas[_z * ny * nx + _y * nx + _x];
}

Spectrum BlackBodyVolume::sigma_t(const Point &p, const Vector &, float) const {
	int x, y, z;
	this->GetVoxelIndex(p, x, y, z);
	//Spectrum result = Density(WorldToVolume(p)) * extinction[z * ny * nx + y * nx + x];
#if TAEKYUDEBUG
	SanityCheck(sigmat[z * ny * nx + y * nx + x]);
#endif
	return sigmat[z * ny * nx + y * nx + x];
}

Spectrum BlackBodyVolume::sigma_t(int _x, int _y, int _z, const Vector & vec) const {
#if TAEKYUDEBUG
	SanityCheck(sigmat[_z * ny * nx + _y * nx + _x]);
#endif
	return sigmat[_z * ny * nx + _y * nx + _x];
}

Spectrum BlackBodyVolume::GetDensity(const Point &p) const
{
	int x, y, z;
	this->GetVoxelIndex(p, x, y, z);
#if TAEKYUDEBUG
	SanityCheck(density[z * ny * nx + y * nx + x]);
#endif
	return density[z * ny * nx + y * nx + x];
}

Spectrum BlackBodyVolume::GetDensity(int _x, int _y, int _z) const
{
#if TAEKYUDEBUG
	SanityCheck(density[_z * ny * nx + _y * nx + _x]);
#endif
	return density[_z * ny * nx + _y * nx + _x];
}

bool BlackBodyVolume::GetPowerStoredInVoxel(int _nLight, int _x, int _y, int _z, Spectrum & outPower)
{
	Spectrum power = lPower[_z * ny * nx + _y * nx + _x][_nLight];
	float rgb[3];
	power.ToRGB(rgb);
	if(rgb[0] == -1.f && rgb[1] == -1.f && rgb[2] == -1.f) {
		return false;
	}
	else {
		outPower = power;
		return true;
	}
}

void BlackBodyVolume::StorePowerInVoxel(int lightIndex, int _x, int _y, int _z, const Spectrum & inPower)
{
	lPower[_z * ny * nx + _y * nx + _x][lightIndex] = inPower;
}

BlackBodyVolume *gBlackBodyVol = NULL;

BlackBodyVolume *CreateBlackbodyGridVolumeRegion(const Transform &volume2world,
	const ParamSet &params, string & absorpFilename, string & extinctFilename) {
		// Initialize common volume region parameters
		Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
		Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
		float g = params.FindOneFloat("g", 0.); // When g = 0 , it means ISO tropic scattering.
		float coefficientScale = params.FindOneFloat("extinctionScale", 1.f);
		
		Spectrum Le = params.FindOneSpectrum("Le", 0.);

		if(sigma_a.At(0) == 0.0 || sigma_s.At(0) == 0.0 || g == 0.0)
		{
			printf("You better not use 0.0 for these.\n");
			exit(-13);
		}

		Point p0 = params.FindOnePoint("p0", Point(0,0,0));
		Point p1 = params.FindOnePoint("p1", Point(1,1,1));
		int nitems;
		const float *data = params.FindFloat("density", &nitems); //nItem == # of density values , data == the bucky volume data pointer
		if (!data) {
			Error("No \"density\" values provided for volume grid?");
			return NULL;
		}
		int nx = params.FindOneInt("nx", 1);
		int ny = params.FindOneInt("ny", 1);
		int nz = params.FindOneInt("nz", 1);
		nxBlackbodyMipmap = nx;
		nyBlackbodyMipmap = ny;
		nzBlackbodyMipmap = nz;
		if (nitems != nx*ny*nz) {
			Error("BlackBodyVolume has %d density values but nx*ny*nz = %d",
				nitems, nx*ny*nz);
			return NULL;
		}

		bool bDataInVoxels = params.FindOneBool("StoreDataInVoxel", false);
		gBlackBodyVol = new BlackBodyVolume(sigma_a, sigma_s, g, Le, BBox(p0, p1),
			volume2world, nx, ny, nz, data, absorpFilename, extinctFilename, coefficientScale, bDataInVoxels);
		return gBlackBodyVol;
}
