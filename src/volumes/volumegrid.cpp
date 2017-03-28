// volumes/volumegrid.cpp*
#include "stdafx.h"
#include "volumes/volumegrid.h"
#include "paramset.h"
#include "../thesis/VoxelBox.h"
#include "../thesis/Preprocessor.h"

int nxForMultipleScatteringMipmap = -1;
int nyForMultipleScatteringMipmap = -1;
int nzForMultipleScatteringMipmap = -1;

VolumeGridDensity::VolumeGridDensity(const Spectrum &sa, const Spectrum &ss, float gg,
            const Spectrum &emit, const BBox &e, const Transform &v2w,
            int x, int y, int z, const float *d, float _stepSize)	// d == data for the volume(bucky)
        : DensityRegion(sa, ss, gg, emit, v2w), nx(x), ny(y), nz(z), extent(e)
{
	density = new float[nx*ny*nz];
	memcpy(density, d, nx*ny*nz*sizeof(float));
	voxels = new VoxelBoxes(x, y, z, (BBox *)&this->extent);
}

// VolumeGridDensity Method Definitions
float VolumeGridDensity::Density(const Point &Pobj) const {
	if (!extent.Inside(Pobj)) {
		return 0.f;
	}
    // Compute voxel coordinates and offsets for _Pobj_
    Vector vox = extent.Offset(Pobj);
    vox.x = vox.x * nx - .5f;
    vox.y = vox.y * ny - .5f;
    vox.z = vox.z * nz - .5f;
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

BBox VolumeGridDensity::GetVoxelInVolume(int x, int y, int z) {	// By Taekyu Shin
	Point coordinates(x, y, z);
	Vector maxVec = (extent.pMax - extent.pMin);
	Point oneVoxelSize(maxVec.x / nx, maxVec.y / ny, maxVec.z / nz);
	Point minP(oneVoxelSize.x * coordinates.x, oneVoxelSize.y * coordinates.y, oneVoxelSize.z * coordinates.z);
	Point maxP = minP + oneVoxelSize;
	return BBox(minP, maxP);
}

VolumeGridDensity *gVolGridDensity = NULL;

VolumeGridDensity *CreateGridVolumeRegion(const Transform &volume2world,
        const ParamSet &params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.); // When g = 0 , it means ISO tropic scattering.
    Spectrum Le = params.FindOneSpectrum("Le", 0.);

	if(sigma_a.At(0) == 0.0 || sigma_s.At(0) == 0.0 || g == 0.0 )
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
	nxForMultipleScatteringMipmap = nx;
	nyForMultipleScatteringMipmap = ny;
	nzForMultipleScatteringMipmap = nz;
    if (nitems != nx*ny*nz) {
        Error("VolumeGridDensity has %d density values but nx*ny*nz = %d",
            nitems, nx*ny*nz);
        return NULL;
    }
	gVolGridDensity = new VolumeGridDensity(sigma_a, sigma_s, g, Le, BBox(p0, p1),
		volume2world, nx, ny, nz, data, -100000000.f);
    return gVolGridDensity;
}
