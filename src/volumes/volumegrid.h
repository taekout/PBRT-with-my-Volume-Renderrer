#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_VOLUMES_VOLUMEGRID_H
#define PBRT_VOLUMES_VOLUMEGRID_H

// volumes/volumegrid.h*
#include "volume.h"
#include "../thesis/VoxelBox.h"

class Preprocessor;

// VolumeGridDensity Declarations
class VolumeGridDensity : public DensityRegion {
public:
    // VolumeGridDensity Public Methods
    VolumeGridDensity(const Spectrum &sa, const Spectrum &ss, float gg,
            const Spectrum &emit, const BBox &e, const Transform &v2w,
            int x, int y, int z, const float *d, float _stepSize); // d == data for the volume(bucky)

    ~VolumeGridDensity() { delete[] density; delete voxels; }
    BBox WorldBound() const { return Inverse(WorldToVolume)(extent); }
    bool IntersectP(const Ray &r, float *t0, float *t1) const {
        Ray ray = WorldToVolume(r);
        return extent.IntersectP(ray, t0, t1); // IntersectP is BBox::IntersectP() The book p.193- 195
    }
    float Density(const Point &Pobj) const;
    float D(int x, int y, int z) const {
        x = Clamp(x, 0, nx-1);
        y = Clamp(y, 0, ny-1);
        z = Clamp(z, 0, nz-1);
        return density[z*nx*ny + y*nx + x];
    }

	Point GetPointInVolume(int x, int y, int z) { // By Taekyu Shin
		Point p(0.,0.,0.); Point coordinates(x, y, z);
		p += (extent.pMax - extent.pMin);
		p = p / Point(nx,ny,nz);
		p.x *= x; p.y *= y; p.z *= z;
		p += extent.pMin;
		return p;
	}

	BBox GetVoxelInVolume(int x, int y, int z);	// By Taekyu Shin

	inline VoxelPoint GetVoxelIndex(const Point &p) const {
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

	inline void GetVoxelIndex(const Point &p, int & outX, int & outY, int & outZ) {
		// Compute voxel coordinates and offsets for _p_
		Vector vox = extent.Offset(WorldToVolume(p));
		vox.x = vox.x * nx;
		vox.y = vox.y * ny;
		vox.z = vox.z * nz;
		int vx = Floor2Int(vox.x), vy = Floor2Int(vox.y), vz = Floor2Int(vox.z);
		outX = (vx == nx) ? vx - 1 : vx; // vx = [0, 32]. vx could be 32(==nx). so make it 31 in that case.
		outY = (vy == ny) ? vy - 1 : vy;
		outZ = (vz == nz) ? vz - 1 : vz;
	}

	inline void GetVoxelIndex(const Point &p, int & outNX, int & outNY, int & outNZ, float & outDX, float & outDY, float & outDZ) {
		// Compute voxel coordinates and offsets for _p_
		Vector vox = extent.Offset(WorldToVolume(p));
		outDX = vox.x * nx;
		outDY = vox.y * ny;
		outDZ = vox.z * nz;
		int vx = Floor2Int(outDX), vy = Floor2Int(outDY), vz = Floor2Int(outDZ);
		// vx = [0, 32]. vx could be 32(==nx). so make it 31 in that case.
		outNX = vx;
		outNY = vy;
		outNZ = vz;
	}

	inline void GetVoxelIndex(const Point &p, float & outX, float & outY, float & outZ) {
		Vector vox = extent.Offset(WorldToVolume(p));
		outX = vox.x * nx;
		outY = vox.y * ny;
		outZ = vox.z * nz;
	}

	inline VoxelBoxes * GetVoxels(void) {
		return voxels;
	}

	inline void GetDimension(int &_nx, int &_ny, int &_nz) { _nx = nx; _ny = ny; _nz = nz; }

	inline bool IsInsideExtent(const Point &p) {
		return extent.Inside(/*WorldToVolume(*/p/*)*/); // For Thesis purposes, it does not need WorldToVolume because I assume all points are with w = 1.0.
	}

	inline void GetExtentSize(float & x, float & y, float & z) {
		x = extent.pMax.x - extent.pMin.x;
		y = extent.pMax.y - extent.pMin.y;
		z = extent.pMax.z - extent.pMin.z;
	}

private:
    // VolumeGridDensity Private Data
    float *density;
    const int nx, ny, nz;
    const BBox extent;
	VoxelBoxes *voxels;

public:
	// Taekyu Shin
	Scene * mScene;
};



VolumeGridDensity *CreateGridVolumeRegion(const Transform &volume2world,
        const ParamSet &params);

#endif // PBRT_VOLUMES_VOLUMEGRID_H
