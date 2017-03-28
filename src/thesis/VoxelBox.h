#pragma once
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "memory.h"

struct VoxelPoint {
public:
	int nx, ny, nz;

	VoxelPoint(void) : nx(-1), ny(-1), nz(-1) {}
	VoxelPoint(int _x, int _y, int _z) : nx(_x), ny(_y), nz(_z) {}
	virtual ~VoxelPoint(void) {}

	inline void SetIndices(int _nx, int _ny, int _nz) {
		nx = _nx; ny = _ny; nz = _nz;
	}

	bool operator ==(const VoxelPoint &vp) {
		return (nx == vp.nx) && (ny == vp.ny) && (nz == vp.nz);
	}

	VoxelPoint operator = (const VoxelPoint &p) {
		nx = p.nx; ny = p.ny; nz = p.nz;
		return VoxelPoint(nx, ny, nz);
	}
};

struct VoxelBox {
	// For now, assume that there is only one light.
	Point center;
	Spectrum Tau;
	Spectrum Intensity;
	Spectrum l_2;

	VoxelBox(void);
	virtual ~VoxelBox(void);

	void SetTau(const Spectrum &_Tau);
	void SetIntensity(const Spectrum &_Intensity);
	inline void SetL2(const Spectrum &_l2) {
		l_2 = _l2;
	}
	void AddIntensity(const Spectrum & _intensity);
	void AddL2(const Spectrum & _l2);

	// Get Functions
	inline Point GetCenter() { return center; }
};

class VoxelBoxes {
protected:
	VoxelBox ***grids;
	int nx, ny, nz;
	MemoryArena arena;
	static const int INVALID_V = -10000;

	BBox * mpBox;

public:
	VoxelBoxes(int _nx, int _ny, int _nz, BBox *_pBox);
	virtual ~VoxelBoxes(void);

	// Get Functions
	inline Point GetCenter(int x, int y, int z) {
		return grids[z][y][x].center;
	}
	inline void GetDimension(int *x, int *y, int *z) {
		*x = nx ; *y = ny; *z = nz; 
	}
	inline VoxelBox * GetVoxel(int x, int y, int z) {
		if(x >= nx || y >= ny || z >= nz || x < 0 || y < 0 || z < 0) {
			printf("Voxel Index outside the voxel extent is accessed.\n%d, %d, %d\n", x, y, z);
			exit(-8);
		}
		return &grids[z][y][x];
	}
	inline VoxelBox * GetVoxel(const VoxelPoint & vp) {
		if(vp.nx < nx || vp.ny < ny || vp.nz < nz) {
			printf("Voxel Index outside the voxel extent is accessed.\n%d, %d, %d\n", vp.nx, vp.ny, vp.nz);
			exit(-9);
		}
		return &grids[vp.nz][vp.ny][vp.nx];
	}
};
