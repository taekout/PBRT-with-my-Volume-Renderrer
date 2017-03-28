#include "VoxelBox.h"

#define ERRORF -1.f

VoxelBox::VoxelBox(void)
{
	float rgb[3] = {ERRORF};
	Tau.FromRGB(rgb);
	Intensity.FromRGB(rgb);
}

VoxelBox::~VoxelBox(void)
{}

void VoxelBox::SetTau(const Spectrum &_Tau)
{
	Tau = _Tau;
}

void VoxelBox::SetIntensity(const Spectrum &_Intensity)
{
	Intensity = _Intensity;
}

void VoxelBox::AddIntensity(const Spectrum & _intensity)
{
	Spectrum prevIntensity = Intensity;
	Intensity += _intensity;
	Assert(prevIntensity.At(0) <= Intensity.At(0) &&
		   prevIntensity.At(1) <= Intensity.At(1) &&
		   prevIntensity.At(2) <= Intensity.At(2));
}

void VoxelBox::AddL2(const Spectrum & _l2)
{
	Spectrum prevL2 = l_2;
	l_2 += _l2;
	Assert(prevL2.At(0) <= l_2.At(0) &&
		prevL2.At(1) <= l_2.At(1) &&
		prevL2.At(2) <= l_2.At(2));
}

VoxelBoxes::VoxelBoxes(int _nx, int _ny, int _nz, BBox * _pBox)
{
	this->nx = _nx;
	this->ny = _ny;
	this->nz = _nz;
	grids = arena.Alloc<VoxelBox **> (_nz);//new VoxelBox ** [_nz];
	for(int i = 0 ; i < _nz ; i++) {
		grids[i] = arena.Alloc<VoxelBox *> (_ny);//new VoxelBox *[_ny];
		for(int j = 0 ; j < _ny ; j++) {
			grids[i][j] = arena.Alloc<VoxelBox > (_nx);//new VoxelBox[_nz];
		}
	}

	mpBox = _pBox;
	Point pMin = mpBox->pMin;
	Point pMax = mpBox->pMax;
	Vector Width = mpBox->pMax - mpBox->pMin;
	Width.x /= _nx;
	Width.y /= _ny;
	Width.z /= _nz; // Now we have width, height of one voxel.
	Vector HalfWidth = Width / 2.f;
	// Compute Voxel Centers for each voxel.
	for(int i = 0 ; i < _nz ; i++) {
		for(int j = 0 ; j < _ny ; j++) {
			for(int k = 0 ; k < _nx ; k++) {
				Point _center;
				_center.x = pMin.x + k * Width.x + HalfWidth.x;
				_center.y = pMin.y + j * Width.y + HalfWidth.y;
				_center.z = pMin.x + i * Width.z + HalfWidth.z;
				this->grids[i][j][k].center = _center;
			}
		}
	}
}

VoxelBoxes::~VoxelBoxes(void)
{
	arena.FreeAll();
}

