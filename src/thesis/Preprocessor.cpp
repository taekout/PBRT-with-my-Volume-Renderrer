#include "Preprocessor.h"
#include "../volumes/volumegrid.h"
#include "scene.h"
#include <fstream>
#include <time.h>

#define DISTANCE_BASED_EXTINCTION		1

Preprocessor::Preprocessor(VolumeGridDensity * _vgd, float _stepSize, bool _bRedoPreprocessing)
{
	this->nPhaseEvents = 0;
	vgd = _vgd;
	stepSize = _stepSize;
	int nx, ny, nz;
	vgd->GetDimension(nx, ny, nz);
	bNotPrintOutPreprocessData = _bRedoPreprocessing;
	sumN = 0.0;
}


Preprocessor::~Preprocessor(void)
{
}

void Preprocessor::PreprocessVolume(TMipmap *mipmap, vector<Light *> & lights)
{
	try {
	ComputeLightAttenuation(vgd, vgd->GetVoxels(), lights);
	BuildVolumePyramid(mipmap, vgd->GetVoxels());
	printf("Precomputation Done.\n");
	}
	catch(char * msg) {
		printf("%s", msg);
		exit(-1);
	}
}

/*
Spectrum Preprocessor::MultipleScatteringLi(const RayDifferential &ray, RNG &rng, Spectrum *T) const {
	float t0, t1;
	if (!this->vgd || !vgd->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
		*T = 1.f;
		return 0.f;
	}
	return 1.f;
}*/

/*
From each voxel, compute light transport from light.
*/
void Preprocessor::ComputeLightAttenuation(VolumeGridDensity *volGrids, VoxelBoxes *vb, vector<Light *> & lights)
{
	// For each voxel, compute Tau().
	int x, y, z;
	vb ->GetDimension(&x, &y, &z);
	std::ifstream inFile;
	inFile.open("preprocessdata.txt");
	if(inFile.is_open() == true && bNotPrintOutPreprocessData == false) {
		// read data from file.
		float rgb[3];
		int i , j , k;
		for(i = 0 ; i < z ; i++) {
			for(j = 0 ; j < y ; j++) {
				for(k = 0 ; k < x ; k++) {
					VoxelBox * oneGrid = vb->GetVoxel(k, j, i);
#if 0
					std::string line;
					std::getline(inFile, line);
					std::string oneFloat;
					int pos = 0, nextpos = 0;
					for(int g = 0 ; g < 3 ; g++) {
						nextpos = line.find(' ', pos);
						oneFloat = line.substr(pos, nextpos - pos);
						pos = nextpos + 1;
						rgb[g] = atof(oneFloat.data());
					}
					Spectrum _Intensity;
					_Intensity.FromRGB(rgb);
					oneGrid->SetIntensity(_Intensity);
					for(int g = 0 ; g < 3 ; g++) {
						nextpos = line.find(' ', pos); // nextpos is set to -1 at the end of the string, but it is still fine and it will work.
						oneFloat = line.substr(pos, nextpos - pos);
						pos = nextpos + 1;
						rgb[g] = atof(oneFloat.data());
					}
					Spectrum _L2;
					_L2.FromRGB(rgb);
					oneGrid->SetL2(_L2);
#endif // 0
					inFile >> rgb[0];
					inFile >> rgb[1];
					inFile >> rgb[2];
					if(rgb[0] == 0.f || rgb[1] == 0.f || rgb[2] == 0.f) {
						printf("Intensity can't be %f\n (%d, %d, %d)", rgb[0], k, j, i);
						exit(-3);
					}
					oneGrid->Intensity = Spectrum::FromRGB(rgb);
					inFile >> rgb[0];
					inFile >> rgb[1];
					inFile >> rgb[2];
					oneGrid->l_2 = Spectrum::FromRGB(rgb);
					if(rgb[0] == 0.f || rgb[1] == 0.f || rgb[2] == 0.f) {
						printf("L_2 can't be %f\n (%d, %d, %d)", rgb[0], k, j, i);
						exit(-3);
					}
				}
			}
		}
		if(i != z || j != y || k != z) {
			printf("Error. They must go through x * y * z\n");
			exit(-18);
		}
	}
	else {
		clock_t start = clock();
		std::ofstream outFile;
		if(bNotPrintOutPreprocessData == false) {
			outFile.open("preprocessdata.txt");
			if(outFile.is_open() == false) {
				printf("File not opened : preprocessdata.txt\n");
				exit(-30);
			}
		}
		Vector garbageVec(1.0f, 1.0f, 1.0f);
		int i , j , k;
		int NoOfLights = lights.size();
		VoxelBox * oneGrid = NULL;
		for(i = 0 ; i < z ; i++) {
			for(j = 0 ; j < y ; j++) {
				for(k = 0 ; k < x ; k++) {
					float fDistance;
					oneGrid = vb ->GetVoxel(k, j, i);
					for(int lightIndex = 0 ; lightIndex < NoOfLights ; lightIndex++) {

						const Point lightPos = ((PointLight *)lights.at(lightIndex)) ->GetPoint();

#if DISTANCE_BASED_EXTINCTION
						Spectrum sig_t = volGrids->sigma_t(lightPos, garbageVec, 1.0f);
						Spectrum sig_s = volGrids->sigma_s(lightPos, garbageVec, 1.0f);
#endif
						fDistance = Distance(oneGrid ->center, lightPos); // center not intiialized!?
						Ray r(lightPos, Normalize(oneGrid->center - lightPos), 0.f, INFINITY); // ray = from Light to Center of Voxel.
						Spectrum Tau_BasedOnRayMarching, l_2;
#if 0
						Spectrum Tau_BasedOnDensity = volGrids->tau(r, stepSize, 1.0f); //step Size should be 2.0 / grid_Dimension.
#endif
						//bool bSuccess = volGrids->tau_all_and_scattering(r, stepSize, 1.0f, l_2, Tau_BasedOnRayMarching);
						// Distance based Intensity works better.(Much better.) Why?
#if DISTANCE_BASED_EXTINCTION
						Tau_BasedOnRayMarching = sig_t * fDistance;
						l_2 = sig_s * fDistance;
#else
						if(volGrids->tau_all_and_scattering(r, stepSize, 1.0f, l_2, Tau_BasedOnRayMarching) == false) {
							throw "Somehow the length of the ray is 0.\n";
						}
#endif
						float fAttenuation[3] = {0.f};
						for(int temp = 0 ; temp < 3 ; temp++)
							fAttenuation[temp] = exp(-Tau_BasedOnRayMarching.At(temp));
						Spectrum multiplicant;
						multiplicant = Spectrum::FromRGB(fAttenuation);
						oneGrid->AddIntensity( ((PointLight *)lights.at(lightIndex))->GetIntensity() * Exp(-Tau_BasedOnRayMarching));
						oneGrid->AddL2(l_2);
					}

					if(bNotPrintOutPreprocessData == false) {
						float rgb[3];
//						outFile << lightIndex << ":";
						outFile << "Distance:" << fDistance;
						outFile << ":Tau(" << i << "," << j << "," << k << ")";
//						outFile << ":" << Tau_BasedOnRayMarching.At(0) << "," << Tau_BasedOnRayMarching.At(1) << "," << Tau_BasedOnRayMarching.At(2);
						//outFile << ":Sigma:" << sig_t.At(0) << "," << sig_t.At(1) << "," << sig_t.At(2);
						outFile << ":VoxelCenter:"<< oneGrid->GetCenter().x << "," << oneGrid->GetCenter().y << "," << oneGrid->GetCenter().z;
//						outFile << ":lightPos:" << lightPos.x << "," << lightPos.y << "," << lightPos.z << ": " << fDistance;

						oneGrid->Intensity.ToRGB(rgb);
						outFile << "Intensity:" << rgb[0] << ' ' << rgb[1] << ' ' << rgb[2];
						/*if(rgb[0] == 0.f || rgb[1] == 0.f || rgb[2] == 0.f) {
						printf("Intensity can't be %f\n (%d, %d, %d)", rgb[0], k, j, i);
						exit(-3);
						}*/
						oneGrid->l_2.ToRGB(rgb);
						outFile << ":L2:" << rgb[0] << " " << rgb[1] << " " << rgb[2] << '\n';
						/*if(rgb[0] == 0.f || rgb[1] == 0.f || rgb[2] == 0.f) {
						printf("L_2 can't be %f\n (%d, %d, %d)", rgb[0], k, j, i);
						exit(-3);
						}*/
					}
				}
			}
		}
		outFile.close();
		if(i != z || j != y || k != x) {
			printf("Error. They must go through x * y * z\n");
			exit(-18);
		}
		clock_t end = clock();
		printf("%f seconds for preprocessing.\n", ((float)(start - end)) / CLOCKS_PER_SEC);
	}
}

void Preprocessor::BuildVolumePyramid(TMipmap *mipmap, VoxelBoxes *vb)
{
	int x, y, z;
	int x2, y2, z2;
	// I assume that the index is divisible by 2.
	vb ->GetDimension(&x, &y, &z);
	RGBColor ***vol; float tmp[3];
	vol = mipmap->levelMaps[0];
	x2 = mipmap->X(); y2 = mipmap->Y(); z2 = mipmap->Z();
	if(x != x2 || y != y2 || z != z2) {
		printf("Volume size didn't match.\n");
		exit(-82);
	}

	for(int i = 0 ; i < z ; i++) {
		for(int j = 0 ; j < y ; j++) {
			for(int k = 0 ; k < x ; k++) {
				vb->GetVoxel(k, j, i)->Intensity.ToRGB(tmp);
				vol[i][j][k].c[0] = tmp[0];
				vol[i][j][k].c[1] = tmp[1];
				vol[i][j][k].c[2] = tmp[2];
			}
		}
	}
	mipmap->InitializeMap(vol);
}

