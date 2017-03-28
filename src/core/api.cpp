#include "stdafx.h"
#include "api.h"
#include "parallel.h"
#include "paramset.h"
#include "spectrum.h"
#include "scene.h"
#include "renderer.h"
#include "film.h"
#include "volume.h"
#include "probes.h"

// API Additional Headers
#include "accelerators/bvh.h"
#include "accelerators/grid.h"
#include "accelerators/kdtreeaccel.h"
#include "cameras/environment.h"
#include "cameras/orthographic.h"
#include "cameras/perspective.h"
#include "film/image.h"
#include "filters/box.h"
#include "filters/gaussian.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"
#include "filters/triangle.h"
#include "integrators/ambientocclusion.h"
#include "integrators/diffuseprt.h"
#include "integrators/dipolesubsurface.h"
#include "integrators/directlighting.h"
#include "integrators/emission.h"
#include "integrators/glossyprt.h"
#include "integrators/igi.h"
#include "integrators/irradiancecache.h"
#include "integrators/path.h"
#include "integrators/photonmap.h"
#include "integrators/single.h"
#include "thesis/MultipleScatteringIntegrator.h"
#include "thesis/BlackBodyMultiScattering.h"
#include "thesis/ForwardMultipleScattering.h"
#include "integrators/useprobes.h"
#include "integrators/whitted.h"
#include "lights/diffuse.h"
#include "lights/distant.h"
#include "lights/goniometric.h"
#include "lights/infinite.h"
#include "lights/point.h"
#include "thesis/VoxelLight.h"
#include "lights/projection.h"
#include "lights/spot.h"
#include "materials/glass.h"
#include "materials/kdsubsurface.h"
#include "materials/matte.h"
#include "materials/measured.h"
#include "materials/metal.h"
#include "materials/mirror.h"
#include "materials/mixmat.h"
#include "materials/plastic.h"
#include "materials/substrate.h"
#include "materials/subsurface.h"
#include "materials/shinymetal.h"
#include "materials/translucent.h"
#include "materials/uber.h"
#include "renderers/aggregatetest.h"
#include "renderers/createprobes.h"
#include "renderers/metropolis.h"
#include "renderers/samplerrenderer.h"
#include "renderers/surfacepoints.h"
#include "samplers/adaptive.h"
#include "samplers/bestcandidate.h"
#include "samplers/halton.h"
#include "samplers/lowdiscrepancy.h"
#include "samplers/random.h"
#include "samplers/stratified.h"
#include "shapes/cone.h"
#include "shapes/cylinder.h"
#include "shapes/disk.h"
#include "shapes/heightfield.h"
#include "shapes/hyperboloid.h"
#include "shapes/loopsubdiv.h"
#include "shapes/nurbs.h"
#include "shapes/paraboloid.h"
#include "shapes/sphere.h"
#include "shapes/trianglemesh.h"
#include "textures/bilerp.h"
#include "textures/checkerboard.h"
#include "textures/constant.h"
#include "textures/dots.h"
#include "textures/fbm.h"
#include "textures/imagemap.h"
#include "textures/marble.h"
#include "textures/mix.h"
#include "textures/scale.h"
#include "textures/uv.h"
#include "textures/windy.h"
#include "textures/wrinkled.h"
#include "volumes/exponential.h"
#include "volumes/homogeneous.h"
#include "volumes/volumegrid.h"
#include "thesis/BlackBodyVolume.h"
#include "../thesis/FireColor.h"
#include <map>

 #if (_MSC_VER >= 1400)
 #include <stdio.h>
 #define snprintf _snprintf
 #endif
using std::map;

// API Global Variables
Options PbrtOptions;

// API Local Classes
#define MAX_TRANSFORMS 2
#define START_TRANSFORM_BITS (1 << 0)
#define END_TRANSFORM_BITS   (1 << 1)
#define ALL_TRANSFORMS_BITS  ((1 << MAX_TRANSFORMS) - 1)
struct TransformSet {
   // TransformSet Public Methods
   Transform &operator[](int i) {
       Assert(i >= 0 && i < MAX_TRANSFORMS);
       return t[i];
   }
   const Transform &operator[](int i) const { Assert(i >= 0 && i < MAX_TRANSFORMS); return t[i]; }
   friend TransformSet Inverse(const TransformSet &ts) {
       TransformSet t2;
       for (int i = 0; i < MAX_TRANSFORMS; ++i)
           t2.t[i] = Inverse(ts.t[i]);
       return t2;
   }
   bool IsAnimated() const {
       for (int i = 0; i < MAX_TRANSFORMS-1; ++i)
           if (t[i] != t[i+1]) return true;
       return false;
   }
private:
    Transform t[MAX_TRANSFORMS];
};


struct RenderOptions {
    // RenderOptions Public Methods
    RenderOptions();
    Scene *MakeScene();
    Camera *MakeCamera() const;
    Renderer *MakeRenderer() const;

    // RenderOptions Public Data
    float transformStartTime, transformEndTime;
    string FilterName;
    ParamSet FilterParams;
    string FilmName;
    ParamSet FilmParams;
    string SamplerName;
    ParamSet SamplerParams;
    string AcceleratorName;
    ParamSet AcceleratorParams;
    string RendererName;
    string SurfIntegratorName, VolIntegratorName;
    ParamSet RendererParams;
    ParamSet SurfIntegratorParams, VolIntegratorParams;
    string CameraName;
    ParamSet CameraParams;
    TransformSet CameraToWorld;
    vector<Light *> lights;
	vector<bool>    IsLightsFireKind; // Fire Light Or Normal Light.
    vector<Reference<Primitive> > primitives;
    mutable vector<VolumeRegion *> volumeRegions;
    map<string, vector<Reference<Primitive> > > instances;
    vector<Reference<Primitive> > *currentInstance;
};


RenderOptions::RenderOptions() {
    // RenderOptions Constructor Implementation
    transformStartTime = 0.f;
    transformEndTime = 1.f;
    FilterName = "box";
    FilmName = "image";
    SamplerName = "lowdiscrepancy";
    AcceleratorName = "bvh";
    RendererName = "sampler";
    SurfIntegratorName = "directlighting";
    VolIntegratorName = "emission";
    CameraName = "perspective";
    currentInstance = NULL;
}


struct GraphicsState {
    // Graphics State Methods
    GraphicsState();
    Reference<Material> CreateMaterial(const ParamSet &params);

    // Graphics State
    map<string, Reference<Texture<float> > > floatTextures;
    map<string, Reference<Texture<Spectrum> > > spectrumTextures;
    ParamSet materialParams;
    string material;
    map<string, Reference<Material> > namedMaterials;
    string currentNamedMaterial;
    ParamSet areaLightParams;
    string areaLight;
    bool reverseOrientation;
};


GraphicsState::GraphicsState() {
    // GraphicsState Constructor Implementation
    material = "matte";
    reverseOrientation = false;
}


class TransformCache {
public:
    // TransformCache Public Methods
    void Lookup(const Transform &t, Transform **tCached,
            Transform **tCachedInverse) {
        map<Transform, std::pair<Transform *, Transform *> >::iterator iter;
        iter = cache.find(t);
        if (iter == cache.end()) {
            Transform *tr = arena.Alloc<Transform>();
            *tr = t;
            Transform *tinv = arena.Alloc<Transform>();
            *tinv = Transform(Inverse(t));
            cache[t] = std::make_pair(tr, tinv);
            iter = cache.find(t);
            PBRT_ALLOCATED_CACHED_TRANSFORM();
        }
        else
            PBRT_FOUND_CACHED_TRANSFORM();
        if (tCached) *tCached = iter->second.first;
        if (tCachedInverse) *tCachedInverse = iter->second.second;
    }
    void Clear() {
        arena.FreeAll();
        cache.erase(cache.begin(), cache.end());
    }
private:
    // TransformCache Private Data
    map<Transform, std::pair<Transform *, Transform *> > cache;
    MemoryArena arena;
};



// API Static Data
#define STATE_UNINITIALIZED  0
#define STATE_OPTIONS_BLOCK  1
#define STATE_WORLD_BLOCK    2
static int currentApiState = STATE_UNINITIALIZED;
static TransformSet curTransform;
static int activeTransformBits = ALL_TRANSFORMS_BITS;
static map<string, TransformSet> namedCoordinateSystems;
static RenderOptions *renderOptions = NULL;
static GraphicsState graphicsState;
static vector<GraphicsState> pushedGraphicsStates;
static vector<TransformSet> pushedTransforms;
static vector<uint32_t> pushedActiveTransformBits;
static TransformCache transformCache;

// API Macros
#define VERIFY_INITIALIZED(func) \
if (currentApiState == STATE_UNINITIALIZED) { \
    Error("pbrtInit() must be before calling \"%s()\". " \
          "Ignoring.", func); \
    return; \
} else /* swallow trailing semicolon */
#define VERIFY_OPTIONS(func) \
VERIFY_INITIALIZED(func); \
if (currentApiState == STATE_WORLD_BLOCK) { \
    Error("Options cannot be set inside world block; " \
          "\"%s\" not allowed.  Ignoring.", func); \
    return; \
} else /* swallow trailing semicolon */
#define VERIFY_WORLD(func) \
VERIFY_INITIALIZED(func); \
if (currentApiState == STATE_OPTIONS_BLOCK) { \
    Error("Scene description must be inside world block; " \
          "\"%s\" not allowed. Ignoring.", func); \
    return; \
} else /* swallow trailing semicolon */
#define FOR_ACTIVE_TRANSFORMS(expr) \
    for (int i = 0; i < MAX_TRANSFORMS; ++i) \
        if (activeTransformBits & (1 << i)) { expr }
#define WARN_IF_ANIMATED_TRANSFORM(func) \
do { if (curTransform.IsAnimated()) \
         Warning("Animated transformations set; ignoring for \"%s\"" \
                 "and using the start transform only", func); \
} while (false)

// Object Creation Function Definitions
Reference<Shape> MakeShape(const string &name,
        const Transform *object2world, const Transform *world2object,
        bool reverseOrientation, const ParamSet &paramSet) {
    Shape *s = NULL;

    if (name == "sphere")
        s = CreateSphereShape(object2world, world2object,
                              reverseOrientation, paramSet);
    // Create remaining _Shape_ types
    else if (name == "cylinder")
        s = CreateCylinderShape(object2world, world2object, reverseOrientation,
                                paramSet);
    else if (name == "disk")
        s = CreateDiskShape(object2world, world2object, reverseOrientation,
                            paramSet);
    else if (name == "cone")
        s = CreateConeShape(object2world, world2object, reverseOrientation,
                            paramSet);
    else if (name == "paraboloid")
        s = CreateParaboloidShape(object2world, world2object, reverseOrientation,
                                  paramSet);
    else if (name == "hyperboloid")
        s = CreateHyperboloidShape(object2world, world2object, reverseOrientation,
                                   paramSet);
    else if (name == "trianglemesh")
        s = CreateTriangleMeshShape(object2world, world2object, reverseOrientation,
                                    paramSet, &graphicsState.floatTextures);
    else if (name == "heightfield")
        s = CreateHeightfieldShape(object2world, world2object, reverseOrientation,
                                   paramSet);
    else if (name == "loopsubdiv")
        s = CreateLoopSubdivShape(object2world, world2object, reverseOrientation,
                                  paramSet);
    else if (name == "nurbs")
        s = CreateNURBSShape(object2world, world2object, reverseOrientation,
                             paramSet);
    else
        Warning("Shape \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return s;
}


Reference<Material> MakeMaterial(const string &name,
        const Transform &mtl2world,
        const TextureParams &mp) {
    Material *material = NULL;
    if (name == "matte")
        material = CreateMatteMaterial(mtl2world, mp);
    else if (name == "plastic")
        material = CreatePlasticMaterial(mtl2world, mp);
    else if (name == "translucent")
        material = CreateTranslucentMaterial(mtl2world, mp);
    else if (name == "glass")
        material = CreateGlassMaterial(mtl2world, mp);
    else if (name == "mirror")
        material = CreateMirrorMaterial(mtl2world, mp);
    else if (name == "mix") {
        string m1 = mp.FindString("namedmaterial1", "");
        string m2 = mp.FindString("namedmaterial2", "");
        Reference<Material> mat1 = graphicsState.namedMaterials[m1];
        Reference<Material> mat2 = graphicsState.namedMaterials[m2];
        if (!mat1) {
            Error("Named material \"%s\" undefined.  Using \"matte\"",
                  m1.c_str());
            mat1 = MakeMaterial("matte", curTransform[0], mp);
        }
        if (!mat2) {
            Error("Named material \"%s\" undefined.  Using \"matte\"",
                  m2.c_str());
            mat2 = MakeMaterial("matte", curTransform[0], mp);
        }

        material = CreateMixMaterial(mtl2world, mp, mat1, mat2);
    }
    else if (name == "metal")
        material = CreateMetalMaterial(mtl2world, mp);
    else if (name == "substrate")
        material = CreateSubstrateMaterial(mtl2world, mp);
    else if (name == "uber")
        material = CreateUberMaterial(mtl2world, mp);
    else if (name == "subsurface")
        material = CreateSubsurfaceMaterial(mtl2world, mp);
    else if (name == "kdsubsurface")
        material = CreateKdSubsurfaceMaterial(mtl2world, mp);
    else if (name == "measured")
        material = CreateMeasuredMaterial(mtl2world, mp);
    else if (name == "shinymetal")
        material = CreateShinyMetalMaterial(mtl2world, mp);
    else
        Warning("Material \"%s\" unknown.", name.c_str());
    mp.ReportUnused();
    if (!material) Error("Unable to create material \"%s\"", name.c_str());
    return material;
}


Reference<Texture<float> > MakeFloatTexture(const string &name,
        const Transform &tex2world, const TextureParams &tp) {
    Texture<float> *tex = NULL;
    if (name == "constant")
        tex = CreateConstantFloatTexture(tex2world, tp);
    else if (name == "scale")
        tex = CreateScaleFloatTexture(tex2world, tp);
    else if (name == "mix")
        tex = CreateMixFloatTexture(tex2world, tp);
    else if (name == "bilerp")
        tex = CreateBilerpFloatTexture(tex2world, tp);
    else if (name == "imagemap")
        tex = CreateImageFloatTexture(tex2world, tp);
    else if (name == "uv")
        tex = CreateUVFloatTexture(tex2world, tp);
    else if (name == "checkerboard")
        tex = CreateCheckerboardFloatTexture(tex2world, tp);
    else if (name == "dots")
        tex = CreateDotsFloatTexture(tex2world, tp);
    else if (name == "fbm")
        tex = CreateFBmFloatTexture(tex2world, tp);
    else if (name == "wrinkled")
        tex = CreateWrinkledFloatTexture(tex2world, tp);
    else if (name == "marble")
        tex = CreateMarbleFloatTexture(tex2world, tp);
    else if (name == "windy")
        tex = CreateWindyFloatTexture(tex2world, tp);
    else
        Warning("Float texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();
    return tex;
}


Reference<Texture<Spectrum> > MakeSpectrumTexture(const string &name,
        const Transform &tex2world, const TextureParams &tp) {
    Texture<Spectrum> *tex = NULL;
    if (name == "constant")
        tex = CreateConstantSpectrumTexture(tex2world, tp);
    else if (name == "scale")
        tex = CreateScaleSpectrumTexture(tex2world, tp);
    else if (name == "mix")
        tex = CreateMixSpectrumTexture(tex2world, tp);
    else if (name == "bilerp")
        tex = CreateBilerpSpectrumTexture(tex2world, tp);
    else if (name == "imagemap")
        tex = CreateImageSpectrumTexture(tex2world, tp);
    else if (name == "uv")
        tex = CreateUVSpectrumTexture(tex2world, tp);
    else if (name == "checkerboard")
        tex = CreateCheckerboardSpectrumTexture(tex2world, tp);
    else if (name == "dots")
        tex = CreateDotsSpectrumTexture(tex2world, tp);
    else if (name == "fbm")
        tex = CreateFBmSpectrumTexture(tex2world, tp);
    else if (name == "wrinkled")
        tex = CreateWrinkledSpectrumTexture(tex2world, tp);
    else if (name == "marble")
        tex = CreateMarbleSpectrumTexture(tex2world, tp);
    else if (name == "windy")
        tex = CreateWindySpectrumTexture(tex2world, tp);
    else
        Warning("Spectrum texture \"%s\" unknown.", name.c_str());
    tp.ReportUnused();
    return tex;
}


Light *MakeLightAutomaticPlacement(const string &name,
        const Transform &light2world, const ParamSet &paramSet, string &filenameTemperature, string &filenameHRRPUV, vector<Light*> *multipleLights, float threshold) {
    Light *light = NULL;
    if (name == "point" || name == "voxel")
	{
		// Taekyu Shin added these lines.
		float fLightScale = paramSet.FindOneFloat("lightScale", 1.0f);
		vector<Point> lightPos;
		vector<BBox> lightVoxels;
		vector<Color3d> lightColor;
		// Here, check volumeRegion in (RenderingOption).
		// Check the volume data and place virtual lights instead of the information from pbrt parsed data
		vector<VolumeRegion *>::iterator it;
		for(it = renderOptions ->volumeRegions.begin() ; it < renderOptions ->volumeRegions.end() ; it++) {
			FireColor fire;
			std::ifstream fileH, fileExtin, fileAbsorp;
			std::ifstream fileTemp;
			fileH.open(filenameHRRPUV);
			fileTemp.open(filenameTemperature);
			if(!fileH.is_open() || !fileH.good()) {
				printf("Smoke Fire HRRPUV file open error.\n");
				exit(-3);
			}
			else if(!fileTemp.is_open() || !fileTemp.good()) {
				printf("Smoke Temperature file open error.\n");
				exit(-3);
			}

			string token;
			int nx, ny, nz;
			fileH >> token; fileH >> token; fileH >> token; fileH >> token; fileH >> nx;
			fileH >> token; fileH >> token; fileH >> ny;
			fileH >> token; fileH >> token; fileH >> nz;
			getline(fileH, token); getline(fileH, token); getline(fileH, token);
			getline(fileTemp, token); getline(fileTemp, token); getline(fileTemp, token);
			getline(fileExtin, token); getline(fileExtin, token); getline(fileExtin, token);
			getline(fileAbsorp, token); getline(fileAbsorp, token); getline(fileAbsorp, token);

			float HRRPUV, Temperature;
			VolumeRegion *vol = *it;
			if(dynamic_cast<BlackBodyVolume * > (vol) == NULL) {
				printf("OMG. VolumeGridDensity and BlackBodyVolume are all confused!\n");
				exit(-10);
			}
			BlackBodyVolume * pvgd = (BlackBodyVolume*)(*it);
			float *HRRPUVs = pvgd->GetHRRPUVs();
			Point *centerPoints = new Point[nz * ny * nx];
			Color3d zeroColor(0.,0.,0.);
			for(int z = 0 ; z < nz ; z++) {
				for(int y = 0 ; y < ny ; y++) {
					for(int x = 0 ; x < nx ; x++) {
						fileH >> HRRPUV;
						fileTemp >> Temperature;					// Kelvin temperature
						// Based on Planckian Locus Approximation, determine the colors.
						if(Temperature > threshold) {
							
							// Here get the light position and then add it to lightpos vectors.
							// color *= HRRPUV or color *= BoltzmanEnergy.
							Color3d color = fire.ComputePlanckLocus(Temperature);
							float power = HRRPUV * fLightScale;//fire.ComputeStefanBoltzmann(Temperature);
// 							Color3d color1 = color * (HRRPUV / 1000.f);
// 							Color3d color2 = color * power / 100000.;

							lightPos.push_back(pvgd ->GetPointInVolume(x, y, z));
							lightVoxels.push_back(pvgd ->GetVoxelInVolume(x, y, z));
							//color *= (power * scale);
							lightColor.push_back(color * power);
							// Check Temperature to determine colors.
						}

						// Save HRRPUV.
						HRRPUVs[z * ny * nx + y * nx + x] = HRRPUV;
						centerPoints[z * ny * nx + y * nx + x] = pvgd->GetPointInVolume(x, y, z);
						
					}
				}
			}

			double minGrad = 100000.f, maxGrad = -1.f;
			int nCount1 = 0;int nCount2 = 0;int nCount3 = 0;
			// Precompute HRRPUV Gradients.
			Vector * gradientVectors = pvgd->GetGradientVectors();
			float * gradients = pvgd->GetGradients();
			for(int z = 0 ; z < nz ; z++) {
				for(int y = 0 ; y < ny ; y++) {
					for(int x = 0 ; x < nx ; x++) {
						// Compute Gradient Vectors from HRRPUV.
						if(x == 0 || y == 0 || z == 0 ||
							x == nx - 1 || y == ny - 1 || z == nz - 1) {

								gradientVectors[z * ny * nx + y * nx + x] = Vector(0, 0, -1);
								gradients[z * ny * nx + y * nx + x] = 0.0f;
						}
						else {
							Vector dirVec[6];
							int index[6] = { z * ny * nx + y * nx + (x + 1), z * ny * nx + y * nx + (x - 1),
											z * ny * nx + (y + 1) * nx + x, z * ny * nx + (y - 1) * nx + x,
											(z + 1) * ny * nx + y * nx + x, (z - 1) * ny * nx + y * nx + x };
							int index_c = z * ny * nx + y * nx + x;
							Vector sumVec(0.f, 0.f, 0.f);
							for(int vi = 0 ; vi < 6 ; vi++) {
								dirVec[vi] = centerPoints[index[vi]] - centerPoints[index_c];
								float cHRRPUVgradient = HRRPUVs[index[vi]] - HRRPUVs[index_c];
								sumVec += (dirVec[vi] * cHRRPUVgradient);
							}
							float length = sumVec.LengthSquared();
							if(length > maxGrad) maxGrad = length;
							if(length < minGrad) minGrad = length;
							if(length > 0. && length < 1000.) nCount1++;
							if(length > 1000. && length < 2000.) nCount2 ++;
							if(length > 2000. && length < 3300.) nCount3 ++;
							gradients[z * ny * nx + y * nx + x] = length;
							if(length != 0.f) {
								float oneOverLength = 1.0 / sqrtf(length);
								sumVec.x *= oneOverLength; sumVec.y *= oneOverLength; sumVec.z *= oneOverLength;
							}
							gradientVectors[z * ny * nx + y * nx + x] = sumVec; // most of the times, gradient vec will be (0,0,0).
						}
					}
				}
			}

			// Also add light at volume center + (eye - volume Center) * dist(eye-volume center) + (0.0, 1.0, 0.0)
			int count = lightColor.size();
			if(count == 0) {
				cout << "No light is generated due to the temperature and etc.\n";
				exit(-14);
			}
			else if(count > MAX_LIT_VOXELS - 1) {
				cout << "Can't allow more than " << MAX_LIT_VOXELS <<" lights. because Blackbody volume has voxelization and datastructure only for the lights.\n";
				exit(-19);
			}
		}
		// Here traverse lightpos vector to create multiple VPLs.
		// First, update light2world, paramSet. -- light2world is just the Identity Matrix. I guess at this point, light2world does not seem to matter.
		// paramSet
		if(renderOptions->volumeRegions.size() != 1) { printf("One volume is allowed Taekyu.\n"); exit(-32); }
		VolumeRegion * pvgd = renderOptions->volumeRegions.at(0);
		if(name == "point") {
			light = CreateFirePointLights(pvgd, light2world, paramSet, lightPos, lightColor, multipleLights);
		}

	}
    else if (name == "spot")
        light = CreateSpotLight(light2world, paramSet);
    else if (name == "goniometric")
        light = CreateGoniometricLight(light2world, paramSet);
    else if (name == "projection")
        light = CreateProjectionLight(light2world, paramSet);
    else if (name == "distant")
        light = CreateDistantLight(light2world, paramSet);
    else if (name == "infinite" || name == "exinfinite")
        light = CreateInfiniteLight(light2world, paramSet);
    else
        Warning("Light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return light;
}


Light *MakeLight(const string &name,
        const Transform &light2world, const ParamSet &paramSet, vector<Light*> *multipleLights) {
    Light *light = NULL;
    if (name == "point") {
		// Taekyu Shin
		vector<Point> lightPos;
		vector<VolumeRegion *>::iterator it;
        light = CreatePointLight(light2world, paramSet, lightPos, multipleLights);
	}
	else if (name == "voxel") {
		vector<BBox> lightBoxes;
		vector<VolumeRegion *>::iterator it;
        light = CreateVoxelLight(light2world, paramSet, lightBoxes, multipleLights);
	}
    else if (name == "spot")
        light = CreateSpotLight(light2world, paramSet);
    else if (name == "goniometric")
        light = CreateGoniometricLight(light2world, paramSet);
    else if (name == "projection")
        light = CreateProjectionLight(light2world, paramSet);
    else if (name == "distant")
        light = CreateDistantLight(light2world, paramSet);
    else if (name == "infinite" || name == "exinfinite")
        light = CreateInfiniteLight(light2world, paramSet);
    else
        Warning("Light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return light;
}


AreaLight *MakeAreaLight(const string &name,
        const Transform &light2world, const ParamSet &paramSet,
        const Reference<Shape> &shape) {
     AreaLight *area = NULL;
    if (name == "area" || name == "diffuse")
        area = CreateDiffuseAreaLight(light2world, paramSet, shape);
    else
        Warning("Area light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return area;
}

VolumeRegion *MakeVolumeRegion(const string &name,
        const Transform &volume2world, const ParamSet &paramSet, string & absorpFilename, string & extinctFilename) {
    VolumeRegion *vr = NULL;
    if (name == "homogeneous")
        vr = CreateHomogeneousVolumeDensityRegion(volume2world, paramSet);
    else if (name == "volumegrid") {
        vr = CreateGridVolumeRegion(volume2world, paramSet);
	}
    else if (name == "exponential")
        vr = CreateExponentialVolumeRegion(volume2world, paramSet);
	else if (name == "blackbody")
		vr = CreateBlackbodyGridVolumeRegion(volume2world, paramSet, absorpFilename, extinctFilename);
    else
        Warning("Volume region \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return vr;
}


SurfaceIntegrator *MakeSurfaceIntegrator(const string &name,
        const ParamSet &paramSet) {
    SurfaceIntegrator *si = NULL;
    if (name == "whitted")
        si = CreateWhittedSurfaceIntegrator(paramSet);
    else if (name == "directlighting")
        si = CreateDirectLightingIntegrator(paramSet);
    else if (name == "path")
        si = CreatePathSurfaceIntegrator(paramSet);
    else if (name == "photonmap" || name == "exphotonmap")
        si = CreatePhotonMapSurfaceIntegrator(paramSet);
    else if (name == "irradiancecache")
        si = CreateIrradianceCacheIntegrator(paramSet);
    else if (name == "igi")
        si = CreateIGISurfaceIntegrator(paramSet);
    else if (name == "dipolesubsurface")
        si = CreateDipoleSubsurfaceIntegrator(paramSet);
    else if (name == "ambientocclusion")
        si = CreateAmbientOcclusionIntegrator(paramSet);
    else if (name == "useprobes")
        si = CreateRadianceProbesSurfaceIntegrator(paramSet);
    else if (name == "diffuseprt")
        si = CreateDiffusePRTIntegratorSurfaceIntegrator(paramSet);
    else if (name == "glossyprt")
        si = CreateGlossyPRTIntegratorSurfaceIntegrator(paramSet);
    else
        Warning("Surface integrator \"%s\" unknown.", name.c_str());

    paramSet.ReportUnused();
    return si;
}


VolumeIntegrator *MakeVolumeIntegrator(const string &name,
        const ParamSet &paramSet) {
    VolumeIntegrator *vi = NULL;
    if (name == "single")
        vi = CreateSingleScatteringIntegrator(paramSet);
	else if(name == "premozemultiple")
		vi = CreateMultipleScatteringIntegrator(paramSet);
    else if (name == "emission")
        vi = CreateEmissionVolumeIntegrator(paramSet);
    else if (name == "blackbodymultiple")
		vi = CreateBlackBodyMultiScattering(paramSet);
	else if (name == "forwardmultiple")
		vi = CreateForwardMultipleScattering(paramSet);
	else {
        Warning("Volume integrator \"%s\" unknown.", name.c_str());
		exit(-10);
	}
    paramSet.ReportUnused();
    return vi;
}


Primitive *MakeAccelerator(const string &name,
        const vector<Reference<Primitive> > &prims,
        const ParamSet &paramSet) {
    Primitive *accel = NULL;
    if (name == "bvh")
        accel = CreateBVHAccelerator(prims, paramSet);
    else if (name == "grid")
        accel = CreateGridAccelerator(prims, paramSet);
    else if (name == "kdtree")
        accel = CreateKdTreeAccelerator(prims, paramSet);
    else
        Warning("Accelerator \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return accel;
}


Camera *MakeCamera(const string &name,
        const ParamSet &paramSet,
        const TransformSet &cam2worldSet, float transformStart,
        float transformEnd, Film *film) {
    Camera *camera = NULL;
    Assert(MAX_TRANSFORMS == 2);
    Transform *cam2world[2];
    transformCache.Lookup(cam2worldSet[0], &cam2world[0], NULL);
    transformCache.Lookup(cam2worldSet[1], &cam2world[1], NULL);
    AnimatedTransform animatedCam2World(cam2world[0], transformStart,
        cam2world[1], transformEnd);
    if (name == "perspective")
        camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film);
    else if (name == "orthographic")
        camera = CreateOrthographicCamera(paramSet, animatedCam2World, film);
    else if (name == "environment")
        camera = CreateEnvironmentCamera(paramSet, animatedCam2World, film);
    else
        Warning("Camera \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return camera;
}


Sampler *MakeSampler(const string &name,
        const ParamSet &paramSet, const Film *film, const Camera *camera) {
    Sampler *sampler = NULL;
    if (name == "adaptive")
        sampler = CreateAdaptiveSampler(paramSet, film, camera);
    else if (name == "bestcandidate")
        sampler = CreateBestCandidateSampler(paramSet, film, camera);
    else if (name == "halton")
        sampler = CreateHaltonSampler(paramSet, film, camera);
    else if (name == "lowdiscrepancy")
        sampler = CreateLowDiscrepancySampler(paramSet, film, camera);
    else if (name == "random")
        sampler = CreateRandomSampler(paramSet, film, camera);
    else if (name == "stratified")
        sampler = CreateStratifiedSampler(paramSet, film, camera);
    else
        Warning("Sampler \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return sampler;
}


Filter *MakeFilter(const string &name,
    const ParamSet &paramSet) {
    Filter *filter = NULL;
    if (name == "box")
        filter = CreateBoxFilter(paramSet);
    else if (name == "gaussian")
        filter = CreateGaussianFilter(paramSet);
    else if (name == "mitchell")
        filter = CreateMitchellFilter(paramSet);
    else if (name == "sinc")
        filter = CreateSincFilter(paramSet);
    else if (name == "triangle")
        filter = CreateTriangleFilter(paramSet);
    else
        Warning("Filter \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return filter;
}


Film *MakeFilm(const string &name,
    const ParamSet &paramSet, Filter *filter) {
    Film *film = NULL;
    if (name == "image")
        film = CreateImageFilm(paramSet, filter);
    else
        Warning("Film \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return film;
}



// API Function Definitions
void pbrtInit(const Options &opt) {
    PbrtOptions = opt;
    // API Initialization
    if (currentApiState != STATE_UNINITIALIZED)
        Error("pbrtInit() has already been called.");
    currentApiState = STATE_OPTIONS_BLOCK;
    renderOptions = new RenderOptions;
    graphicsState = GraphicsState();
    SampledSpectrum::Init();
}


void pbrtCleanup() {
    ProbesCleanup();
    // API Cleanup
    if (currentApiState == STATE_UNINITIALIZED)
        Error("pbrtCleanup() called without pbrtInit().");
    else if (currentApiState == STATE_WORLD_BLOCK)
        Error("pbrtCleanup() called while inside world block.");
    currentApiState = STATE_UNINITIALIZED;
    delete renderOptions;
    renderOptions = NULL;
}


void pbrtIdentity() {
    VERIFY_INITIALIZED("Identity");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform();)
}


void pbrtTranslate(float dx, float dy, float dz) {
    VERIFY_INITIALIZED("Translate");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] =
        curTransform[i] * Translate(Vector(dx, dy, dz));)
}


void pbrtTransform(float tr[16]) {
    VERIFY_INITIALIZED("Transform");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform(Matrix4x4(
        tr[0], tr[4], tr[8], tr[12],
        tr[1], tr[5], tr[9], tr[13],
        tr[2], tr[6], tr[10], tr[14],
        tr[3], tr[7], tr[11], tr[15]));)
}


void pbrtConcatTransform(float tr[16]) {
    VERIFY_INITIALIZED("ConcatTransform");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Transform(
                Matrix4x4(tr[0], tr[4], tr[8], tr[12],
                          tr[1], tr[5], tr[9], tr[13],
                          tr[2], tr[6], tr[10], tr[14],
                          tr[3], tr[7], tr[11], tr[15]));)
}


void pbrtRotate(float angle, float dx, float dy, float dz) {
    VERIFY_INITIALIZED("Rotate");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Rotate(angle, Vector(dx, dy, dz));)
}


void pbrtScale(float sx, float sy, float sz) {
    VERIFY_INITIALIZED("Scale");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Scale(sx, sy, sz);)
}


void pbrtLookAt(float ex, float ey, float ez, float lx, float ly, // This is where to set eye position(ex,ey,ez) and lookAt position(lx,ly,lz) then of course camera orientation vector.
        float lz, float ux, float uy, float uz) {
    VERIFY_INITIALIZED("LookAt");
    FOR_ACTIVE_TRANSFORMS({ Warning("This version of pbrt fixes a bug in the LookAt transformation.\n"
                                    "If your rendered images unexpectedly change, add a \"Scale -1 1 1\"\n"
                                    "to the start of your scene file."); break; })
    FOR_ACTIVE_TRANSFORMS(curTransform[i] =
        curTransform[i] * LookAt(Point(ex, ey, ez), Point(lx, ly, lz), Vector(ux, uy, uz));)
}


void pbrtCoordinateSystem(const string &name) {
    VERIFY_INITIALIZED("CoordinateSystem");
    namedCoordinateSystems[name] = curTransform;
}


void pbrtCoordSysTransform(const string &name) {
    VERIFY_INITIALIZED("CoordSysTransform");
    if (namedCoordinateSystems.find(name) !=
        namedCoordinateSystems.end())
        curTransform = namedCoordinateSystems[name];
    else
        Warning("Couldn't find named coordinate system \"%s\"",
                name.c_str());
}


void pbrtActiveTransformAll() {
    activeTransformBits = ALL_TRANSFORMS_BITS;
}


void pbrtActiveTransformEndTime() {
    activeTransformBits = END_TRANSFORM_BITS;
}


void pbrtActiveTransformStartTime() {
    activeTransformBits = START_TRANSFORM_BITS;
}


void pbrtTransformTimes(float start, float end) {
    VERIFY_OPTIONS("TransformTimes");
    renderOptions->transformStartTime = start;
    renderOptions->transformEndTime = end;
}


void pbrtPixelFilter(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("PixelFilter");
    renderOptions->FilterName = name;
    renderOptions->FilterParams = params;
}


void pbrtFilm(const string &type, const ParamSet &params) {
    VERIFY_OPTIONS("Film");
    renderOptions->FilmParams = params;
    renderOptions->FilmName = type;
}


void pbrtSampler(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Sampler");
    renderOptions->SamplerName = name;
    renderOptions->SamplerParams = params;
}


void pbrtAccelerator(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Accelerator");
    renderOptions->AcceleratorName = name;
    renderOptions->AcceleratorParams = params;
}


void pbrtSurfaceIntegrator(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("SurfaceIntegrator");
    renderOptions->SurfIntegratorName = name;
    renderOptions->SurfIntegratorParams = params;
}


void pbrtVolumeIntegrator(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("VolumeIntegrator");
    renderOptions->VolIntegratorName = name;
    renderOptions->VolIntegratorParams = params;
}


void pbrtRenderer(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Renderer");
    renderOptions->RendererName = name;
    renderOptions->RendererParams = params;
}


void pbrtCamera(const string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Camera");
    renderOptions->CameraName = name;
    renderOptions->CameraParams = params;
    renderOptions->CameraToWorld = Inverse(curTransform);
    namedCoordinateSystems["camera"] = renderOptions->CameraToWorld;
}


void pbrtWorldBegin() {
    VERIFY_OPTIONS("WorldBegin");
    currentApiState = STATE_WORLD_BLOCK;
    for (int i = 0; i < MAX_TRANSFORMS; ++i)
        curTransform[i] = Transform();
    activeTransformBits = ALL_TRANSFORMS_BITS;
    namedCoordinateSystems["world"] = curTransform;
}


void pbrtAttributeBegin() {
    VERIFY_WORLD("AttributeBegin");
    pushedGraphicsStates.push_back(graphicsState);
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
}


void pbrtAttributeEnd() {
    VERIFY_WORLD("AttributeEnd");
    if (!pushedGraphicsStates.size()) {
        Error("Unmatched pbrtAttributeEnd() encountered. "
              "Ignoring it.");
        return;
    }
    graphicsState = pushedGraphicsStates.back();
    pushedGraphicsStates.pop_back();
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
}


void pbrtTransformBegin() {
    VERIFY_WORLD("TransformBegin");
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
}


void pbrtTransformEnd() {
    VERIFY_WORLD("TransformEnd");
    if (!pushedTransforms.size()) {
        Error("Unmatched pbrtTransformEnd() encountered. "
            "Ignoring it.");
        return;
    }
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
}


void pbrtTexture(const string &name, const string &type,
                 const string &texname, const ParamSet &params) {
    VERIFY_WORLD("Texture");
    TextureParams tp(params, params, graphicsState.floatTextures,
                     graphicsState.spectrumTextures);
    if (type == "float")  {
        // Create _float_ texture and store in _floatTextures_
        if (graphicsState.floatTextures.find(name) !=
            graphicsState.floatTextures.end())
            Info("Texture \"%s\" being redefined", name.c_str());
        WARN_IF_ANIMATED_TRANSFORM("Texture");
        Reference<Texture<float> > ft = MakeFloatTexture(texname,
                                                         curTransform[0], tp);
        if (ft) graphicsState.floatTextures[name] = ft;
    }
    else if (type == "color" || type == "spectrum")  {
        // Create _color_ texture and store in _spectrumTextures_
        if (graphicsState.spectrumTextures.find(name) != graphicsState.spectrumTextures.end())
            Info("Texture \"%s\" being redefined", name.c_str());
        WARN_IF_ANIMATED_TRANSFORM("Texture");
        Reference<Texture<Spectrum> > st = MakeSpectrumTexture(texname,
            curTransform[0], tp);
        if (st) graphicsState.spectrumTextures[name] = st;
    }
    else
        Error("Texture type \"%s\" unknown.", type.c_str());
}


void pbrtMaterial(const string &name, const ParamSet &params) {
    VERIFY_WORLD("Material");
    graphicsState.material = name;
    graphicsState.materialParams = params;
    graphicsState.currentNamedMaterial = "";
}


void pbrtMakeNamedMaterial(const string &name,
        const ParamSet &params) {
    VERIFY_WORLD("MakeNamedMaterial");
    // error checking, warning if replace, what to use for transform?
    TextureParams mp(params, graphicsState.materialParams,
                     graphicsState.floatTextures,
                     graphicsState.spectrumTextures);
    string matName = mp.FindString("type");
    WARN_IF_ANIMATED_TRANSFORM("MakeNamedMaterial");
    if (matName == "") Error("No parameter string \"type\" found in MakeNamedMaterial");
    else {
        Reference<Material> mtl = MakeMaterial(matName, curTransform[0], mp);
        if (mtl) graphicsState.namedMaterials[name] = mtl;
    }
}



void pbrtNamedMaterial(const string &name) {
    VERIFY_WORLD("NamedMaterial");
    graphicsState.currentNamedMaterial = name;
}

/*
15,6,9 = (z,y,x) ==> (0.9381,0.381, 0.566)
15,20,27 ==> (.9381, 1.2475, 1.68)
26,20,10 --> (1.618, 1.2475, 0628)
These are the brightest voxel coordinates in the volume. 3 points.
*/

Light * gSunLight = NULL;
extern float fVoxelSize;
void pbrtLightSource(const string &name, const ParamSet &params) {
    VERIFY_WORLD("LightSource");
    WARN_IF_ANIMATED_TRANSFORM("LightSource");
	vector<Light *> multipleFireLights;
	bool bFireLight = params.FindOneBool("FireLight", false);
	bool bAutomaticLightPlacement = params.FindOneBool("AutomaticLightPlacement", false);
	Light *lt;
	if(bAutomaticLightPlacement == true) {
		float threshold = params.FindOneFloat("thresholdForLightPlacement", 1667.f);
		std::string fileHRRPUV, fileTemperature;
		std::string fileAbsorption, fileExtinction;
		fileHRRPUV = params.FindOneFilename("hrrpuvFile", "SmokeHRRPUV.pbrt");
		if(fileHRRPUV.substr(fileHRRPUV.length() - 8, 3) == "200") {
			fVoxelSize = 0.0099004975124378f;
		}
		else {
			fVoxelSize = 0.0245679012345679f;
		}
		fileTemperature = params.FindOneFilename("temperatureFile", "SmokeKelvinTemperature.pbrt");

		
		lt = MakeLightAutomaticPlacement(name, curTransform[0], params, fileTemperature, fileHRRPUV, &multipleFireLights, threshold);
	}
	else
		lt = MakeLight(name, curTransform[0], params, &multipleFireLights);
    if (lt == NULL)
        Error("pbrtLightSource: light type \"%s\" unknown.", name.c_str());
    else
	{
		if(multipleFireLights.size() == 0) {
			renderOptions->lights.push_back(lt);
			renderOptions->IsLightsFireKind.push_back(bFireLight);
		}
		else {
// 			renderOptions->lights.push_back(lt); // In case of VoxelLight, lights[0] is still a point light.
// 			renderOptions->IsLightsFireKind.push_back(bFireLight);
			gSunLight = lt;
			PointLight * pl = dynamic_cast<PointLight* > (lt);
			if(pl != NULL)
				pl->SetSampleID(0);

			vector<Light *>::iterator it;
			for(it = multipleFireLights.begin() ; it < multipleFireLights.end() ; it++) {
				renderOptions ->lights.push_back((*it));
				renderOptions ->IsLightsFireKind.push_back(true); 
			}
		}
	}
}


void pbrtAreaLightSource(const string &name,
                         const ParamSet &params) {
    VERIFY_WORLD("AreaLightSource");
    graphicsState.areaLight = name;
    graphicsState.areaLightParams = params;
}


void pbrtShape(const string &name, const ParamSet &params) {
    VERIFY_WORLD("Shape");
    Reference<Primitive> prim;
    AreaLight *area = NULL;
    if (!curTransform.IsAnimated()) {
        // Create primitive for static shape
        Transform *obj2world, *world2obj;
        transformCache.Lookup(curTransform[0], &obj2world, &world2obj);
        Reference<Shape> shape = MakeShape(name, obj2world, world2obj,
            graphicsState.reverseOrientation, params);
        if (!shape) return;
        Reference<Material> mtl = graphicsState.CreateMaterial(params);
        params.ReportUnused();

        // Possibly create area light for shape
        if (graphicsState.areaLight != "") {
            area = MakeAreaLight(graphicsState.areaLight, curTransform[0],
                                 graphicsState.areaLightParams, shape);
        }
        prim = new GeometricPrimitive(shape, mtl, area);
    } else {
        // Create primitive for animated shape

        // Create initial _Shape_ for animated shape
        if (graphicsState.areaLight != "")
            Warning("Ignoring currently set area light when creating "
                    "animated shape");
        Transform *identity;
        transformCache.Lookup(Transform(), &identity, NULL);
        Reference<Shape> shape = MakeShape(name, identity, identity,
            graphicsState.reverseOrientation, params);
        if (!shape) return;
        Reference<Material> mtl = graphicsState.CreateMaterial(params);
        params.ReportUnused();

        // Get _animatedWorldToObject_ transform for shape
        Assert(MAX_TRANSFORMS == 2);
        Transform *world2obj[2];
        transformCache.Lookup(curTransform[0], NULL, &world2obj[0]);
        transformCache.Lookup(curTransform[1], NULL, &world2obj[1]);
        AnimatedTransform
             animatedWorldToObject(world2obj[0], renderOptions->transformStartTime,
                                   world2obj[1], renderOptions->transformEndTime);
        Reference<Primitive> baseprim = new GeometricPrimitive(shape, mtl, NULL);
        if (!baseprim->CanIntersect()) {
            // Refine animated shape and create BVH if more than one shape created
            vector<Reference<Primitive> > refinedPrimitives;
            baseprim->FullyRefine(refinedPrimitives);
            if (refinedPrimitives.size() == 0) return;
            if (refinedPrimitives.size() > 1)
                baseprim = new BVHAccel(refinedPrimitives);
            else
                baseprim = refinedPrimitives[0];
        }
        prim = new TransformedPrimitive(baseprim, animatedWorldToObject);
    }
    // Add primitive to scene or current instance
    if (renderOptions->currentInstance) {
        if (area)
            Warning("Area lights not supported with object instancing");
        renderOptions->currentInstance->push_back(prim);
    }
    
    else {
        renderOptions->primitives.push_back(prim);
        if (area != NULL) {
            renderOptions->lights.push_back(area);
        }
    }
}


Reference<Material> GraphicsState::CreateMaterial(const ParamSet &params) {
    TextureParams mp(params, materialParams,
                     floatTextures,
                     spectrumTextures);
    Reference<Material> mtl;
    if (currentNamedMaterial != "" &&
        namedMaterials.find(currentNamedMaterial) != namedMaterials.end())
        mtl = namedMaterials[graphicsState.currentNamedMaterial];
    if (!mtl)
        mtl = MakeMaterial(material, curTransform[0], mp);
    if (!mtl)
        mtl = MakeMaterial("matte", curTransform[0], mp);
    if (!mtl)
        Severe("Unable to create \"matte\" material?!");
    return mtl;
}


void pbrtReverseOrientation() {
    VERIFY_WORLD("ReverseOrientation");
    graphicsState.reverseOrientation =
        !graphicsState.reverseOrientation;
}


void pbrtVolume(const string &name, const ParamSet &params) {
    VERIFY_WORLD("Volume");
    WARN_IF_ANIMATED_TRANSFORM("Volume");

	std::string fileAbsorption = params.FindOneFilename("absorptionFile", "");
	std::string fileExtinction = params.FindOneFilename("extinctionFile", "");
	if(fileAbsorption == "" || fileExtinction == "")  {
		printf("Volume file names are required.\n");
		exit(-1);
	}

    VolumeRegion *vr = MakeVolumeRegion(name, curTransform[0], params, fileAbsorption, fileExtinction);
    if (vr) renderOptions->volumeRegions.push_back(vr); // volumeRegions == the type of VolumeRegion. DensityRegion derives from the class.
	// So Later, I need to check this data to see where to place lights.
}


void pbrtObjectBegin(const string &name) {
    VERIFY_WORLD("ObjectBegin");
    pbrtAttributeBegin();
    if (renderOptions->currentInstance)
        Error("ObjectBegin called inside of instance definition");
    renderOptions->instances[name] = vector<Reference<Primitive> >();
    renderOptions->currentInstance = &renderOptions->instances[name];
}


void pbrtObjectEnd() {
    VERIFY_WORLD("ObjectEnd");
    if (!renderOptions->currentInstance)
        Error("ObjectEnd called outside of instance definition");
    renderOptions->currentInstance = NULL;
    pbrtAttributeEnd();
}


void pbrtObjectInstance(const string &name) {
    VERIFY_WORLD("ObjectInstance");
    // Object instance error checking
    if (renderOptions->currentInstance) {
        Error("ObjectInstance can't be called inside instance definition");
        return;
    }
    if (renderOptions->instances.find(name) == renderOptions->instances.end()) {
        Error("Unable to find instance named \"%s\"", name.c_str());
        return;
    }
    vector<Reference<Primitive> > &in = renderOptions->instances[name];
    if (in.size() == 0) return;
    if (in.size() > 1 || !in[0]->CanIntersect()) {
        // Refine instance _Primitive_s and create aggregate
        Reference<Primitive> accel =
             MakeAccelerator(renderOptions->AcceleratorName,
                             in, renderOptions->AcceleratorParams);
        if (!accel) accel = MakeAccelerator("bvh", in, ParamSet());
        if (!accel) Severe("Unable to create \"bvh\" accelerator");
        in.erase(in.begin(), in.end());
        in.push_back(accel);
    }
    Assert(MAX_TRANSFORMS == 2);
    Transform *world2instance[2];
    transformCache.Lookup(curTransform[0], NULL, &world2instance[0]);
    transformCache.Lookup(curTransform[1], NULL, &world2instance[1]);
    AnimatedTransform animatedWorldToInstance(world2instance[0],
        renderOptions->transformStartTime,
        world2instance[1], renderOptions->transformEndTime);
    Reference<Primitive> prim =
        new TransformedPrimitive(in[0], animatedWorldToInstance);
    renderOptions->primitives.push_back(prim);
}


void pbrtWorldEnd() {
    VERIFY_WORLD("WorldEnd");
    // Ensure there are no pushed graphics states
    while (pushedGraphicsStates.size()) {
        Warning("Missing end to pbrtAttributeBegin()");
        pushedGraphicsStates.pop_back();
        pushedTransforms.pop_back();
    }
    while (pushedTransforms.size()) {
        Warning("Missing end to pbrtTransformBegin()");
        pushedTransforms.pop_back();
    }

    // Create scene and render
    Renderer *renderer = renderOptions->MakeRenderer();
    Scene *scene = renderOptions->MakeScene();
    if (scene && renderer) renderer->Render(scene); // Most Important Part! Here the rendering happens, but the Integrator will add up all the samples. Actual rendering may happen in SampleRenderTask
    TasksCleanup();
    delete renderer;
    delete scene;

    // Clean up after rendering
    graphicsState = GraphicsState();
    transformCache.Clear();
    currentApiState = STATE_OPTIONS_BLOCK;
    ProbesPrint(stdout);
    for (int i = 0; i < MAX_TRANSFORMS; ++i)
        curTransform[i] = Transform();
    activeTransformBits = ALL_TRANSFORMS_BITS;
    namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
                                 namedCoordinateSystems.end());
    ImageTexture<float, float>::ClearCache();
    ImageTexture<RGBSpectrum, Spectrum>::ClearCache();
}


Scene *RenderOptions::MakeScene() {
    // Initialize _volumeRegion_ from volume region(s)
    VolumeRegion *volumeRegion;
    if (volumeRegions.size() == 0)
        volumeRegion = NULL;
    else if (volumeRegions.size() == 1) {
        volumeRegion = volumeRegions[0];
	}
    else
        volumeRegion = new AggregateVolume(volumeRegions);
    Primitive *accelerator = MakeAccelerator(AcceleratorName,
        primitives, AcceleratorParams);
    if (!accelerator)
        accelerator = MakeAccelerator("bvh", primitives, ParamSet());
    if (!accelerator)
        Severe("Unable to create \"bvh\" accelerator.");
	Scene *scene = new Scene(accelerator, lights, IsLightsFireKind, volumeRegion);
    // Erase primitives, lights, and volume regions from _RenderOptions_
    primitives.erase(primitives.begin(), primitives.end());
    lights.erase(lights.begin(), lights.end());
    volumeRegions.erase(volumeRegions.begin(), volumeRegions.end());
	if( volumeRegion != NULL && dynamic_cast<VolumeGridDensity *> (volumeRegion) ) {
		((VolumeGridDensity *)volumeRegion)->mScene = scene;
	}
	else {
		if( volumeRegion != NULL && dynamic_cast<BlackBodyVolume *> (volumeRegion) ) {
			((BlackBodyVolume *)volumeRegion)->mScene = scene;
		}
	}
    return scene;
}


Renderer *RenderOptions::MakeRenderer() const {
    Renderer *renderer = NULL;
    Camera *camera = MakeCamera();
    if (RendererName == "metropolis") {
		// From here, creation of a samplerRenderer.
		bool visIds = RendererParams.FindOneBool("visualizeobjectids", false);
        RendererParams.ReportUnused();
        Sampler *sampler = MakeSampler(SamplerName, SamplerParams, camera->film, camera);
        if (!sampler) Severe("Unable to create sampler.");
        // Create surface and volume integrators
        SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
            SurfIntegratorParams);
        if (!surfaceIntegrator) Severe("Unable to create surface integrator.");
        VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
            VolIntegratorParams);
        if (!volumeIntegrator) Severe("Unable to create volume integrator.");
        SamplerRenderer *sampleRenderer = new SamplerRenderer(sampler, camera, surfaceIntegrator, // add metropolis
                                       volumeIntegrator, visIds, NULL); // pMetroMediaRenderer is always NULL.
        // Warn if no light sources are defined
        if (lights.size() == 0)
            Warning("No light sources defined in scene; "
                "possibly rendering a black image.");

		// Metropolis renderer creation
		renderer = CreateMetropolisRenderer(RendererParams, camera, sampleRenderer);
        RendererParams.ReportUnused();
        // Warn if no light sources are defined
        if (lights.size() == 0)
            Warning("No light sources defined in scene; "
                "possibly rendering a black image.");
    }
    // Create remaining _Renderer_ types
    else if (RendererName == "createprobes") {
        // Create surface and volume integrators
        SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
            SurfIntegratorParams);
        if (!surfaceIntegrator) Severe("Unable to create surface integrator.");
        VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
            VolIntegratorParams);
        if (!volumeIntegrator) Severe("Unable to create volume integrator.");
        renderer = CreateRadianceProbesRenderer(camera, surfaceIntegrator, volumeIntegrator, RendererParams);
        RendererParams.ReportUnused();
        // Warn if no light sources are defined
        if (lights.size() == 0)
            Warning("No light sources defined in scene; "
                "possibly rendering a black image.");
    }
    else if (RendererName == "aggregatetest") {
        renderer = CreateAggregateTestRenderer(RendererParams, primitives);
        RendererParams.ReportUnused();
    }
    else if (RendererName == "surfacepoints") {
        Point pCamera = camera->CameraToWorld(camera->shutterOpen, Point(0, 0, 0));
        renderer = CreateSurfacePointsRenderer(RendererParams, pCamera, camera->shutterOpen);
        RendererParams.ReportUnused();
    }
    else { // renderer Name == sampler
        if (RendererName != "sampler" && RendererName != "metropolisMedia")
            Warning("Renderer type \"%s\" unknown.  Using \"sampler\".",
                    RendererName.c_str());
        bool visIds = RendererParams.FindOneBool("visualizeobjectids", false);
        RendererParams.ReportUnused();
		// Taekyu - sampler is responsible for choosing the points on the image plane from which rays are traced.
        Sampler *sampler = MakeSampler(SamplerName, SamplerParams, camera->film, camera);
        if (!sampler) Severe("Unable to create sampler.");
        // Create surface and volume integrators
        SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
            SurfIntegratorParams);
        if (!surfaceIntegrator) Severe("Unable to create surface integrator.");
        VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
            VolIntegratorParams);
        if (!volumeIntegrator) Severe("Unable to create volume integrator.");
		float nRaymarchingEyeDistanceScale = RendererParams.FindOneFloat("raymarchingEyeDistanceScale", 10.0f);
		float nRaymarchingLightDistanceScale = RendererParams.FindOneFloat("raymarchingLightDistanceScale", 40.f);
		renderer = new SamplerRenderer(sampler, camera, surfaceIntegrator, // add metropolis
                                       volumeIntegrator, visIds, NULL); // pMetroMediaRenderer is always NULL.
        // Warn if no light sources are defined
        if (lights.size() == 0)
            Warning("No light sources defined in scene; "
                "possibly rendering a black image.");
    }
    return renderer;
}


Camera *RenderOptions::MakeCamera() const {
    Filter *filter = MakeFilter(FilterName, FilterParams);
    Film *film = MakeFilm(FilmName, FilmParams, filter);
    if (!film) Severe("Unable to create film.");
    Camera *camera = ::MakeCamera(CameraName, CameraParams,
        CameraToWorld, renderOptions->transformStartTime,
        renderOptions->transformEndTime, film);
    if (!camera) Severe("Unable to create camera.");
    return camera;
}


