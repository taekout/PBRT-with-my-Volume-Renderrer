#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_API_H
#define PBRT_CORE_API_H

// core/api.h*
#include "pbrt.h"

// API Function Declarations
void pbrtInit(const Options &opt);
void pbrtCleanup();
void pbrtIdentity();
void pbrtTranslate(float dx, float dy, float dz);
void pbrtRotate(float angle, float ax, float ay, float az);
void pbrtScale(float sx, float sy, float sz);
void pbrtLookAt(float ex, float ey, float ez,
                float lx, float ly, float lz,
                float ux, float uy, float uz);
void pbrtConcatTransform(float transform[16]);
void pbrtTransform(float transform[16]);
void pbrtCoordinateSystem(const string &);
void pbrtCoordSysTransform(const string &);
void pbrtActiveTransformAll();
void pbrtActiveTransformEndTime();
void pbrtActiveTransformStartTime();
void pbrtTransformTimes(float start, float end);
void pbrtPixelFilter(const string &name, const ParamSet &params);
void pbrtFilm(const string &type, const ParamSet &params);
void pbrtSampler(const string &name, const ParamSet &params);
void pbrtAccelerator(const string &name, const ParamSet &params);
void pbrtSurfaceIntegrator(const string &name, const ParamSet &params);
void pbrtVolumeIntegrator(const string &name, const ParamSet &params);
void pbrtRenderer(const string &name, const ParamSet &params);
void pbrtCamera(const string &, const ParamSet &cameraParams);
void pbrtWorldBegin();
void pbrtAttributeBegin();
void pbrtAttributeEnd();
void pbrtTransformBegin();
void pbrtTransformEnd();
void pbrtTexture(const string &name, const string &type,
    const string &texname, const ParamSet &params);
void pbrtMaterial(const string &name, const ParamSet &params);
void pbrtMakeNamedMaterial(const string &name, const ParamSet &params);
void pbrtNamedMaterial(const string &name);
void pbrtLightSource(const string &name, const ParamSet &params);
void pbrtAreaLightSource(const string &name, const ParamSet &params);
void pbrtShape(const string &name, const ParamSet &params);
void pbrtReverseOrientation();
void pbrtVolume(const string &name, const ParamSet &params);
void pbrtObjectBegin(const string &name);
void pbrtObjectEnd();
void pbrtObjectInstance(const string &name);
void pbrtWorldEnd();


#endif // PBRT_CORE_API_H
