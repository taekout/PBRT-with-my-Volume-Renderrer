
#include "pbrt.h"
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"
#include "paramset.h"
#include "pbrtparse.cpp"

void PrevMain()
{
	Options opt;
	pbrtInit(opt);
	ParamSet params;
    InitParamSet(params, SPECTRUM_REFLECTANCE);
	string vol = "volume";
	pbrtLookAt(3.3, 4, -4, 1.2, .7, 0, 0, 1, 0);
	FreeArgs();
    InitParamSet(params, SPECTRUM_REFLECTANCE);
    pbrtCamera("perspective", params);
    FreeArgs();
	//pbrtVolume(vol, params);
	InitParamSet(params, SPECTRUM_REFLECTANCE); // The problem is that it did not add properly all the more necessary information such as "integer xresolution" [800] "integer yresolution" [800]
												//"string filename" "smoke-2.exr"
    pbrtFilm("image", params);
    FreeArgs();

}