// header file for building precompiled headers under windows
// a no-op on other architectures
 
#if defined(_MSC_VER)

#include "pbrt.h"
#include "camera.h"
#include "scene.h"
#include "imageio.h"
#include "intersection.h"
#include "montecarlo.h"
#include "sampler.h"
#include "texture.h"
#include "integrator.h"
#include <iostream>
#include <fstream>

#define MAX_LIGHT	5
#define MAX_STEP	20

using namespace std;

#endif // _MSC_VER
