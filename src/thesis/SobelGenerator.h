#pragma once
#include <stdlib.h>
#include "geometry.h"

struct OutgoingDistribution{
        double **points;
        float *pdf;

        OutgoingDistribution(double **_points, float *_pdf) { this ->points = _points; pdf = _pdf;}
        OutgoingDistribution() { this ->points = NULL; pdf = NULL;}
        ~OutgoingDistribution() { if(points != NULL) delete [] points; if(pdf != NULL) delete [] pdf;}
} ;

class SobelGenerator
{
public:
        SobelGenerator(void);
        ~SobelGenerator(void);

        void DistributedVectorGenerator(unsigned N, char *dir_file, float power); // power == 1 --> well distributed over the sphere.
        double **sobol_points(unsigned N, unsigned D, char *dir_file);
        void TransformDirectionIntoCurrentRaySpace(Vector &z, vector<Vector *> & largeRayDirectionInIdentitySpace
        , vector<Vector *> & resultRayDirectionInRaySpace); 
        void PickRandomPoints(const RNG &rng, int nSamples, vector<int> &VectorIndices, vector<Vector *> &largeRayDirections);
        float PhaseHG(const Vector &w, const Vector &wp, float g);
        inline Vector GetPoint(int index) {
                Vector vec(outDistr ->points[index][0], outDistr ->points[index][1], outDistr ->points[index][2]);
                return vec;
        }
        inline double **GetAllPoints(void) {
                return outDistr ->points;
        }

private:
        OutgoingDistribution *outDistr ;
};