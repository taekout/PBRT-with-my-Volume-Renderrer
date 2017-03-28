#include "SobelGenerator.h"

#include "stdafx.h"

using namespace std;

SobelGenerator::SobelGenerator(void)
{
        outDistr = new OutgoingDistribution();
}


SobelGenerator::~SobelGenerator(void)
{
        if(outDistr != NULL) delete outDistr;
}

/*
** z_basis is the scatter vector direction. x and y are generated with cross product.
(x_basis, y_basis, z_basis) * (0,0,1)' = (IdentityMatrix) * (SobelSequence)';
*/
void SobelGenerator::TransformDirectionIntoCurrentRaySpace(Vector &scatterDir, vector<Vector *> & largeRayDirectionInIdentitySpace,
                                                                        vector<Vector *> & resultRayDirectionInRaySpace)
{
        Vector arbitrary(0.f, 1.f, 0.f);
        Vector x = Cross(scatterDir, arbitrary);
        Vector y = Cross(x, scatterDir);
        vector<Vector *>::iterator it = largeRayDirectionInIdentitySpace.begin();
        vector<Vector *>::iterator itend = largeRayDirectionInIdentitySpace.end();
        float dirx, diry, dirz;
        for( ; it < itend ; it++) {
                Vector *vec = (*it);
                dirx = x.x * vec ->x + y.x * vec ->y + scatterDir.x * vec ->z;
                diry = x.y * vec ->x + y.y * vec ->y + scatterDir.y * vec ->z;
                dirz = x.z * vec ->x + y.z * vec ->y + scatterDir.z * vec ->z;
                resultRayDirectionInRaySpace.push_back(new Vector(dirx, diry, dirz));
        }
}

// Now I have random sobel sequences. I use rng to pick random n number of samples. and return it in nVectors.
// Also probability is computed based on PhaseHG.
// pdf[i] = phaseHG(...);

float SobelGenerator::PhaseHG(const Vector &w, const Vector &wp, float g) {
    float costheta = Dot(w, wp);
    return 1.f / (4.f * M_PI) *
        (1.f - g*g) / powf(1.f + g*g - 2.f * g * costheta, 1.5f);
}

void SobelGenerator::PickRandomPoints(const RNG &rng, int nSamples,
        vector<int> &VectorIndices, vector<Vector *> &largeRayDirections)
{
        for(int i = 0 ; i < nSamples ; i++) {
                unsigned int index = rng.RandomUInt() % 1024;
                VectorIndices.push_back(index);

                // Cannot compute PDFs here. PDF is based on the Wi and Wo. Thus I must compute it on the fly.
                Vector vec(outDistr ->points[index][0], outDistr ->points[index][1], outDistr ->points[index][2]);
                Vector *pVec = new Vector(vec);
                largeRayDirections.push_back(pVec);
        }
}

double ** SobelGenerator::sobol_points(unsigned N, unsigned D, char *dir_file)
{
  ifstream infile(dir_file,ios::in);
  if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
    exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 

  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [N];
  C[0] = 1;
  for (unsigned i=1;i<=N-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  double **POINTS = new double * [N];
  for (unsigned i=0;i<=N-1;i++) POINTS[i] = new double [D];
  for (unsigned j=0;j<=D-1;j++) POINTS[0][j] = 0; 

  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  unsigned *V = new unsigned [L+1]; 
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [N];
  X[0] = 0;
  for (unsigned i=1;i<=N-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    POINTS[i][0] = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }
  
  // Clean up
  delete [] V;
  delete [] X;
  
  
  // ----- Compute the remaining dimensions -----
  for (unsigned j=1;j<=D-1;j++) {
    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    unsigned *V = new unsigned [L+1];
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
        V[i] = V[i-s] ^ (V[i-s] >> s); 
        for (unsigned k=1;k<=s-1;k++) 
          V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [N];
    X[0] = 0;
    for (unsigned i=1;i<=N-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
      POINTS[i][j] = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
   }
    
    // Clean up
    delete [] m;
    delete [] V;
    delete [] X;
  }
  delete [] C;
  
  return POINTS;
}

void SobelGenerator::DistributedVectorGenerator(unsigned N, char *dir_file, float power) // power == 1 --> well distributed over the sphere.
{
        unsigned int D = 3;
        double **P = sobol_points(N, 2, dir_file);
        double **points;
        points = new double*[N];
        for(int i = 0 ; i < N ; i++)  {
                (points)[i] = new double[D];
                for(int j = 0;  j < D ; j++) {
                        (points)[i][j] = 0.0;
                }
        }
        float *pdfArray = new float[N];
        for(int i = 0 ; i < N ; i++) {
                (pdfArray)[i] = 0.f;
        }

        for (unsigned i=0;i<=N-1;i++) {
                        // P[i][j]. i = [0, 1023], j = [0,1];
                        /*
distributed by cos^n with
 angle = 2*pi*u
 z = pow(v, 1/(n+1))
 r = sqrt(1 - z*z)
 H = vec3(r*cos(a), r*sin(a), z)

From H, you can get a direction over the entire sphere by reflecting (0,0,1)
 D = 2*H.z*H - (0,0,1)

According to pbrt, in core/reflection.cpp, the probability of this sample is
 pdf = (n+1)*pow(z,n) / (8*pi*H.z)

For n=1, you get a distribution over the entire sphere. As n gets bigger, the distribution narrows dramatically, with few, if any, points far from the center.
                        */
                float angle = 2 * M_PI * P[i][0];
                float z = pow((float)P[i][1], 1/(power+1));
                float r = sqrt(1 - z * z);
                Vector *H = new Vector(r*cos(angle), r * sin(angle), z);
                *H *= H ->z * 2.0f;
                Vector D(0.f, 0.f, -1.f);
                D += *H;
                // D is the direction now save it in (*points)[i][j].
                (points)[i][0] = D.x;
                (points)[i][1] = D.y;
                (points)[i][2] = D.z;

                //compute pdf.
                float pdf = (power+1) * pow(z, power) / ( 8.f * M_PI * H ->z );
                (pdfArray)[i] = pdf;
        }
        this ->outDistr ->points = points;
        this ->outDistr ->pdf = pdfArray;
        delete [] P;
}