// Copyright (C) 2016 Thomas Schrefl, Dieter Suess, Florian Bruckner
// Last modified by Florian Bruckner, 2016-06-17

#include <math.h>
#include <sstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "basic.h"
#include "avector.h"
#include "surface3d.h"
#include "laplacebem3d.h"
#include "h2compression.h"


class Lindholm_C {
  public:
    Lindholm_C();
    int geometry_from_file(std::string infile);
    int geometry_from_array(unsigned int N, double coordinates[][3], unsigned int NE, int elements[][3]);
    int setup();
    int matvec(unsigned int N, double x[], double b[]);
    ~Lindholm_C();

  private:
    // variable definitions
    psurface3d gr;
    uint vertices;
    uint q_reg, q_sing;
    pbem3d bem;
    uint clf;
    pcluster root;
    real eta;
    pblock broot;
    phmatrix Kh;
    real eps_aca;
    real eps_recomp;
    ptruncmode tm;
    ph2matrix Kh2;
    pstopwatch sw;

    // method definitions
};


