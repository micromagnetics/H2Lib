// Copyright (C) 2016 Thomas Schrefl, Dieter Suess, Florian Bruckner
// Last modified by Florian Bruckner, 2016-06-17

#include <math.h>
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
    ~Lindholm_C();
    
  private:
    psurface3d gr;
    uint vertices;
    uint q_reg, q_sing;
    pbem3d bem;
    uint clf;
    //pamatrix K;
    pcluster root;
    real eta;
    pblock broot;
    phmatrix Kh;
    real eps_aca;
    real eps_recomp;
    ptruncmode tm;
    ph2matrix Kh2;
    pstopwatch sw;
};

