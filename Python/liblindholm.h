// Copyright (C) 2016 Thomas Schrefl, Dieter Suess, Florian Bruckner
// Last modified by Florian Bruckner, 2016-06-17

#include <math.h>
#include <sstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif

extern "C" {
  #include "basic.h"
  #include "avector.h"
  #include "surface3d.h"
  #include "laplacebem3d.h"
  #include "h2compression.h"
}

class Lindholm_C {
  public:
    Lindholm_C();
    int geometry_from_file(std::string infile);
    int geometry_from_array(unsigned int N, double coordinates[][3], unsigned int NE, int elements[][3]);
    int setup_HCA();
    int setup_GCA();
    int matvec(unsigned int N, double x[], double b[]);
    int get_size();
    ~Lindholm_C();

  private:
    int _setup_root();
    psurface3d gr;
    uint vertices;
    uint q_reg, q_sing;
    pbem3d bem;
    uint clf;
    pcluster root;
    pclusterbasis rb, cb;
    real eta;
    pblock broot;
    real eps_aca;
    real eps_recomp;
    real delta;
    uint m;
    ptruncmode tm;
    ph2matrix Kh2;
};


