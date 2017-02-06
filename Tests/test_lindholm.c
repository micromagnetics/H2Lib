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

/* math utility functions */
static inline void cross(double a[3], double b[3], double axb[3]) {
  axb[0] = a[1] * b[2] - a[2] * b[1];
  axb[1] = a[2] * b[0] - a[0] * b[2];
  axb[2] = a[0] * b[1] - a[1] * b[0];
}

static inline double dot(double a[3], double b[3]) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

static inline double norm3(double a[3]) {
  return (sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]));
}

static inline void normalize(double a[3]) {
  double norm = norm3(a);
  a[0] /= norm;
  a[1] /= norm;
  a[2] /= norm;
}

static inline double solidangle_rech(double bro1, double bro2, double bro3,
    double ro1[3], double ro2[3], double ro3[3]) {
  int i;
  double nro1[3], nro2[3], nro3[3];
  double z1, z2, z3;
  double n1, n2, n3;
  double z, n;

  if ((bro1 <= 0.) || (bro2 <= 0.) || (bro3 <= 0.)) {
    return (0.);
  }

  for (i = 0; i < 3; i++) {
    nro1[i] = ro1[i] / bro1;
    nro2[i] = ro2[i] / bro2;
    nro3[i] = ro3[i] / bro3;
  }

  z1 = nro1[0] * (nro2[1] * nro3[2] - nro2[2] * nro3[1]);
  z2 = nro1[1] * (nro2[2] * nro3[0] - nro2[0] * nro3[2]);
  z3 = nro1[2] * (nro2[0] * nro3[1] - nro2[1] * nro3[0]);
  z = z1 + z2 + z3;

  n1 = nro2[0] * nro3[0] + nro2[1] * nro3[1] + nro2[2] * nro3[2];
  n2 = nro3[0] * nro1[0] + nro3[1] * nro1[1] + nro3[2] * nro1[2];
  n3 = nro1[0] * nro2[0] + nro1[1] * nro2[1] + nro1[2] * nro2[2];
  n = 1. + n1 + n2 + n3;

  if ((n == 0.) && (z == 0.)) {
    return (0.);
  }
  return (2. * atan2(z, n));
}

/* potential formulas  */
void doublelayer_lindholm_C(double x[3], double x1[3], double x2[3],
    double x3[3], double res[3]) {
  int i;
  double side1[3], side2[3], side3[3];
  double s1, s2, s3;
  double area;
  double normt[3];
  double eta1[3], eta2[3], eta3[3];
  double gamma[3][3];
  double ro1[3], ro2[3], ro3[3];
  double bro1, bro2, bro3;
  double zeta;
  double peta1, peta2, peta3;
  double omega;
  double P1, P2, P3;
  double g1, g2, g3;

  // s1, s2, s3
  for (i = 0; i < 3; i++) {
    side1[i] = x2[i] - x1[i];
    side2[i] = x3[i] - x2[i];
    side3[i] = x1[i] - x3[i];
  }
  s1 = norm3(side1);
  s2 = norm3(side2);
  s3 = norm3(side3);

  { // area
    double tmp[3];
    cross(side1, side3, tmp);
    area = 0.5 * norm3(tmp);
  }

  // normalize side vectors
  for (i = 0; i < 3; i++) {
    side1[i] = side1[i] / s1;
    side2[i] = side2[i] / s2;
    side3[i] = side3[i] / s3;
  }

  // normt
  cross(side3, side1, normt);
  normalize(normt);

  // eta1, eta2, eta3
  cross(normt, side1, eta1);
  cross(normt, side2, eta2);
  cross(normt, side3, eta3);
  normalize(eta1);
  normalize(eta2);
  normalize(eta3);

  // gamma
  gamma[0][0] = dot(side2, side1);
  gamma[0][1] = dot(side2, side2);
  gamma[0][2] = dot(side2, side3);
  gamma[1][0] = dot(side3, side1);
  gamma[1][1] = dot(side3, side2);
  gamma[1][2] = dot(side3, side3);
  gamma[2][0] = dot(side1, side1);
  gamma[2][1] = dot(side1, side2);
  gamma[2][2] = dot(side1, side3);

  // ro1
  for (i = 0; i < 3; i++) {
    ro1[i] = x1[i] - x[i];
    ro2[i] = x2[i] - x[i];
    ro3[i] = x3[i] - x[i];
  }

  // bro1
  bro1 = norm3(ro1);
  bro2 = norm3(ro2);
  bro3 = norm3(ro3);

  // zeta
  zeta = dot(normt, ro1);

  // peta
  peta1 = dot(eta1, ro1);
  peta2 = dot(eta2, ro2);
  peta3 = dot(eta3, ro3);

  // omega (solidangle rech))
  omega = solidangle_rech(bro1, bro2, bro3, ro1, ro2, ro3);

  // P
  P1 = (bro1 + bro2 - s1) / (bro1 + bro2 + s1);
  P2 = (bro2 + bro3 - s2) / (bro2 + bro3 + s2);
  P3 = (bro3 + bro1 - s3) / (bro3 + bro1 + s3);

  P1 = (P1 > 0.) ? -log(P1) : 0.;
  P2 = (P2 > 0.) ? -log(P2) : 0.;
  P3 = (P3 > 0.) ? -log(P3) : 0.;

  // g = L_tria_i
  g1 = (s2 / (8. * M_PI * area))
      * (peta2 * omega - zeta
          * (gamma[0][0] * P1 + gamma[0][1] * P2 + gamma[0][2] * P3));

  g2 = (s3 / (8. * M_PI * area))
      * (peta3 * omega - zeta
          * (gamma[1][0] * P1 + gamma[1][1] * P2 + gamma[1][2] * P3));
  g3 = (s1 / (8. * M_PI * area))
      * (peta1 * omega - zeta
          * (gamma[2][0] * P1 + gamma[2][1] * P2 + gamma[2][2] * P3));

  res[0] = g1;
  res[1] = g2;
  res[2] = g3;

}

void assemble_doublelayer_lindholm_C(int NB, double nodes[][3], int NEB,
    int triangles[][3], double M[]) {
  int nb, neb, i;
  double res[3];

  for (neb = 0; neb < NEB; neb++) {
    for (nb = 0; nb < NB; nb++) {
      doublelayer_lindholm_C(nodes[nb], nodes[triangles[neb][0]],
          nodes[triangles[neb][1]], nodes[triangles[neb][2]], res);
      for (i = 0; i < 3; i++) {
        M[nb * NB + triangles[neb][i]] += res[i];
      }
    }
  }

  return;
}

static void fill_dlp_l_collocation_near_laplacebem3d(const uint * ridx,
    const uint * cidx, pcbem3d bem, bool ntrans, pamatrix N) {
  pcsurface3d gr = bem->gr;
  real (*gr_x)[3] = (real (*)[3]) gr->x;
  uint (*gr_t)[3] = (uint (*)[3]) gr->t;
  plistnode *v2t = bem->v2t;
  uint triangles = gr->triangles;
  uint vertices = gr->vertices;
  field *aa = N->a;
  uint rows = ntrans ? N->cols : N->rows;
  uint cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  ptri_list tl_r;

#ifdef USE_OPENMP
#pragma omp parallel if(!omp_in_parallel() && (cols > 64*(1 << max_pardepth))) num_threads(1 << max_pardepth)
  {
#endif

  ptri_list tl1_r;
  pvert_list vl_r;
  plistnode v;
  uint i, ii, j, jj, l, rj, vv;
  uint *tri_i;
  real res[3];

#ifdef USE_OPENMP
#pragma omp single
  {
#endif
  clear_amatrix(N);

  tl_r = NULL;

  rj = 0;
  for (i = 0; i < rows; ++i) {
    ii = (ridx == NULL ? i : ridx[i]);
    for (v = v2t[ii], vv = v->data; v->next != NULL;
        v = v->next, vv = v->data) {

      tl1_r = tl_r;
      while (tl1_r && tl1_r->t != vv) {
        tl1_r = tl1_r->next;
      }

      if (tl1_r == NULL) {
        tl1_r = tl_r = new_tri_list(tl_r);
        tl_r->t = vv;
        rj++;
      }

      tl1_r->vl = new_vert_list(tl1_r->vl);
      tl1_r->vl->v = i;
    }
  }

#ifdef USE_OPENMP
}
#pragma omp for
#endif
  for (j = 0; j < cols; ++j) {
    jj = (cidx == NULL ? j : cidx[j]);
    assert(jj < vertices);
    for (tl1_r = tl_r; tl1_r != NULL; tl1_r = tl1_r->next) {
      ii = tl1_r->t;
      assert(ii < triangles);

      tri_i = gr_t[ii];

      doublelayer_lindholm_C(gr_x[jj], gr_x[tri_i[0]], gr_x[tri_i[1]],
          gr_x[tri_i[2]], res);

      vl_r = tl1_r->vl;
      while (vl_r) {
        i = vl_r->v;
        if (i < rows) {
          ii = (ridx == NULL ? i : ridx[i]);
          for (l = 0; l < 3; ++l) {
            if (ii == tri_i[l]) {
              if (ntrans) {
                aa[j + i * ld] += res[l];
              } else {
                aa[i + j * ld] += res[l];
              }
            }
          }
        }
        vl_r = vl_r->next;
      }
    }
  }

#ifdef USE_OPENMP
}
#endif

  del_tri_list(tl_r);

}

pbem3d new_dlp_collocation_laplace_bem3d(pcsurface3d gr, uint q_regular,
    uint q_singular, basisfunctionbem3d basis_neumann,
    basisfunctionbem3d basis_dirichlet) {

  pbem3d bem;

  bem = new_bem3d(gr);

  if (basis_neumann == BASIS_LINEAR_BEM3D || basis_dirichlet
      == BASIS_LINEAR_BEM3D) {
    setup_vertex_to_triangle_map_bem3d(bem);
  }

  bem->sq = build_singquad2d(gr, q_regular, q_singular);

  bem->basis_neumann = basis_neumann;
  bem->basis_dirichlet = basis_dirichlet;

  bem->N_neumann = gr->vertices;
  bem->N_dirichlet = gr->vertices;

  assert(
      basis_neumann == BASIS_LINEAR_BEM3D && basis_dirichlet
          == BASIS_LINEAR_BEM3D);

  bem->nearfield = fill_dlp_l_collocation_near_laplacebem3d;
  bem->nearfield_far = fill_dlp_l_collocation_near_laplacebem3d;

  return bem;
}

int main(int argc, char **argv) {
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
  pavector x, b;

  init_h2lib(&argc, &argv);
  sw = new_stopwatch();

  q_reg = 2;
  q_sing = 1;
  clf = 32;
  eta = 2.0;
  eps_aca = 1.0e-4;
  tm = new_releucl_truncmode();
  eps_recomp = 1.0e-4;

//   gr = build_interactive_surface3d();
  gr = read_gmsh_surface3d("model.msh");

  printf("Geometry has %d vertices, %d edges, %d triangles\n", gr->vertices,
      gr->edges, gr->triangles);
  printf("================================\n");
  vertices = gr->vertices;

  bem = new_dlp_collocation_laplace_bem3d(gr, q_reg, q_sing, BASIS_LINEAR_BEM3D,
      BASIS_LINEAR_BEM3D);

//  K = new_zero_amatrix(vertices, vertices);

//  printf("Fill dense matrix K:\n");
  printf("Dense matrix K:\n");
//  start_stopwatch(sw);
////  assemble_doublelayer_lindholm_C(vertices, gr->x, gr->triangles,
////      (int (*)[3]) gr->t, K->a);
//  bem->nearfield(NULL, NULL, bem, false, K);
//  printf("  %.2f s\n", stop_stopwatch(sw));
//  printf("  %.3f MB\n", getsize_amatrix(K) / 1024.0 / 1024.0);
  printf("  %.3lf MB\n", (double)(vertices)*vertices * 8. / 1024.0 / 1024.0);
  printf("================================\n");

  printf("Setup and fill H-matrix K:\n");
  start_stopwatch(sw);
  root = build_bem3d_cluster(bem, clf, BASIS_LINEAR_BEM3D);
  printf("  %d clusters created\n", root->desc);
  broot = build_strict_block(root, root, &eta, admissible_2_cluster);
  printf("  %d blocks created\n", broot->desc);

  Kh = build_from_block_hmatrix(broot, 0);

//  setup_hmatrix_aprx_paca_bem3d(bem, root, root, broot, eps_aca);
  setup_hmatrix_aprx_aca_bem3d(bem, root, root, broot, eps_aca);
  setup_hmatrix_recomp_bem3d(bem, true, eps_recomp, true, eps_recomp);

  assemble_bem3d_hmatrix(bem, broot, Kh);

  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n", getsize_hmatrix(Kh) / 1024.0 / 1024.0);

//  printf("rel error: %.5e\n",
//      norm2diff_amatrix_hmatrix(Kh, K) / norm2_amatrix(K));
  printf("================================\n");

  printf("Convert H-matrix Kh to H2-matrix kh2:\n");
  start_stopwatch(sw);

  Kh2 = compress_hmatrix_h2matrix(Kh, tm, eps_recomp);

  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n",
      (getsize_h2matrix(Kh2) + getsize_clusterbasis(Kh2->rb)
       + getsize_clusterbasis(Kh2->cb))
      / 1024.0 / 1024.0);

//  printf("rel error: %.5e\n",
//      norm2diff_amatrix_h2matrix(Kh2, K) / norm2_amatrix(K));


  // Test H2Matrix - Vector - Multiplication
  printf("TEST Matrix Vector Multiplication\n");
  x = new_avector(vertices);
  random_avector(x);
  b = new_avector(vertices);
  clear_avector(b);
  mvm_h2matrix_avector(1., false, Kh2, x, b);

//  del_amatrix(K);
  del_cluster(root);
  del_block(broot);
  del_hmatrix(Kh);
  del_h2matrix(Kh2);
  del_bem3d(bem);
  del_surface3d(gr);
  del_stopwatch(sw);

  uninit_h2lib();

  return 0;
}
