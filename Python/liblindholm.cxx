// Copyright (C) 2016 Thomas Schrefl, Dieter Suess, Florian Bruckner
// Last modified by Florian Bruckner, 2016-06-17

#include <math.h>
#include <sstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "liblindholm.h"
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
static void doublelayer_lindholm_C(double x[3], double x1[3], double x2[3],
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
                                     int triangles[][3], double M[])
{
  int       nb, neb, i;
  double    res[3];

  for (neb = 0; neb < NEB; neb++) {
    for (nb = 0; nb < NB; nb++) {
      doublelayer_lindholm_C(nodes[nb], nodes[triangles[neb][0]],
			     nodes[triangles[neb][1]],
			     nodes[triangles[neb][2]], res);
      for (i = 0; i < 3; i++) {
	M[nb * NB + triangles[neb][i]] += res[i];
      }
    }
  }

  return;
}

/*** H2Lib interface ***/
static void
fill_dlp_l_collocation_near_laplacebem3d(const uint * ridx,
					 const uint * cidx, pcbem3d bem,
					 bool ntrans, pamatrix N)
{
  pcsurface3d gr = bem->gr;
  real(*gr_x)[3] = (real(*)[3]) gr->x;
  uint(*gr_t)[3] = (uint(*)[3]) gr->t;
  plistnode *v2t = bem->v2t;
  uint      triangles = gr->triangles;
  uint      vertices = gr->vertices;
  field    *aa = N->a;
  uint      rows = ntrans ? N->cols : N->rows;
  uint      cols = ntrans ? N->rows : N->cols;
  longindex ld = N->ld;

  ptri_list tl_c;

#ifdef USE_OPENMP
#pragma omp parallel if(!omp_in_parallel() && (cols > 64*(1 << max_pardepth))) num_threads(1 << max_pardepth)
  {
#endif

    ptri_list tl1_c;
    pvert_list vl_c;
    plistnode v;
    uint      j, jj, i, ii, l, cj, vv;
    uint     *tri_j;
    real      res[3];

#ifdef USE_OPENMP
#pragma omp single
    {
#endif
      clear_amatrix(N);

      tl_c = NULL;

      cj = 0;
      for (j = 0; j < cols; ++j) {
	jj = (cidx == NULL ? j : cidx[j]);
	for (v = v2t[jj], vv = v->data; v->next != NULL;
	     v = v->next, vv = v->data) {

	  tl1_c = tl_c;
	  while (tl1_c && tl1_c->t != vv) {
	    tl1_c = tl1_c->next;
	  }

	  if (tl1_c == NULL) {
	    tl1_c = tl_c = new_tri_list(tl_c);
	    tl_c->t = vv;
	    cj++;
	  }

	  tl1_c->vl = new_vert_list(tl1_c->vl);
	  tl1_c->vl->v = j;
	}
      }

#ifdef USE_OPENMP
    }
#pragma omp for
#endif
    for (i = 0; i < rows; ++i) {
      ii = (ridx == NULL ? i : ridx[i]);
      assert(ii < vertices);
      for (tl1_c = tl_c; tl1_c != NULL; tl1_c = tl1_c->next) {
	jj = tl1_c->t;
	assert(jj < triangles);

	tri_j = gr_t[jj];

	doublelayer_lindholm_C(gr_x[ii], gr_x[tri_j[0]], gr_x[tri_j[1]],
			       gr_x[tri_j[2]], res);

	vl_c = tl1_c->vl;
	while (vl_c) {
	  j = vl_c->v;
	  if (j < cols) {
	    jj = (cidx == NULL ? j : cidx[j]);
	    for (l = 0; l < 3; ++l) {
	      if (jj == tri_j[l]) {
		if (ntrans) {
		  aa[j + i * ld] += res[l];
		}
		else {
		  aa[i + j * ld] += res[l];
		}
	      }
	    }
	  }
	  vl_c = vl_c->next;
	}
      }
    }

#ifdef USE_OPENMP
  }
#endif

  del_tri_list(tl_c);

}

static void assemble_fundamental_collocation_row_bem3d(const uint * idx, const real(*Z)[3], pcbem3d bem, bool trans, pamatrix A)
{
  real(*X)[3] = (real(*)[3]) allocreal(3 * A->rows);
  uint      i, ii;

  for (i = 0; i < A->rows; ++i) {
    ii = (idx != NULL ? idx[i] : i);
    X[i][0] = bem->gr->x[ii][0];
    X[i][1] = bem->gr->x[ii][1];
    X[i][2] = bem->gr->x[ii][2];
  }

  bem->kernels->fundamental(bem, X, Z, A);

  freemem(X);
}

static void assemble_dnz_fundamental_collocation_row_bem3d(const uint * idx, const real(*Z)[3], const real(*N)[3], pcbem3d bem, bool trans, pamatrix A)
{
  real(*X)[3] = (real(*)[3]) allocreal(3 * A->rows);
  uint      i, ii;

  for (i = 0; i < A->rows; ++i) {
    ii = (idx != NULL ? idx[i] : i);
    X[i][0] = bem->gr->x[ii][0];
    X[i][1] = bem->gr->x[ii][1];
    X[i][2] = bem->gr->x[ii][2];
  }

  bem->kernels->dny_fundamental(bem, X, Z, N, A);

  freemem(X);
}

static void assemble_lagrange_collocation_row_bem3d(const uint * idx, pcrealavector px, pcrealavector py, pcrealavector pz, pcbem3d bem, pamatrix A)
{
  real(*X)[3] = (real(*)[3]) allocreal(3 * A->rows);
  uint      i, ii;

  for (i = 0; i < A->rows; ++i) {
    ii = (idx != NULL ? idx[i] : i);
    X[i][0] = bem->gr->x[ii][0];
    X[i][1] = bem->gr->x[ii][1];
    X[i][2] = bem->gr->x[ii][2];
  }

  assemble_bem3d_lagrange_amatrix(X, px, py, pz, bem, A);

  freemem(X);
}

static pbem3d new_dlp_collocation_laplace_bem3d(pcsurface3d gr, uint q_regular, uint q_singular, basisfunctionbem3d basis_neumann, basisfunctionbem3d basis_dirichlet)
{

  pbem3d    bem;

  bem = new_dlp_laplace_bem3d(gr, q_regular, q_singular, basis_neumann,
			      basis_dirichlet, 0.5);

  assert(basis_neumann == BASIS_LINEAR_BEM3D && basis_dirichlet
	 == BASIS_LINEAR_BEM3D);

  bem->nearfield = fill_dlp_l_collocation_near_laplacebem3d;
  bem->nearfield_far = fill_dlp_l_collocation_near_laplacebem3d;

  bem->kernels->fundamental_row = assemble_fundamental_collocation_row_bem3d;
  bem->kernels->dnz_fundamental_row = assemble_dnz_fundamental_collocation_row_bem3d;
  bem->kernels->lagrange_row = assemble_lagrange_collocation_row_bem3d;

  return bem;
}

int Lindholm_C::geometry_from_file(std::string infile)
{
  gr = read_gmsh_surface3d(infile.c_str());
  printf("Geometry has %d vertices, %d edges, %d triangles\n", gr->vertices,
      gr->edges, gr->triangles);
  printf("================================\n");
  vertices = gr->vertices;
  return(0);
}

int Lindholm_C::geometry_from_array(unsigned int N, double coordinates[][3], unsigned int NE, int elements[][3])
{
  gr = create_surface3d_fromarray(N, coordinates, NE, elements);

  printf("Geometry has %d vertices, %d edges, %d triangles\n", gr->vertices,
      gr->edges, gr->triangles);
  printf("================================\n");
  vertices = gr->vertices;
  return(0);
}

int Lindholm_C::_setup_init()
{
  bem = new_dlp_collocation_laplace_bem3d(gr, q_reg, q_sing, BASIS_LINEAR_BEM3D, BASIS_LINEAR_BEM3D);

  printf("Dense matrix K:\n");
  printf("  %.3lf MB\n", (double)(vertices)*vertices * 8. / 1024.0 / 1024.0);
  printf("================================\n");

  root = build_bem3d_cluster(bem, clf, BASIS_LINEAR_BEM3D);
  printf("  %d clusters created\n", root->desc);
  broot = build_strict_block(root, root, &eta, admissible_2_cluster);
  printf("  %d blocks created\n", broot->desc);
  return(0);
}

int Lindholm_C::_setup_cleanup()
{ // TODO: check what can be free'd
  if (broot != NULL) del_block(broot);
  if (bem != NULL) del_bem3d(bem);
  return(0);
}

int Lindholm_C::setup_HCA()
{
  pstopwatch sw = new_stopwatch();
  _setup_init();
  printf("Setup and fill H-matrix K:\n");
  start_stopwatch(sw);
  phmatrix Kh;
  Kh = build_from_block_hmatrix(broot, 0);

//  setup_hmatrix_aprx_paca_bem3d(bem, root, root, broot, eps_aca);
//  setup_hmatrix_aprx_aca_bem3d(bem, root, root, broot, eps_aca);
  setup_hmatrix_aprx_hca_bem3d(bem, root, root, broot, 5, eps_aca);
//  setup_hmatrix_aprx_inter_row_bem3d(bem, root, root, broot, 5);
//  setup_hmatrix_recomp_bem3d(bem, true, eps_recomp, true, eps_recomp); // TODO: enable this?
  assemble_bem3d_hmatrix(bem, broot, Kh);

  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n", getsize_hmatrix(Kh) / 1024.0 / 1024.0);
  printf("================================\n");
  printf("Convert H-matrix Kh to H2-matrix kh2:\n");
  start_stopwatch(sw);
  Kh2 = compress_hmatrix_h2matrix(Kh, tm, eps_recomp);
  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n",
      (getsize_h2matrix(Kh2) + getsize_clusterbasis(Kh2->rb)
       + getsize_clusterbasis(Kh2->cb))
      / 1024.0 / 1024.0);
  del_hmatrix(Kh); // TODO: free H-matrix
  del_stopwatch(sw);
  _setup_cleanup();
  return(0);
}

int Lindholm_C::setup_GCA()
{
  pstopwatch sw = new_stopwatch();
  pclusterbasis rb, cb;
  _setup_init();
  rb = build_from_cluster_clusterbasis(root); // TODO: need to be freed?
  cb = build_from_cluster_clusterbasis(root);

  setup_h2matrix_aprx_greenhybrid_bem3d(bem, rb, cb, broot, m, 1, delta,
					eps_aca, build_bem3d_cube_quadpoints);

  printf("assemble row basis:\n");
  start_stopwatch(sw);
  assemble_bem3d_h2matrix_row_clusterbasis(bem, rb);
  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n", getsize_clusterbasis(rb) / 1024.0 / 1024.0);
  printf("================================\n");

  printf("assemble col basis:\n");
  start_stopwatch(sw);
  assemble_bem3d_h2matrix_col_clusterbasis(bem, cb);
  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n", getsize_clusterbasis(cb) / 1024.0 / 1024.0);
  printf("================================\n");

  Kh2 = build_from_block_h2matrix(broot, rb, cb);

  printf("Assemble H2-matrix:\n");
  start_stopwatch(sw);
  assemble_bem3d_h2matrix(bem, Kh2);
  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n", getsize_h2matrix(Kh2) / 1024.0 / 1024.0);
  printf("================================\n");

  printf("Recompress H2-matrix:\n");
  start_stopwatch(sw);
  recompress_inplace_h2matrix(Kh2, tm, eps_recomp);
  printf("  %.2f s\n", stop_stopwatch(sw));
  printf("  %.3f MB\n", getsize_h2matrix(Kh2) / 1024.0 / 1024.0);
  printf("================================\n");
  del_stopwatch(sw);
  _setup_cleanup();
  return(0);
}

int Lindholm_C::matvec(unsigned int N, double x[], double b[])
{
  printf("Matrix Vector Multiplication\n");
  pavector x_vec, b_vec;
  x_vec = new_pointer_avector((double *)x, N);
  b_vec = new_pointer_avector((double *)b, N);
  mvm_h2matrix_avector(1., false, Kh2, x_vec, b_vec);
  del_avector(x_vec);
  del_avector(b_vec);
  return(0);
}

int Lindholm_C::get_size()
{
  int size = -1;

  if ((Kh2 != NULL) && (Kh2->rb != NULL) && (Kh2->cb != 0)) {
    size = getsize_h2matrix(Kh2) + getsize_clusterbasis(Kh2->rb) + getsize_clusterbasis(Kh2->cb);
  }
  return(size);
}


Lindholm_C::Lindholm_C()
{
  q_reg = 3;
  q_sing = q_reg + 2;
  clf = 32;
  eta = 2.0;
  m = 2;
  delta = 2.0;
  eps_aca = 1.0e-8;
  tm = new_releucl_truncmode();
  eps_recomp = 2.0e-3;

  root = NULL;
  broot = NULL;
  Kh2 = NULL;
  bem = NULL;
  gr = NULL;
}

Lindholm_C::~Lindholm_C()
{
  if (root != NULL) del_cluster(root);
  if (Kh2 != NULL) del_h2matrix(Kh2);
  if (gr != NULL) del_surface3d(gr);
}

int main( void )
{
  Lindholm_C h2lib;
  h2lib.geometry_from_file("model.msh");
//  h2lib.setup_HCA();
  h2lib.setup_GCA();

//  printf("TEST Matrix Vector Multiplication\n");
//  pavector x, b;
//  x = new_avector(vertices);
//  random_avector(x);
//  b = new_avector(vertices);
//  clear_avector(b);
//  mvm_h2matrix_avector(1., false, h2lib.Kh2, x, b);
  printf("SIZE: %d\n", h2lib.get_size());
}

// TODO: init / uninit?
