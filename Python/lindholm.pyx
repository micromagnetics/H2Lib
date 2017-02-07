import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libcpp.string cimport string
from warnings import warn

### C++ definitions
cdef extern from "liblindholm.h":
  cdef cppclass Lindholm_C:
    Lindholm_C()
    int geometry_from_file(string infile) except +
    int geometry_from_array(unsigned int N, double coordinates[][3], unsigned int NE, int cells[][3]) except +
    int setup_HCA() except +
    int setup_GCA() except +
    int matvec(unsigned int N, double x[], double b[]) except +
    int get_size()

### Python wrappers
cdef class Lindholm:
  cdef Lindholm_C *cobj

  def __cinit__(self):
    self.cobj = new Lindholm_C()

  def geometry_from_file(self, string infile):
    return self.cobj.geometry_from_file(infile)

  def geometry_from_array(self, np.ndarray[np.float64_t,ndim=2] nodes, np.ndarray[np.uint32_t,ndim=2] cells):
    cdef int N = nodes.shape[0]
    cdef int NE = cells.shape[0]
    nodes = np.ascontiguousarray(nodes)
    cells = np.ascontiguousarray(cells)
    return self.cobj.geometry_from_array(N, <double(*)[3]> nodes.data, NE, <int(*)[3]> cells.data)

  def setup(self, method="GCA"):
    if method == "GCA":
      return self.cobj.setup_GCA()
    elif method == "HCA":
      return self.cobj.setup_HCA()
    else:
      raise RuntimeError("H2Matrix assembly method '%s' not known (use 'GCA' or 'HCA')" % method)

  def matvec(self, np.ndarray[np.float64_t,ndim=1] u1):
    cdef int N = u1.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] u2 = np.zeros(N)
    if not u1.flags['C_CONTIGUOUS']:
      warn("ndarray mag should be continguous in order to increase performance!")
      u1 = np.ascontiguousarray(u1)

    self.cobj.matvec(N, <double(*)> u1.data, <double(*)> u2.data)
    return u2

  def get_size(self):
    return self.cobj.get_size()

  def __dealloc__(self):
    if self.cobj != NULL:
      del self.cobj
