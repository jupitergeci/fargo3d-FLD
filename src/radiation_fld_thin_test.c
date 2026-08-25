//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadiationSetThinTest_cpu(real rho0,real slope,real r0,Field* Erad,Field* Rho) {

//<USER_DEFINED>
  OUTPUT(Erad);
  OUTPUT(Rho);
//<\USER_DEFINED>

//<EXTERNAL>
  real* erad=Erad->field_cpu;
  real* rho=Rho->field_cpu;
  int pitch=Pitch_cpu;
  int stride=Stride_cpu;
  int size_x=Nx+2*NGHX;
  int size_y=Ny+2*NGHY;
  int size_z=Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  real r;
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
//<\CONSTANT>

//<MAIN_LOOP>
  i=j=k=0;

#ifdef Z
  for (k=0;k<size_z;k++) {
#endif
#ifdef Y
    for (j=0;j<size_y;j++) {
#endif
#ifdef X
      for (i=0;i<size_x;i++) {
#endif
//<#>

        r=ymed(j);

        rho[l]=rho0;

        /*
         * Monotonic exponential radiation profile.
         *
         * |grad(E_R)|/E_R = slope
         *
         * for a purely radial field, giving approximately
         *
         * R = slope/(rho kappa_R)
         *
         * throughout the domain.
         */
        erad[l]=exp(-slope*(r-r0));

//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}