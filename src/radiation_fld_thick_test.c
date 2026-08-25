//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadiationSetThickTest_cpu(real rho0,real amp,real sigma,real rc,Field* Erad,Field* Rho) {

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
  real x;
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
        x=(r-rc)/sigma;

        rho[l]=rho0;

        /*
         * Smooth, low-amplitude radiation perturbation.
         *
         * The constant background keeps E_R positive while the
         * large optical depth forces R << 1.
         */
        erad[l]=1.0+amp*exp(-0.5*x*x);

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