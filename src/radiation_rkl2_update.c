//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadiationRKL2Update_cpu(real mu,real nu,real c0,real mut,real gamt,real dt,Field* Ym1,Field* Ym2,Field* Y0,Field* Lm1,Field* L0,Field* Ynew) {

//<USER_DEFINED>
  INPUT(Ym1);
  INPUT(Ym2);
  INPUT(Y0);
  INPUT(Lm1);
  INPUT(L0);
  OUTPUT(Ynew);
//<\USER_DEFINED>

//<EXTERNAL>
  real* ym1=Ym1->field_cpu;
  real* ym2=Ym2->field_cpu;
  real* y0=Y0->field_cpu;
  real* lm1=Lm1->field_cpu;
  real* l0=L0->field_cpu;
  real* ynew=Ynew->field_cpu;
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
//<\INTERNAL>

//<MAIN_LOOP>
  i=j=k=0;

#ifdef Z
  for (k=NGHZ;k<size_z-NGHZ;k++) {
#endif
#ifdef Y
    for (j=NGHY;j<size_y-NGHY;j++) {
#endif
#ifdef X
      for (i=NGHX;i<size_x-NGHX;i++) {
#endif
//<#>

        ynew[l]=mu*ym1[l]
               +nu*ym2[l]
               +c0*y0[l]
               +dt*(mut*lm1[l]+gamt*l0[l]);

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