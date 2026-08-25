//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadiationSetVariableTest_cpu(real amp,real d0,Field* Erad,Field* Diff) {

//<USER_DEFINED>
  OUTPUT(Erad);
  OUTPUT(Diff);
//<\USER_DEFINED>

//<EXTERNAL>
  real* erad=Erad->field_cpu;
  real* diff=Diff->field_cpu;
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

        erad[l]=amp*r*r;
        diff[l]=d0*r;

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