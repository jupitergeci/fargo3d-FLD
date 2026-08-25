//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadiationDiffusionFLDStep_cpu(real dt,Field* Erad,Field* Diff,Field* EradNew) {

//<USER_DEFINED>
  INPUT(Erad);
  INPUT(Diff);
  OUTPUT(EradNew);
//<\USER_DEFINED>

//<EXTERNAL>
  real* erad=Erad->field_cpu;
  real* diff=Diff->field_cpu;
  real* eradnew=EradNew->field_cpu;
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
  real d0;
  real dm;
  real dp;
  real Dm;
  real Dp;
  real denom;
  real gradm;
  real gradp;
  real divdiff;
//<\INTERNAL>

//<CONSTANT>
// real Sxi(Nx+2*NGHX);
// real Sxj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syj(Ny+2*NGHY);
// real Syk(Nz+2*NGHZ);
// real Szj(Ny+2*NGHY);
// real Szk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvDiffXmed(Nx+2*NGHX);
//<\CONSTANT>

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

        divdiff=0.0;
        d0=diff[l];

#ifdef X
        dm=diff[lxm];
        dp=diff[lxp];

        Dm=0.0;
        if ((d0 > 0.0) && (dm > 0.0)) {
          denom=d0+dm;
          if (denom > 0.0)
            Dm=2.0*d0*dm/denom;
        }

        Dp=0.0;
        if ((d0 > 0.0) && (dp > 0.0)) {
          denom=d0+dp;
          if (denom > 0.0)
            Dp=2.0*d0*dp/denom;
        }

        gradm=(erad[l]-erad[lxm])*Inv_zone_size_xmed(i,j,k);
        gradp=(erad[lxp]-erad[l])*Inv_zone_size_xmed(ixp,j,k);

        divdiff+=(Dp*gradp-Dm*gradm)*SurfX(j,k);
#endif

#ifdef Y
        dm=diff[lym];
        dp=diff[lyp];

        Dm=0.0;
        if ((d0 > 0.0) && (dm > 0.0)) {
          denom=d0+dm;
          if (denom > 0.0)
            Dm=2.0*d0*dm/denom;
        }

        Dp=0.0;
        if ((d0 > 0.0) && (dp > 0.0)) {
          denom=d0+dp;
          if (denom > 0.0)
            Dp=2.0*d0*dp/denom;
        }

        gradm=(erad[l]-erad[lym])/(ymed(j)-ymed(j-1));
        gradp=(erad[lyp]-erad[l])/(ymed(j+1)-ymed(j));

        divdiff+=Dp*gradp*SurfY(i,j+1,k)
                -Dm*gradm*SurfY(i,j,k);
#endif

#ifdef Z
        dm=diff[lzm];
        dp=diff[lzp];

        Dm=0.0;
        if ((d0 > 0.0) && (dm > 0.0)) {
          denom=d0+dm;
          if (denom > 0.0)
            Dm=2.0*d0*dm/denom;
        }

        Dp=0.0;
        if ((d0 > 0.0) && (dp > 0.0)) {
          denom=d0+dp;
          if (denom > 0.0)
            Dp=2.0*d0*dp/denom;
        }

        gradm=(erad[l]-erad[lzm])/
              (ymed(j)*(zmed(k)-zmed(k-1)));

        gradp=(erad[lzp]-erad[l])/
              (ymed(j)*(zmed(k+1)-zmed(k)));

        divdiff+=Dp*gradp*SurfZ(i,j,k+1)
                -Dm*gradm*SurfZ(i,j,k);
#endif

        eradnew[l]=erad[l]+dt*divdiff*InvVol(i,j,k);

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