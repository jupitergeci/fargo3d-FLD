//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadiationFLDFields_cpu(real clight,Field* Erad,Field* Rho,Field* KappaR,Field* RadR,Field* RadLambda,Field* RadDiff) {

//<USER_DEFINED>
  INPUT(Erad);
  INPUT(Rho);
  INPUT(KappaR);
  OUTPUT(RadR);
  OUTPUT(RadLambda);
  OUTPUT(RadDiff);
//<\USER_DEFINED>

//<EXTERNAL>
  real* erad=Erad->field_cpu;
  real* rho=Rho->field_cpu;
  real* kappar=KappaR->field_cpu;
  real* radr=RadR->field_cpu;
  real* radlambda=RadLambda->field_cpu;
  real* raddiff=RadDiff->field_cpu;
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
  real gradx;
  real grady;
  real gradz;
  real gradm;
  real gradp;
  real gradmag;
  real eradloc;
  real rholoc;
  real kappaloc;
  real chi;
  real R;
  real lambda;
//<\INTERNAL>

//<CONSTANT>
// real InvDiffXmed(Nx+2*NGHX);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
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

        gradx=0.0;
        grady=0.0;
        gradz=0.0;

#ifdef X
        gradm=(erad[l]-erad[lxm])*Inv_zone_size_xmed(i,j,k);
        gradp=(erad[lxp]-erad[l])*Inv_zone_size_xmed(ixp,j,k);
        gradx=0.5*(gradm+gradp);
#endif

#ifdef Y
        gradm=(erad[l]-erad[lym])/(ymed(j)-ymed(j-1));
        gradp=(erad[lyp]-erad[l])/(ymed(j+1)-ymed(j));
        grady=0.5*(gradm+gradp);
#endif

#ifdef Z
        gradm=(erad[l]-erad[lzm])/(ymed(j)*(zmed(k)-zmed(k-1)));
        gradp=(erad[lzp]-erad[l])/(ymed(j)*(zmed(k+1)-zmed(k)));
        gradz=0.5*(gradm+gradp);
#endif

        gradmag=sqrt(gradx*gradx+grady*grady+gradz*gradz);

        eradloc=erad[l];
        rholoc=rho[l];
        kappaloc=kappar[l];

        if (eradloc < 1e-30)
          eradloc=1e-30;

        if (rholoc < 1e-30)
          rholoc=1e-30;

        if (kappaloc < 1e-30)
          kappaloc=1e-30;

        chi=kappaloc*rholoc;

        if (chi < 1e-30)
          chi=1e-30;

        R=gradmag/(chi*eradloc);

        if (R < 0.0)
          R=0.0;

        /*
         * Levermore-Pomraning flux limiter:
         *
         * lambda -> 1/3 for R -> 0
         * lambda -> 1/R for R -> infinity
         */
        if (R < 1e8)
          lambda=(2.0+R)/(6.0+3.0*R+R*R);
        else
          lambda=1.0/R;

        radr[l]=R;
        radlambda[l]=lambda;
        raddiff[l]=clight*lambda/chi;

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