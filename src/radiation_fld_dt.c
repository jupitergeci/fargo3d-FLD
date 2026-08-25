//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadiationFLDDtField_cpu(real kappa,real clight,Field* Erad,Field* Rho,Field* DtField) {

//<USER_DEFINED>
  INPUT(Erad);
  INPUT(Rho);
  OUTPUT(DtField);
//<\USER_DEFINED>

//<EXTERNAL>
  real* erad=Erad->field_cpu;
  real* rho=Rho->field_cpu;
  real* dtfield=DtField->field_cpu;
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

  real EL;
  real ER;
  real rhoL;
  real rhoR;
  real Eface;
  real rhoface;
  real chi;

  real gn;
  real gt1;
  real gt2;
  real gLm;
  real gLp;
  real gRm;
  real gRp;
  real gL;
  real gR;
  real gradmag;

  real R;
  real lambda;
  real Dface;

  real dist;
  real den;
  real rate;

  real wL;
  real wR;
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

        rate=0.0;


/* ================================================================
 * X- FACE
 * ================================================================ */
#ifdef X

        EL=erad[lxm];
        ER=erad[l];

        rhoL=rho[lxm];
        rhoR=rho[l];

        if (EL < 1e-30) EL=1e-30;
        if (ER < 1e-30) ER=1e-30;
        if (rhoL < 1e-30) rhoL=1e-30;
        if (rhoR < 1e-30) rhoR=1e-30;

        Eface=0.5*(EL+ER);
        rhoface=0.5*(rhoL+rhoR);

        if (Eface < 1e-30) Eface=1e-30;
        if (rhoface < 1e-30) rhoface=1e-30;

        chi=kappa*rhoface;

        if (chi < 1e-30)
          chi=1e-30;

        gn=(ER-EL)*Inv_zone_size_xmed(i,j,k);

        gt1=0.0;
        gt2=0.0;

#ifdef Y
        gLm=(erad[lxm]-erad[lxm-pitch])/
            (ymed(j)-ymed(j-1));

        gLp=(erad[lxm+pitch]-erad[lxm])/
            (ymed(j+1)-ymed(j));

        gL=0.5*(gLm+gLp);

        gRm=(erad[l]-erad[lym])/
            (ymed(j)-ymed(j-1));

        gRp=(erad[lyp]-erad[l])/
            (ymed(j+1)-ymed(j));

        gR=0.5*(gRm+gRp);

        gt1=0.5*(gL+gR);
#endif

#ifdef Z
        gLm=(erad[lxm]-erad[lxm-stride])/
            (ymed(j)*(zmed(k)-zmed(k-1)));

        gLp=(erad[lxm+stride]-erad[lxm])/
            (ymed(j)*(zmed(k+1)-zmed(k)));

        gL=0.5*(gLm+gLp);

        gRm=(erad[l]-erad[lzm])/
            (ymed(j)*(zmed(k)-zmed(k-1)));

        gRp=(erad[lzp]-erad[l])/
            (ymed(j)*(zmed(k+1)-zmed(k)));

        gR=0.5*(gRm+gRp);

        gt2=0.5*(gL+gR);
#endif

        gradmag=sqrt(gn*gn+gt1*gt1+gt2*gt2);

        R=gradmag/(chi*Eface);

        if (R < 0.0)
          R=0.0;

        if (R < 1e8)
          lambda=(2.0+R)/(6.0+3.0*R+R*R);
        else
          lambda=1.0/R;

        Dface=clight*lambda/chi;

        dist=1.0/Inv_zone_size_xmed(i,j,k);

        if (dist > 0.0)
          rate+=Dface*SurfX(j,k)/dist;


/* ================================================================
 * X+ FACE
 * ================================================================ */

        EL=erad[l];
        ER=erad[lxp];

        rhoL=rho[l];
        rhoR=rho[lxp];

        if (EL < 1e-30) EL=1e-30;
        if (ER < 1e-30) ER=1e-30;
        if (rhoL < 1e-30) rhoL=1e-30;
        if (rhoR < 1e-30) rhoR=1e-30;

        Eface=0.5*(EL+ER);
        rhoface=0.5*(rhoL+rhoR);

        if (Eface < 1e-30) Eface=1e-30;
        if (rhoface < 1e-30) rhoface=1e-30;

        chi=kappa*rhoface;

        if (chi < 1e-30)
          chi=1e-30;

        gn=(ER-EL)*Inv_zone_size_xmed(ixp,j,k);

        gt1=0.0;
        gt2=0.0;

#ifdef Y
        gLm=(erad[l]-erad[lym])/
            (ymed(j)-ymed(j-1));

        gLp=(erad[lyp]-erad[l])/
            (ymed(j+1)-ymed(j));

        gL=0.5*(gLm+gLp);

        gRm=(erad[lxp]-erad[lxp-pitch])/
            (ymed(j)-ymed(j-1));

        gRp=(erad[lxp+pitch]-erad[lxp])/
            (ymed(j+1)-ymed(j));

        gR=0.5*(gRm+gRp);

        gt1=0.5*(gL+gR);
#endif

#ifdef Z
        gLm=(erad[l]-erad[lzm])/
            (ymed(j)*(zmed(k)-zmed(k-1)));

        gLp=(erad[lzp]-erad[l])/
            (ymed(j)*(zmed(k+1)-zmed(k)));

        gL=0.5*(gLm+gLp);

        gRm=(erad[lxp]-erad[lxp-stride])/
            (ymed(j)*(zmed(k)-zmed(k-1)));

        gRp=(erad[lxp+stride]-erad[lxp])/
            (ymed(j)*(zmed(k+1)-zmed(k)));

        gR=0.5*(gRm+gRp);

        gt2=0.5*(gL+gR);
#endif

        gradmag=sqrt(gn*gn+gt1*gt1+gt2*gt2);

        R=gradmag/(chi*Eface);

        if (R < 0.0)
          R=0.0;

        if (R < 1e8)
          lambda=(2.0+R)/(6.0+3.0*R+R*R);
        else
          lambda=1.0/R;

        Dface=clight*lambda/chi;

        dist=1.0/Inv_zone_size_xmed(ixp,j,k);

        if (dist > 0.0)
          rate+=Dface*SurfX(j,k)/dist;

#endif


/* ================================================================
 * Y- FACE
 * ================================================================ */
#ifdef Y

        EL=erad[lym];
        ER=erad[l];

        rhoL=rho[lym];
        rhoR=rho[l];

        if (EL < 1e-30) EL=1e-30;
        if (ER < 1e-30) ER=1e-30;
        if (rhoL < 1e-30) rhoL=1e-30;
        if (rhoR < 1e-30) rhoR=1e-30;

        den=ymed(j)-ymed(j-1);

        wL=(ymed(j)-ymin(j))/den;
        wR=(ymin(j)-ymed(j-1))/den;

        Eface=wL*EL+wR*ER;
        rhoface=wL*rhoL+wR*rhoR;

        if (Eface < 1e-30) Eface=1e-30;
        if (rhoface < 1e-30) rhoface=1e-30;

        chi=kappa*rhoface;

        if (chi < 1e-30)
          chi=1e-30;

        gn=(ER-EL)/den;

        gt1=0.0;
        gt2=0.0;

#ifdef X
        gLm=(erad[lym]-erad[lxm-pitch])*
            Inv_zone_size_xmed(i,j-1,k);

        gLp=(erad[lxp-pitch]-erad[lym])*
            Inv_zone_size_xmed(ixp,j-1,k);

        gL=0.5*(gLm+gLp);

        gRm=(erad[l]-erad[lxm])*
            Inv_zone_size_xmed(i,j,k);

        gRp=(erad[lxp]-erad[l])*
            Inv_zone_size_xmed(ixp,j,k);

        gR=0.5*(gRm+gRp);

        gt1=wL*gL+wR*gR;
#endif

#ifdef Z
        gLm=(erad[lym]-erad[lym-stride])/
            (ymed(j-1)*(zmed(k)-zmed(k-1)));

        gLp=(erad[lym+stride]-erad[lym])/
            (ymed(j-1)*(zmed(k+1)-zmed(k)));

        gL=0.5*(gLm+gLp);

        gRm=(erad[l]-erad[lzm])/
            (ymed(j)*(zmed(k)-zmed(k-1)));

        gRp=(erad[lzp]-erad[l])/
            (ymed(j)*(zmed(k+1)-zmed(k)));

        gR=0.5*(gRm+gRp);

        gt2=wL*gL+wR*gR;
#endif

        gradmag=sqrt(gn*gn+gt1*gt1+gt2*gt2);

        R=gradmag/(chi*Eface);

        if (R < 0.0)
          R=0.0;

        if (R < 1e8)
          lambda=(2.0+R)/(6.0+3.0*R+R*R);
        else
          lambda=1.0/R;

        Dface=clight*lambda/chi;

        if (den > 0.0)
          rate+=Dface*SurfY(i,j,k)/den;


/* ================================================================
 * Y+ FACE
 * ================================================================ */

        EL=erad[l];
        ER=erad[lyp];

        rhoL=rho[l];
        rhoR=rho[lyp];

        if (EL < 1e-30) EL=1e-30;
        if (ER < 1e-30) ER=1e-30;
        if (rhoL < 1e-30) rhoL=1e-30;
        if (rhoR < 1e-30) rhoR=1e-30;

        den=ymed(j+1)-ymed(j);

        wL=(ymed(j+1)-ymin(j+1))/den;
        wR=(ymin(j+1)-ymed(j))/den;

        Eface=wL*EL+wR*ER;
        rhoface=wL*rhoL+wR*rhoR;

        if (Eface < 1e-30) Eface=1e-30;
        if (rhoface < 1e-30) rhoface=1e-30;

        chi=kappa*rhoface;

        if (chi < 1e-30)
          chi=1e-30;

        gn=(ER-EL)/den;

        gt1=0.0;
        gt2=0.0;

#ifdef X
        gLm=(erad[l]-erad[lxm])*
            Inv_zone_size_xmed(i,j,k);

        gLp=(erad[lxp]-erad[l])*
            Inv_zone_size_xmed(ixp,j,k);

        gL=0.5*(gLm+gLp);

        gRm=(erad[lyp]-erad[lxm+pitch])*
            Inv_zone_size_xmed(i,j+1,k);

        gRp=(erad[lxp+pitch]-erad[lyp])*
            Inv_zone_size_xmed(ixp,j+1,k);

        gR=0.5*(gRm+gRp);

        gt1=wL*gL+wR*gR;
#endif

#ifdef Z
        gLm=(erad[l]-erad[lzm])/
            (ymed(j)*(zmed(k)-zmed(k-1)));

        gLp=(erad[lzp]-erad[l])/
            (ymed(j)*(zmed(k+1)-zmed(k)));

        gL=0.5*(gLm+gLp);

        gRm=(erad[lyp]-erad[lyp-stride])/
            (ymed(j+1)*(zmed(k)-zmed(k-1)));

        gRp=(erad[lyp+stride]-erad[lyp])/
            (ymed(j+1)*(zmed(k+1)-zmed(k)));

        gR=0.5*(gRm+gRp);

        gt2=wL*gL+wR*gR;
#endif

        gradmag=sqrt(gn*gn+gt1*gt1+gt2*gt2);

        R=gradmag/(chi*Eface);

        if (R < 0.0)
          R=0.0;

        if (R < 1e8)
          lambda=(2.0+R)/(6.0+3.0*R+R*R);
        else
          lambda=1.0/R;

        Dface=clight*lambda/chi;

        if (den > 0.0)
          rate+=Dface*SurfY(i,j+1,k)/den;

#endif


/* ================================================================
 * Z- FACE
 * ================================================================ */
#ifdef Z

        EL=erad[lzm];
        ER=erad[l];

        rhoL=rho[lzm];
        rhoR=rho[l];

        if (EL < 1e-30) EL=1e-30;
        if (ER < 1e-30) ER=1e-30;
        if (rhoL < 1e-30) rhoL=1e-30;
        if (rhoR < 1e-30) rhoR=1e-30;

        den=zmed(k)-zmed(k-1);

        wL=(zmed(k)-zmin(k))/den;
        wR=(zmin(k)-zmed(k-1))/den;

        Eface=wL*EL+wR*ER;
        rhoface=wL*rhoL+wR*rhoR;

        if (Eface < 1e-30) Eface=1e-30;
        if (rhoface < 1e-30) rhoface=1e-30;

        chi=kappa*rhoface;

        if (chi < 1e-30)
          chi=1e-30;

        gn=(ER-EL)/(ymed(j)*den);

        gt1=0.0;
        gt2=0.0;

#ifdef X
        gLm=(erad[lzm]-erad[lxm-stride])*
            Inv_zone_size_xmed(i,j,k-1);

        gLp=(erad[lxp-stride]-erad[lzm])*
            Inv_zone_size_xmed(ixp,j,k-1);

        gL=0.5*(gLm+gLp);

        gRm=(erad[l]-erad[lxm])*
            Inv_zone_size_xmed(i,j,k);

        gRp=(erad[lxp]-erad[l])*
            Inv_zone_size_xmed(ixp,j,k);

        gR=0.5*(gRm+gRp);

        gt1=wL*gL+wR*gR;
#endif

#ifdef Y
        gLm=(erad[lzm]-erad[lzm-pitch])/
            (ymed(j)-ymed(j-1));

        gLp=(erad[lzm+pitch]-erad[lzm])/
            (ymed(j+1)-ymed(j));

        gL=0.5*(gLm+gLp);

        gRm=(erad[l]-erad[lym])/
            (ymed(j)-ymed(j-1));

        gRp=(erad[lyp]-erad[l])/
            (ymed(j+1)-ymed(j));

        gR=0.5*(gRm+gRp);

        gt2=wL*gL+wR*gR;
#endif

        gradmag=sqrt(gn*gn+gt1*gt1+gt2*gt2);

        R=gradmag/(chi*Eface);

        if (R < 0.0)
          R=0.0;

        if (R < 1e8)
          lambda=(2.0+R)/(6.0+3.0*R+R*R);
        else
          lambda=1.0/R;

        Dface=clight*lambda/chi;

        dist=ymed(j)*den;

        if (dist > 0.0)
          rate+=Dface*SurfZ(i,j,k)/dist;


/* ================================================================
 * Z+ FACE
 * ================================================================ */

        EL=erad[l];
        ER=erad[lzp];

        rhoL=rho[l];
        rhoR=rho[lzp];

        if (EL < 1e-30) EL=1e-30;
        if (ER < 1e-30) ER=1e-30;
        if (rhoL < 1e-30) rhoL=1e-30;
        if (rhoR < 1e-30) rhoR=1e-30;

        den=zmed(k+1)-zmed(k);

        wL=(zmed(k+1)-zmin(k+1))/den;
        wR=(zmin(k+1)-zmed(k))/den;

        Eface=wL*EL+wR*ER;
        rhoface=wL*rhoL+wR*rhoR;

        if (Eface < 1e-30) Eface=1e-30;
        if (rhoface < 1e-30) rhoface=1e-30;

        chi=kappa*rhoface;

        if (chi < 1e-30)
          chi=1e-30;

        gn=(ER-EL)/(ymed(j)*den);

        gt1=0.0;
        gt2=0.0;

#ifdef X
        gLm=(erad[l]-erad[lxm])*
            Inv_zone_size_xmed(i,j,k);

        gLp=(erad[lxp]-erad[l])*
            Inv_zone_size_xmed(ixp,j,k);

        gL=0.5*(gLm+gLp);

        gRm=(erad[lzp]-erad[lxm+stride])*
            Inv_zone_size_xmed(i,j,k+1);

        gRp=(erad[lxp+stride]-erad[lzp])*
            Inv_zone_size_xmed(ixp,j,k+1);

        gR=0.5*(gRm+gRp);

        gt1=wL*gL+wR*gR;
#endif

#ifdef Y
        gLm=(erad[l]-erad[lym])/
            (ymed(j)-ymed(j-1));

        gLp=(erad[lyp]-erad[l])/
            (ymed(j+1)-ymed(j));

        gL=0.5*(gLm+gLp);

        gRm=(erad[lzp]-erad[lzp-pitch])/
            (ymed(j)-ymed(j-1));

        gRp=(erad[lzp+pitch]-erad[lzp])/
            (ymed(j+1)-ymed(j));

        gR=0.5*(gRm+gRp);

        gt2=wL*gL+wR*gR;
#endif

        gradmag=sqrt(gn*gn+gt1*gt1+gt2*gt2);

        R=gradmag/(chi*Eface);

        if (R < 0.0)
          R=0.0;

        if (R < 1e8)
          lambda=(2.0+R)/(6.0+3.0*R+R*R);
        else
          lambda=1.0/R;

        Dface=clight*lambda/chi;

        dist=ymed(j)*den;

        if (dist > 0.0)
          rate+=Dface*SurfZ(i,j,k+1)/dist;

#endif


        /*
         * Conservative explicit diffusion stability rate:
         *
         * Gamma_i =
         *
         *   1/V_i sum_f D_f A_f/d_f
         */
        rate*=InvVol(i,j,k);

        if (rate > 0.0)
          dtfield[l]=0.20/rate;
        else
          dtfield[l]=1e30;

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