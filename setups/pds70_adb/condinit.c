#include "fargo3d.h"

void Init() {
  int i,j,k;
  real *v1,*v2,*v3,*e,*rho;

#ifdef RADIATION
  real *erad;
#endif

  real r,omega,h_over_r,h_geom,keplerian_velocity;
  real xi,beta;

#ifdef RADDIFFTEST
  const real rad_rc = 1.0;
  const real rad_sigma0 = 0.05;
  const real rad_amp = 1.0;
#endif

  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;


#ifdef RADIATION
  EnergyradNew = CreateField("energyradnew",0,0,0,0);
  RadR         = CreateField("radr",0,0,0,0);
  RadLambda    = CreateField("radlambda",0,0,0,0);
  RadDiff      = CreateField("raddiff",0,0,0,0);
#endif


#ifdef RADIATION
  OUTPUT(Energyrad);
  OUTPUT(RadR);
  OUTPUT(RadLambda);
  OUTPUT(RadDiff);
  erad = Energyrad->field_cpu;
#endif

  /*
   * FARGO3D linear index l is defined internally from i,j,k.
   * Do not declare or assign l manually.
   */
  for (k=0;k<Nz+2*NGHZ;k++) {
    for (j=0;j<Ny+2*NGHY;j++) {
      for (i=0;i<Nx+2*NGHX;i++) {

        r = Ymed(j);

        /* ----------------------------------------------------
         * Hydrodynamic initial condition
         * ---------------------------------------------------- */

        v2[l] = 0.0;
        v3[l] = 0.0;

        omega = sqrt(G*MSTAR/(r*r*r));
        v1[l] = omega*r;

#ifndef CYLINDRICAL
        xi = SIGMASLOPE + 1.0 + FLARINGINDEX;
        beta = 1.0 - 2.0*FLARINGINDEX;

        h_over_r = ASPECTRATIO*pow(r/R0,FLARINGINDEX);
        h_geom = h_over_r*r;

        if (FLARINGINDEX == 0.0) {
          rho[l] =
            SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*
            pow(r/R0,-xi)*
            pow(
              sin(Zmed(k)),
              -beta-xi+1.0/(h_geom*h_geom)
            );
        } else {
          rho[l] =
            SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*
            pow(r/R0,-xi)*
            pow(sin(Zmed(k)),-xi-beta)*
            exp(
              (1.0-pow(sin(Zmed(k)),-2.0*FLARINGINDEX))/
              (2.0*FLARINGINDEX*h_over_r*h_over_r)
            );
        }

        v1[l] *= sqrt(
          pow(sin(Zmed(k)),-2.0*FLARINGINDEX)
          -(beta+xi)*h_over_r*h_over_r
        );

        v1[l] -= OMEGAFRAME*r*sin(Zmed(k));
#endif

#ifdef ISOTHERMAL
        keplerian_velocity = sqrt(G*MSTAR/r);
        e[l] =
          ASPECTRATIO*
          pow(r/R0,FLARINGINDEX)*
          keplerian_velocity;
#else
        e[l] =
          rho[l]*
          h_geom*h_geom*
          omega*omega/
          (GAMMA-1.0);
#endif

        /* ----------------------------------------------------
         * Radiation field
         * ---------------------------------------------------- */

#ifdef RADIATION

#ifdef RADDIFFTEST
        /*
         * Spherically symmetric diffusion test.
         *
         * For u(r,t)=r E_r(r,t), the spherical diffusion
         * equation reduces to the Cartesian 1-D diffusion
         * equation. Therefore initialize
         *
         *   u(r,0) = A exp[-(r-r_c)^2/(2 sigma_0^2)]
         *
         * which implies
         *
         *   E_r(r,0) = u(r,0)/r.
         */
        {
          real dr;

          dr = r-rad_rc;

          erad[l] =
            (rad_amp/r)*
            exp(
              -0.5*dr*dr/
              (rad_sigma0*rad_sigma0)
            );
        }

        /*
         * Disable advection during the pure diffusion test.
         */
        v1[l] = 0.0;
        v2[l] = 0.0;
        v3[l] = 0.0;

#else
        /*
         * Default radiation initialization outside RADDIFFTEST.
         * Retains the previous azimuthal test profile.
         */
        {
          real dphi = Xmed(i);
          real sigma = 0.15;

          erad[l] =
            1.0+
            exp(
              -0.5*dphi*dphi/
              (sigma*sigma)
            );
        }
#endif

#endif /* RADIATION */

      }
    }
  }
}

void CondInit() {
  Fluids[0] = CreateFluid("gas",GAS);
  SelectFluid(0);
  Init();
}