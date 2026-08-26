#include "fargo3d.h"

#ifdef RADIATION

/*
 * ================================================================
 * Temporary FLD validation parameters
 * ================================================================
 */

#ifndef RADFLD_THICK_RHO
#define RADFLD_THICK_RHO 1.0
#endif

#ifndef RADFLD_THICK_KAPPA
#define RADFLD_THICK_KAPPA 1.0e6
#endif

#ifndef RADFLD_THIN_RHO
#define RADFLD_THIN_RHO 1.0
#endif

#ifndef RADFLD_THIN_KAPPA
#define RADFLD_THIN_KAPPA 1.0e-4
#endif

#ifndef RADFLD_THIN_SLOPE
#define RADFLD_THIN_SLOPE 1.0
#endif

#ifndef RADFLD_CONS_RHO
#define RADFLD_CONS_RHO 1.0
#endif

#ifndef RADFLD_CONS_KAPPA
#define RADFLD_CONS_KAPPA 1.0
#endif

#ifndef RADFLD_CONS_AMP
#define RADFLD_CONS_AMP 0.1
#endif

#ifndef RADFLD_CONS_SIGMA
#define RADFLD_CONS_SIGMA 0.05
#endif

#ifndef RADFLD_CONS_RC
#define RADFLD_CONS_RC 1.0
#endif

#ifndef RADFLD_CLIGHT_TEST
#define RADFLD_CLIGHT_TEST 1.0e-5
#endif

/*
 * Physical speed of light in FARGO3D code units:
 *
 * v_unit = sqrt(G_CGS MSTAR_CGS / R0_CGS)
 * c_code = c_CGS / v_unit
 */
#ifndef RAD_CLIGHT_CODE
#define RAD_CLIGHT_CODE (2.99792458e10/sqrt(G_CGS*MSTAR_CGS/R0_CGS))
#endif

#ifdef RADFLDTHICKTEST
#ifndef RADFLD_KAPPA_TEST
#define RADFLD_KAPPA_TEST RADFLD_THICK_KAPPA
#endif
#elif defined(RADFLDTHINTEST)
#ifndef RADFLD_KAPPA_TEST
#define RADFLD_KAPPA_TEST RADFLD_THIN_KAPPA
#endif
#elif defined(RADFLDCONSERVTEST)
#ifndef RADFLD_KAPPA_TEST
#define RADFLD_KAPPA_TEST RADFLD_CONS_KAPPA
#endif
#else
#ifndef RADFLD_KAPPA_TEST
#define RADFLD_KAPPA_TEST 1.0
#endif
#endif

/*
 * ================================================================
 * Manufactured variable-D validation parameters
 * ================================================================
 */

#ifndef RADFLD_VAR_AMP
#define RADFLD_VAR_AMP 1.0
#endif

#ifndef RADFLD_VAR_D0
#define RADFLD_VAR_D0 1.0e-4
#endif

#ifndef RADFLD_VAR_DT
#define RADFLD_VAR_DT 1.0e-3
#endif

/*
 * ================================================================
 * Validation-mode consistency
 * ================================================================
 */

#if (defined(RADFLDVARTEST) + defined(RADFLDOPTEST) + defined(RADFLDFACETEST)) > 1
#error "RADFLDVARTEST, RADFLDOPTEST and RADFLDFACETEST are mutually exclusive."
#endif

#if (defined(RADFLDTHICKTEST) + defined(RADFLDTHINTEST) + defined(RADFLDCONSERVTEST)) > 1
#error "RADFLDTHICKTEST, RADFLDTHINTEST and RADFLDCONSERVTEST are mutually exclusive."
#endif

#if (defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST) || defined(RADFLDCONSERVTEST)) && !defined(RADFLDFACETEST)
#error "Dynamic FLD validation tests require RADFLDFACETEST."
#endif

/*
 * ================================================================
 * Dynamic FLD explicit timestep
 * ================================================================
 */

real RadiationFLDFaceDt(real clight,Field* Erad,Field* Rho,Field* KappaR) {
  real dtlocal,dtglobal;

  if (clight <= 0.0) return 1e30;

  FARGO_SAFE(RadiationFLDDtField(clight,Erad,Rho,KappaR,RadDt));

  dtlocal=reduction_full_MIN(RadDt,NGHY,Ny+NGHY,NGHZ,Nz+NGHZ);

#ifdef PARALLEL
#ifdef FLOAT
  MPI_Allreduce(&dtlocal,&dtglobal,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
#else
  MPI_Allreduce(&dtlocal,&dtglobal,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
#else
  dtglobal=dtlocal;
#endif

  if (dtglobal <= 0.0 || !isfinite(dtglobal)) return 1e30;

  return dtglobal;
}

/*
 * ================================================================
 * Validated constant-D timestep
 * ================================================================
 */

real RadiationDiffusionDt(real diffcoef) {
  int i,j,k;
  real dx,dy,dz,rate;
  real ratemax_local=0.0,ratemax_global=0.0;

  if (diffcoef <= 0.0) return 1e30;

#ifdef Z
  for (k=NGHZ;k<Nz+NGHZ;k++) {
#else
  k=0;
#endif

#ifdef Y
    for (j=NGHY;j<Ny+NGHY;j++) {
#else
    j=0;
#endif

#ifdef X
      for (i=NGHX;i<Nx+NGHX;i++) {
#else
      i=0;
#endif

        rate=0.0;

#ifdef X
        dx=zone_size_x(i,j,k);
        if (dx > 0.0) rate+=1.0/(dx*dx);
#endif

#ifdef Y
        dy=zone_size_y(j,k);
        if (dy > 0.0) rate+=1.0/(dy*dy);
#endif

#ifdef Z
        dz=zone_size_z(j,k);
        if (dz > 0.0) rate+=1.0/(dz*dz);
#endif

        rate*=diffcoef;

        if (rate > ratemax_local) ratemax_local=rate;

#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif

#ifdef PARALLEL
#ifdef FLOAT
  MPI_Allreduce(&ratemax_local,&ratemax_global,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
#else
  MPI_Allreduce(&ratemax_local,&ratemax_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
#else
  ratemax_global=ratemax_local;
#endif

  if (ratemax_global <= 0.0) return 1e30;

  return 0.20/ratemax_global;
}

/*
 * ================================================================
 * Radiation diffusion driver
 * ================================================================
 */

void RadiationDiffusion(real dt) {
  int n,nsub,gas_index=-1;
  real dtdiff,dtsub;
  real clight;

  if (dt <= 0.0) return;

#if !defined(RADFLDVARTEST) && !defined(RADFLDFACETEST)
  if (RADDIFFCOEF <= 0.0) return;
#endif

  /*
   * Restore the gas fluid explicitly.
   */
  for (n=0;n<NFLUIDS;n++) {
    if (Fluids[n]->Fluidtype == GAS) {
      gas_index=n;
      break;
    }
  }

  if (gas_index < 0)
    prs_error("Radiation diffusion requires a gas fluid.");

  FluidIndex=gas_index;
  SelectFluid(FluidIndex);

  /*
   * ==============================================================
   * Select radiation propagation speed
   * ==============================================================
   *
   * Controlled validation tests retain their artificial test speed.
   * The general FLD solver uses the physical speed of light expressed
   * in FARGO3D code units.
   */
#if defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST) || defined(RADFLDCONSERVTEST)
  clight=RADFLD_CLIGHT_TEST;
#else
  clight=RAD_CLIGHT_CODE;
#endif

  /*
   * ==============================================================
   * Controlled dynamic-FLD initial conditions
   * ==============================================================
   */

#if defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST) || defined(RADFLDCONSERVTEST)

  if (Timestepcount == 0) {

#ifdef RADFLDTHICKTEST
    FARGO_SAFE(
      RadiationSetThickTest(
        RADFLD_THICK_RHO,
        1.0e-3,
        0.10,
        1.0,
        Energyrad,
        Density
      )
    );
#endif

#ifdef RADFLDTHINTEST
    FARGO_SAFE(
      RadiationSetThinTest(
        RADFLD_THIN_RHO,
        RADFLD_THIN_SLOPE,
        1.0,
        Energyrad,
        Density
      )
    );
#endif

#ifdef RADFLDCONSERVTEST
    FARGO_SAFE(
      RadiationSetThickTest(
        RADFLD_CONS_RHO,
        RADFLD_CONS_AMP,
        RADFLD_CONS_SIGMA,
        RADFLD_CONS_RC,
        Energyrad,
        Density
      )
    );
#endif

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));
  }

#endif

  /*
   * ==============================================================
   * Rosseland mean opacity
   * ==============================================================
   *
   * Controlled validation benchmarks retain their prescribed
   * constant opacity. RADOPACITYTEST activates the physical Semenov
   * Rosseland opacity evaluated from the current gas density and
   * internal energy.
   */

#if defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST) || defined(RADFLDCONSERVTEST)

  FARGO_SAFE(
    RadiationSetDiffConstant(
      RADFLD_KAPPA_TEST,
      RadKappaR
    )
  );

#elif defined(RADOPACITYTEST)

  FARGO_SAFE(FillGhosts(PrimitiveVariables()));

  FARGO_SAFE(
    RadiationOpacityTable(
      Density,
      Energy,
      RadKappaR
    )
  );

#else

  FARGO_SAFE(
    RadiationSetDiffConstant(
      RADFLD_KAPPA_TEST,
      RadKappaR
    )
  );

#endif

  /*
   * ==============================================================
   * Manufactured variable-D test
   * ==============================================================
   */

#ifdef RADFLDVARTEST

  FARGO_SAFE(
    RadiationSetVariableTest(
      RADFLD_VAR_AMP,
      RADFLD_VAR_D0,
      Energyrad,
      RadDiff
    )
  );

  FARGO_SAFE(
    RadiationDiffusionFLDStep(
      RADFLD_VAR_DT,
      Energyrad,
      RadDiff,
      EnergyradNew
    )
  );

  FARGO_SAFE(copy_field(Energyrad,EnergyradNew));

  return;

#endif

  /*
   * ==============================================================
   * Dynamic face-centered FLD
   * ==============================================================
   */

#ifdef RADFLDFACETEST

#ifdef RADFLDCONSERVTEST

  if (Timestepcount == 0) {

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    dtdiff=RadiationFLDFaceDt(
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density,
      RadKappaR
    );

    if (dtdiff <= 0.0 || !isfinite(dtdiff))
      prs_error("Invalid FLD conservation-test timestep.");

    dtsub=0.10*dtdiff;

    if (dtsub <= 0.0 || !isfinite(dtsub))
      prs_error("Invalid FLD conservation-test substep.");

    if (CPU_Master) {
      printf(
        "FLD_CONSERV_SINGLE_STEP step=%d time=%.12e "
        "dt_global=%.12e dt_fld=%.12e dt_test=%.12e\n",
        Timestepcount,
        (double)PhysicalTime,
        (double)dt,
        (double)dtdiff,
        (double)dtsub
      );
      fflush(stdout);
    }

    FARGO_SAFE(
      RadiationFLDFaceStep(
        dtsub,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        RadKappaR,
        EnergyradNew
      )
    );

    FARGO_SAFE(copy_field(Energyrad,EnergyradNew));

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    FARGO_SAFE(
      RadiationFLDFields(
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        RadKappaR,
        RadR,
        RadLambda,
        RadDiff
      )
    );

    if (CPU_Master) {
      printf("FLD_CONSERV_SINGLE_STEP completed.\n");
      fflush(stdout);
    }
  }

  return;

#endif

#ifdef RADFLDTHINTEST

  if (Timestepcount == 0) {

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    dtdiff=RadiationFLDFaceDt(
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density,
      RadKappaR
    );

    if (dtdiff <= 0.0 || !isfinite(dtdiff))
      prs_error("Invalid FLD thin-test timestep.");

    dtsub=0.10*dtdiff;

    if (dtsub <= 0.0 || !isfinite(dtsub))
      prs_error("Invalid FLD thin-test substep.");

    if (CPU_Master) {
      printf(
        "FLD_THIN_SINGLE_STEP step=%d time=%.12e "
        "dt_global=%.12e dt_fld=%.12e dt_test=%.12e\n",
        Timestepcount,
        (double)PhysicalTime,
        (double)dt,
        (double)dtdiff,
        (double)dtsub
      );
      fflush(stdout);
    }

    FARGO_SAFE(
      RadiationFLDFaceStep(
        dtsub,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        RadKappaR,
        EnergyradNew
      )
    );

    FARGO_SAFE(copy_field(Energyrad,EnergyradNew));

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    FARGO_SAFE(
      RadiationFLDFields(
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        RadKappaR,
        RadR,
        RadLambda,
        RadDiff
      )
    );

    if (CPU_Master) {
      printf("FLD_THIN_SINGLE_STEP completed.\n");
      fflush(stdout);
    }
  }

  return;

#endif

  /*
   * General dynamic face-centered FLD evolution using RKL2.
   *
   * RadiationFLDFaceDt provides the forward-Euler stability limit.
   * RadiationRKL2 then chooses the minimum number of Legendre stages
   * required to cover the full hydrodynamic timestep.
   */
  {
    int rkl_stages;

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    dtdiff=RadiationFLDFaceDt(
      clight,
      Energyrad,
      Density,
      RadKappaR
    );

    if (dtdiff <= 0.0 || !isfinite(dtdiff))
      prs_error("Invalid FLD forward-Euler timestep.");

    rkl_stages=RadiationRKL2(
      dt,
      dtdiff,
      clight,
      Energyrad,
      Density,
      RadKappaR
    );

    if (rkl_stages <= 0)
      prs_error("Invalid RKL2 stage count.");

    if (CPU_Master &&
        (Timestepcount < 10 || Timestepcount%100 == 0)) {

      printf(
        "FLD_RKL2_DIAG step=%d time=%.10e "
        "dt=%.10e dtfe=%.10e ratio=%.4e "
        "clight=%.10e stages=%d\n",
        Timestepcount,
        (double)PhysicalTime,
        (double)dt,
        (double)dtdiff,
        (double)(dt/dtdiff),
        (double)clight,
        rkl_stages
      );

      fflush(stdout);
    }
  }

  FARGO_SAFE(FillGhosts(PrimitiveVariables()));

  FARGO_SAFE(
    RadiationFLDFields(
      clight,
      Energyrad,
      Density,
      RadKappaR,
      RadR,
      RadLambda,
      RadDiff
    )
  );

  return;

#endif

  /*
   * ==============================================================
   * Validated constant-D paths
   * ==============================================================
   */

  dtdiff=RadiationDiffusionDt(RADDIFFCOEF);

  if (dtdiff <= 0.0 || !isfinite(dtdiff))
    prs_error("Invalid constant-D radiation diffusion timestep.");

  nsub=(int)ceil(dt/dtdiff);

  if (nsub < 1)
    nsub=1;

  dtsub=dt/(real)nsub;

  for (n=0;n<nsub;n++) {

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

#ifdef RADFLDOPTEST

    FARGO_SAFE(
      RadiationSetDiffConstant(
        RADDIFFCOEF,
        RadDiff
      )
    );

    FARGO_SAFE(
      RadiationDiffusionFLDStep(
        dtsub,
        Energyrad,
        RadDiff,
        EnergyradNew
      )
    );

#else

    /*
     * Keep the historical validation light speed in this legacy
     * constant-D path. This branch is not the physical dynamic
     * face-centered FLD solver.
     */
    FARGO_SAFE(
      RadiationFLDFields(
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        RadKappaR,
        RadR,
        RadLambda,
        RadDiff
      )
    );

    FARGO_SAFE(
      RadiationDiffusionStep(
        dtsub,
        RADDIFFCOEF,
        Energyrad,
        EnergyradNew
      )
    );

#endif

    FARGO_SAFE(copy_field(Energyrad,EnergyradNew));
  }

  FARGO_SAFE(FillGhosts(PrimitiveVariables()));

  FARGO_SAFE(
    RadiationFLDFields(
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density,
      RadKappaR,
      RadR,
      RadLambda,
      RadDiff
    )
  );

#ifdef RADFLDOPTEST

  FARGO_SAFE(
    RadiationSetDiffConstant(
      RADDIFFCOEF,
      RadDiff
    )
  );

#endif
}

#endif