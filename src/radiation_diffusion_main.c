#include "fargo3d.h"

#ifdef RADIATION

/*
 * ================================================================
 * Temporary FLD validation parameters
 * ================================================================
 *
 * These values are expressed in code units and are used only during
 * the current FLD validation stage. They are not yet the physical
 * Rosseland opacity or physical speed of light.
 */

/*
 * Optically thick test.
 */
#ifndef RADFLD_THICK_RHO
#define RADFLD_THICK_RHO 1.0
#endif

#ifndef RADFLD_THICK_KAPPA
#define RADFLD_THICK_KAPPA 1.0e6
#endif

/*
 * Optically thin test.
 */
#ifndef RADFLD_THIN_RHO
#define RADFLD_THIN_RHO 1.0
#endif

#ifndef RADFLD_THIN_KAPPA
#define RADFLD_THIN_KAPPA 1.0e-4
#endif

#ifndef RADFLD_THIN_SLOPE
#define RADFLD_THIN_SLOPE 1.0
#endif

/*
 * Conservative FLD single-step test.
 */
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

/*
 * Reduced light speed used only for current explicit validation.
 */
#ifndef RADFLD_CLIGHT_TEST
#define RADFLD_CLIGHT_TEST 1.0e-5
#endif

/*
 * Select the opacity appropriate to the active validation problem.
 */
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
 *
 * E_R(r) = A r^2
 * D(r)   = D0 r
 *
 * Exact spherical operator:
 *
 * div(D grad E_R) = 8 A D0 r.
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

/*
 * These are mutually exclusive evolution/operator paths.
 */
#if (defined(RADFLDVARTEST) + defined(RADFLDOPTEST) + defined(RADFLDFACETEST)) > 1
#error "RADFLDVARTEST, RADFLDOPTEST and RADFLDFACETEST are mutually exclusive."
#endif

/*
 * Only one physical face-FLD benchmark may be active at once.
 */
#if (defined(RADFLDTHICKTEST) + defined(RADFLDTHINTEST) + defined(RADFLDCONSERVTEST)) > 1
#error "RADFLDTHICKTEST, RADFLDTHINTEST and RADFLDCONSERVTEST are mutually exclusive."
#endif

/*
 * Thick, thin and conservation tests require the dynamic face solver.
 */
#if (defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST) || defined(RADFLDCONSERVTEST)) && !defined(RADFLDFACETEST)
#error "Dynamic FLD validation tests require RADFLDFACETEST."
#endif


/*
 * ================================================================
 * Dynamic FLD explicit timestep
 * ================================================================
 *
 * RadiationFLDDtField computes
 *
 *   dt_i = Cdiff/Gamma_i
 *
 * where
 *
 *   Gamma_i = (1/V_i) sum_f (D_f A_f/d_f)
 *
 * and Cdiff=0.20.
 *
 * The cell-wise calculation is performed on the selected
 * architecture. In GPU mode the full 3D radiation and density
 * fields remain on the device. FARGO3D's native reduction is then
 * used to obtain the minimum timestep.
 */

real RadiationFLDFaceDt(real kappa,real clight,Field* Erad,Field* Rho) {
  real dtlocal;
  real dtglobal;

  if (kappa <= 0.0 || clight <= 0.0)
    return 1e30;

  FARGO_SAFE(
    RadiationFLDDtField(
      kappa,
      clight,
      Erad,
      Rho,
      RadDt
    )
  );

  dtlocal=reduction_full_MIN(
    RadDt,
    NGHY,
    Ny+NGHY,
    NGHZ,
    Nz+NGHZ
  );

#ifdef PARALLEL

#ifdef FLOAT
  MPI_Allreduce(
    &dtlocal,
    &dtglobal,
    1,
    MPI_FLOAT,
    MPI_MIN,
    MPI_COMM_WORLD
  );
#else
  MPI_Allreduce(
    &dtlocal,
    &dtglobal,
    1,
    MPI_DOUBLE,
    MPI_MIN,
    MPI_COMM_WORLD
  );
#endif

#else

  dtglobal=dtlocal;

#endif

  if (dtglobal <= 0.0 || !isfinite(dtglobal))
    return 1e30;

  return dtglobal;
}


/*
 * ================================================================
 * Validated constant-D timestep
 * ================================================================
 *
 * Retained unchanged for:
 *
 *   1. the original constant-D radiation solver;
 *   2. RADFLDOPTEST.
 */

real RadiationDiffusionDt(real diffcoef) {
  int i,j,k;
  real dx,dy,dz,rate;
  real ratemax_local=0.0;
  real ratemax_global=0.0;

  if (diffcoef <= 0.0)
    return 1e30;

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

        if (dx > 0.0)
          rate+=1.0/(dx*dx);
#endif

#ifdef Y
        dy=zone_size_y(j,k);

        if (dy > 0.0)
          rate+=1.0/(dy*dy);
#endif

#ifdef Z
        dz=zone_size_z(j,k);

        if (dz > 0.0)
          rate+=1.0/(dz*dz);
#endif

        rate*=diffcoef;

        if (rate > ratemax_local)
          ratemax_local=rate;

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
  MPI_Allreduce(
    &ratemax_local,
    &ratemax_global,
    1,
    MPI_FLOAT,
    MPI_MAX,
    MPI_COMM_WORLD
  );
#else
  MPI_Allreduce(
    &ratemax_local,
    &ratemax_global,
    1,
    MPI_DOUBLE,
    MPI_MAX,
    MPI_COMM_WORLD
  );
#endif

#else

  ratemax_global=ratemax_local;

#endif

  if (ratemax_global <= 0.0)
    return 1e30;

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

  if (dt <= 0.0)
    return;

  /*
   * The original constant-D paths require RADDIFFCOEF.
   *
   * RADFLDVARTEST and RADFLDFACETEST construct their own spatially
   * varying coefficients and therefore do not require RADDIFFCOEF.
   */
#if !defined(RADFLDVARTEST) && !defined(RADFLDFACETEST)

  if (RADDIFFCOEF <= 0.0)
    return;

#endif


  /*
   * Restore the gas fluid explicitly.
   *
   * MULTIFLUID() may leave FluidIndex=NFLUIDS after completion.
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
   * Controlled dynamic-FLD initial conditions
   * ==============================================================
   */

#if defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST) || defined(RADFLDCONSERVTEST)

  if (Timestepcount == 0) {

#ifdef RADFLDTHICKTEST

    /*
     * Smooth low-amplitude Gaussian perturbation in a very
     * optically thick medium.
     */
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

    /*
     * Monotonic radial exponential field:
     *
     *   E_R = exp[-a(r-r0)].
     *
     * This produces a controlled optically thin/free-streaming
     * regime with |grad E|/E approximately constant.
     */
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

    /*
     * Localized Gaussian perturbation used for global conservation:
     *
     *   rho = rho0
     *
     *   E_R =
     *     1
     *     + A exp[-(r-rc)^2/(2 sigma^2)].
     *
     * RadiationSetThickTest is reused only as an initializer.
     * The actual opacity for this test is RADFLD_CONS_KAPPA.
     */
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


    /*
     * All face-based limiters require valid halos before the
     * first FLD operator evaluation.
     */
    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );
  }

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

  FARGO_SAFE(
    copy_field(
      Energyrad,
      EnergyradNew
    )
  );

  return;

#endif


  /*
   * ==============================================================
   * Dynamic face-centered FLD
   * ==============================================================
   */

#ifdef RADFLDFACETEST


  /*
   * --------------------------------------------------------------
   * Global conservation single-step test
   * --------------------------------------------------------------
   *
   * The localized Gaussian perturbation is sufficiently far from
   * the domain boundaries that the physical boundary flux is
   * negligible over one small FLD step.
   *
   * The post-step radiation field is frozen after Timestepcount=0,
   * allowing it to be written at output 1 without further evolution.
   */

#ifdef RADFLDCONSERVTEST

  if (Timestepcount == 0) {

    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );

    dtdiff=RadiationFLDFaceDt(
      RADFLD_KAPPA_TEST,
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density
    );

    if (dtdiff <= 0.0 || !isfinite(dtdiff))
      prs_error("Invalid FLD conservation-test timestep.");

    /*
     * Small controlled step.
     */
    dtsub=0.10*dtdiff;

    if (dtsub <= 0.0 || !isfinite(dtsub))
      prs_error("Invalid FLD conservation-test substep.");

    if (CPU_Master) {
      printf(
        "FLD_CONSERV_SINGLE_STEP "
        "step=%d "
        "time=%.12e "
        "dt_global=%.12e "
        "dt_fld=%.12e "
        "dt_test=%.12e\n",
        Timestepcount,
        (double)PhysicalTime,
        (double)dt,
        (double)dtdiff,
        (double)dtsub
      );
      fflush(stdout);
    }

    /*
     * Perform exactly one dynamic FLD update.
     */
    FARGO_SAFE(
      RadiationFLDFaceStep(
        dtsub,
        RADFLD_KAPPA_TEST,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        EnergyradNew
      )
    );

    FARGO_SAFE(
      copy_field(
        Energyrad,
        EnergyradNew
      )
    );

    /*
     * Synchronize the post-step state.
     */
    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );

    /*
     * Recompute cell-centered FLD diagnostics from the final field.
     */
    FARGO_SAFE(
      RadiationFLDFields(
        RADFLD_KAPPA_TEST,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        RadR,
        RadLambda,
        RadDiff
      )
    );

    if (CPU_Master) {
      printf(
        "FLD_CONSERV_SINGLE_STEP completed.\n"
      );
      fflush(stdout);
    }
  }

  /*
   * For all subsequent global timesteps the radiation field remains
   * frozen. FARGO3D continues only until the requested output is
   * written.
   */
  return;

#endif


  /*
   * --------------------------------------------------------------
   * Optically thin single-step test
   * --------------------------------------------------------------
   */

#ifdef RADFLDTHINTEST

  if (Timestepcount == 0) {

    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );

    dtdiff=RadiationFLDFaceDt(
      RADFLD_KAPPA_TEST,
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density
    );

    if (dtdiff <= 0.0 || !isfinite(dtdiff))
      prs_error("Invalid FLD thin-test timestep.");

    dtsub=0.10*dtdiff;

    if (dtsub <= 0.0 || !isfinite(dtsub))
      prs_error("Invalid FLD thin-test substep.");

    if (CPU_Master) {
      printf(
        "FLD_THIN_SINGLE_STEP "
        "step=%d "
        "time=%.12e "
        "dt_global=%.12e "
        "dt_fld=%.12e "
        "dt_test=%.12e\n",
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
        RADFLD_KAPPA_TEST,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        EnergyradNew
      )
    );

    FARGO_SAFE(
      copy_field(
        Energyrad,
        EnergyradNew
      )
    );

    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );

    FARGO_SAFE(
      RadiationFLDFields(
        RADFLD_KAPPA_TEST,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
        RadR,
        RadLambda,
        RadDiff
      )
    );

    if (CPU_Master) {
      printf(
        "FLD_THIN_SINGLE_STEP completed.\n"
      );
      fflush(stdout);
    }
  }

  return;

#endif


  /*
   * --------------------------------------------------------------
   * General dynamic face-centered FLD evolution
   * --------------------------------------------------------------
   *
   * Used, for example, by the thick-limit evolution test.
   */

  {
    real remaining=dt;
    int fld_substep=0;

    while (remaining > 0.0) {

      FARGO_SAFE(
        FillGhosts(
          PrimitiveVariables()
        )
      );

      dtdiff=RadiationFLDFaceDt(
        RADFLD_KAPPA_TEST,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density
      );

      if (CPU_Master && fld_substep == 0) {
        printf(
          "FLD_PERF_DIAG "
          "step=%d "
          "time=%.10e "
          "dt=%.10e "
          "dtfld=%.10e "
          "ratio=%.4e\n",
          Timestepcount,
          (double)PhysicalTime,
          (double)dt,
          (double)dtdiff,
          (double)(dt/dtdiff)
        );
        fflush(stdout);
      }

      if (dtdiff <= 0.0 || !isfinite(dtdiff))
        prs_error("Invalid FLD diffusion timestep.");

      dtsub=(
        remaining < dtdiff
        ? remaining
        : dtdiff
      );

      if (dtsub <= 0.0 || !isfinite(dtsub))
        prs_error("Invalid FLD diffusion substep.");

      FARGO_SAFE(
        RadiationFLDFaceStep(
          dtsub,
          RADFLD_KAPPA_TEST,
          RADFLD_CLIGHT_TEST,
          Energyrad,
          Density,
          EnergyradNew
        )
      );

      FARGO_SAFE(
        copy_field(
          Energyrad,
          EnergyradNew
        )
      );

      remaining-=dtsub;

      if (remaining < 1e-14*dt)
        remaining=0.0;

      fld_substep++;

      if (fld_substep > 10000000)
        prs_error(
          "Too many explicit FLD diffusion substeps."
        );
    }
  }

  /*
   * Synchronize final state and compute diagnostics.
   */
  FARGO_SAFE(
    FillGhosts(
      PrimitiveVariables()
    )
  );

  FARGO_SAFE(
    RadiationFLDFields(
      RADFLD_KAPPA_TEST,
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density,
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

  dtdiff=RadiationDiffusionDt(
    RADDIFFCOEF
  );

  if (dtdiff <= 0.0 || !isfinite(dtdiff))
    prs_error(
      "Invalid constant-D radiation diffusion timestep."
    );

  nsub=(int)ceil(
    dt/dtdiff
  );

  if (nsub < 1)
    nsub=1;

  dtsub=dt/(real)nsub;


  for (n=0;n<nsub;n++) {

    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );


#ifdef RADFLDOPTEST

    /*
     * Constant-D regression of the conservative variable-D operator.
     */
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
     * Standard diagnostic FLD path.
     *
     * The actual evolution in this branch still uses the validated
     * constant-D radiation solver.
     */
    FARGO_SAFE(
      RadiationFLDFields(
        RADFLD_KAPPA_TEST,
        RADFLD_CLIGHT_TEST,
        Energyrad,
        Density,
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


    FARGO_SAFE(
      copy_field(
        Energyrad,
        EnergyradNew
      )
    );
  }


  /*
   * Synchronize final radiation state.
   */
  FARGO_SAFE(
    FillGhosts(
      PrimitiveVariables()
    )
  );


  /*
   * Recompute diagnostics from the final state.
   */
  FARGO_SAFE(
    RadiationFLDFields(
      RADFLD_KAPPA_TEST,
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density,
      RadR,
      RadLambda,
      RadDiff
    )
  );


#ifdef RADFLDOPTEST

  /*
   * In RADFLDOPTEST the actual evolved coefficient is the fixed
   * RADDIFFCOEF. Restore RadDiff before output so the diagnostic file
   * records the coefficient used by that regression.
   */
  FARGO_SAFE(
    RadiationSetDiffConstant(
      RADDIFFCOEF,
      RadDiff
    )
  );

#endif
}

#endif