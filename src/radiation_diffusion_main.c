#include "fargo3d.h"

#ifdef RADIATION

/*
 * Temporary FLD validation parameters.
 *
 * These values are expressed in code units and are used only during
 * the current dynamic-FLD validation stage. They are not yet the
 * physical Rosseland opacity or physical speed of light.
 */
#ifdef RADFLDTHICKTEST

#ifndef RADFLD_KAPPA_TEST
#define RADFLD_KAPPA_TEST RADFLD_THICK_KAPPA
#endif

#elif defined(RADFLDTHINTEST)

#ifndef RADFLD_KAPPA_TEST
#define RADFLD_KAPPA_TEST RADFLD_THIN_KAPPA
#endif

#else

#ifndef RADFLD_KAPPA_TEST
#define RADFLD_KAPPA_TEST 1.0
#endif

#endif

#ifndef RADFLD_THICK_RHO
#define RADFLD_THICK_RHO 1.0
#endif

#ifndef RADFLD_THICK_KAPPA
#define RADFLD_THICK_KAPPA 1.0e6
#endif

#ifndef RADFLD_CLIGHT_TEST
#define RADFLD_CLIGHT_TEST 1.0e-5
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


/*
 * Manufactured variable-D validation parameters.
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
 * Only one FLD validation path may be active at a time.
 */
#if (defined(RADFLDVARTEST) + defined(RADFLDOPTEST) + defined(RADFLDFACETEST)) > 1
#error "RADFLDVARTEST, RADFLDOPTEST and RADFLDFACETEST are mutually exclusive."
#endif

#if defined(RADFLDTHICKTEST) && defined(RADFLDTHINTEST)
#error "RADFLDTHICKTEST and RADFLDTHINTEST cannot be enabled simultaneously."
#endif

#if (defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST)) && !defined(RADFLDFACETEST)
#error "RADFLDTHICKTEST and RADFLDTHINTEST require RADFLDFACETEST."
#endif
/*
 * Explicit FLD stability timestep.
 *
 * RadiationFLDDtField computes, for every active cell,
 *
 *   dt_i = Cdiff / Gamma_i
 *
 * with
 *
 *   Gamma_i = (1/V_i) sum_f (D_f A_f / d_f),
 *
 * and Cdiff = 0.20.
 *
 * RadiationFLDDtField is executed on the selected architecture.
 * In the GPU configuration, Energyrad and Density remain on the
 * device and only the reduced timestep information is transferred
 * back to the host.
 */
real RadiationFLDFaceDt(real kappa,real clight,Field* Erad,Field* Rho) {
  real dtlocal;
  real dtglobal;

  if (kappa <= 0.0 || clight <= 0.0)
    return 1e30;

  /*
   * Compute the cell-wise FLD timestep field.
   *
   * On GPU builds this operates directly on device-resident
   * Energyrad and Density.
   */
  FARGO_SAFE(
    RadiationFLDDtField(
      kappa,
      clight,
      Erad,
      Rho,
      RadDt
    )
  );

  /*
   * Native FARGO3D reduction.
   *
   * reduction_full_MIN performs the first reduction stage using the
   * selected Reduction architecture. With Reduction=GPU, the large
   * 3D field is reduced on the device before the small reduced array
   * is inspected on the host.
   */
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
 * Validated constant-coefficient diffusion timestep.
 *
 * This function is retained unchanged for:
 *
 *   1. the original constant-D diffusion solver;
 *   2. RADFLDOPTEST.
 *
 * RADFLDFACETEST does not use this timestep.
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


void RadiationDiffusion(real dt) {
  int n,nsub,gas_index=-1;
  real dtdiff,dtsub;

  if (dt <= 0.0)
    return;


  /*
   * The original constant-D paths require RADDIFFCOEF.
   *
   * RADFLDVARTEST and RADFLDFACETEST construct their own spatially
   * dependent diffusion coefficients and therefore do not require
   * RADDIFFCOEF.
   */
#if !defined(RADFLDVARTEST) && !defined(RADFLDFACETEST)

  if (RADDIFFCOEF <= 0.0)
    return;

#endif


  /*
   * Restore the gas fluid explicitly.
   *
   * MULTIFLUID() may leave FluidIndex=NFLUIDS after completion.
   * FillGhosts() and radiation routines must therefore operate with
   * the gas fluid explicitly selected.
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

// #ifdef RADFLDTHICKTEST

//   if (Timestepcount == 0) {
//     FARGO_SAFE(
//       RadiationSetThickTest(
//         RADFLD_THICK_RHO,
//         1.0e-3,
//         0.10,
//         1.0,
//         Energyrad,
//         Density
//       )
//     );

//     FARGO_SAFE(
//       FillGhosts(
//         PrimitiveVariables()
//       )
//     );
//   }

// #endif

#if defined(RADFLDTHICKTEST) || defined(RADFLDTHINTEST)

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

    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );
  }

#endif


#ifdef RADFLDVARTEST

  /*
   * Manufactured variable-D operator test.
   *
   * E_R(r) = A r^2
   * D(r)   = D0 r
   *
   * Exact spherical result:
   *
   * div(D grad E_R) = 8 A D0 r.
   *
   * RadiationSetVariableTest fills both the active domain and the
   * required ghost zones explicitly for this controlled test.
   */
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


#ifdef RADFLDFACETEST

#ifdef RADFLDTHINTEST

  /*
   * Optically thin single-step validation.
   *
   * The thin benchmark is intentionally evolved only once.
   * Continuing the explicit FLD evolution for many global
   * timesteps causes the profile to flatten, D_FLD to increase,
   * and the parabolic explicit timestep to become progressively
   * restrictive. That behavior is not required to validate the
   * free-streaming limit.
   */
  if (Timestepcount == 0) {

    /*
     * Energyrad and Density must have valid halos because both
     * the multidimensional face limiter and the timestep use
     * neighboring cells.
     */
    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );

    /*
     * Compute the stable explicit timestep for the initial
     * optically thin state.
     */
    dtdiff=RadiationFLDFaceDt(
      RADFLD_KAPPA_TEST,
      RADFLD_CLIGHT_TEST,
      Energyrad,
      Density
    );

    if (dtdiff <= 0.0 || !isfinite(dtdiff))
      prs_error("Invalid FLD thin-test timestep.");

    /*
     * Take only 10% of the explicit stability limit.
     *
     * This is deliberately small: the purpose is to verify the
     * thin-limit FLD flux, not to evolve the profile over a long
     * radiation-transport time.
     */
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
     * Synchronize the resulting radiation state.
     */
    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );

    /*
     * Evaluate the already validated cell-centered diagnostics
     * on the post-step radiation state.
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
      printf("FLD_THIN_SINGLE_STEP completed.\n");
      fflush(stdout);
    }
  }

  /*
   * For Timestepcount > 0 the radiation state is intentionally
   * frozen. FARGO3D may continue advancing PhysicalTime until the
   * next requested output, but no additional thin-test diffusion
   * is applied.
   */
  return;

#endif


  /*
   * General dynamic face-FLD evolution.
   *
   * This branch is used when RADFLDTHINTEST is not active.
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
          "FLD_PERF_DIAG step=%d time=%.10e dt=%.10e dtfld=%.10e ratio=%.4e\n",
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
        prs_error("Too many explicit FLD diffusion substeps.");
    }
  }

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
   * ================================================================
   * Validated constant-D paths
   * ================================================================
   */

  dtdiff=RadiationDiffusionDt(RADDIFFCOEF);

  if (dtdiff <= 0.0 || !isfinite(dtdiff))
    prs_error("Invalid constant-D radiation diffusion timestep.");

  nsub=(int)ceil(dt/dtdiff);

  if (nsub < 1)
    nsub=1;

  dtsub=dt/(real)nsub;


  for (n=0;n<nsub;n++) {

    /*
     * Energyrad and Density require valid halo values before
     * gradients or FLD diagnostics are evaluated.
     */
    FARGO_SAFE(
      FillGhosts(
        PrimitiveVariables()
      )
    );


#ifdef RADFLDOPTEST

    /*
     * Constant-D regression of the conservative variable-D operator.
     *
     * The complete local RadDiff field is forced to RADDIFFCOEF so
     * that RadiationDiffusionFLDStep must reduce exactly to the
     * validated constant-coefficient solver.
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
     * Compute the already validated cell-centered quantities:
     *
     *   R
     *   lambda
     *   D_FLD
     *
     * The evolution itself remains the original constant-D solver in
     * this branch.
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
   * Synchronize the final radiation state.
   */
  FARGO_SAFE(
    FillGhosts(
      PrimitiveVariables()
    )
  );


  /*
   * Recompute diagnostics from the final E_R state so that all output
   * fields correspond to the same radiation snapshot.
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
   * During RADFLDOPTEST the actual evolved coefficient was the fixed
   * RADDIFFCOEF rather than the diagnostic FLD coefficient.
   *
   * Restore RadDiff before output so raddiff*.dat records the
   * coefficient used by the regression operator.
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