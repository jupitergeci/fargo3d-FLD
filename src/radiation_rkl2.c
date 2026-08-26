#include "fargo3d.h"

#ifdef RADIATION

#ifndef RADRKL2_MAX_STAGES
#define RADRKL2_MAX_STAGES 10000
#endif

static real RadiationRKL2b(int j) {
  real jr;

  if (j < 2)
    return 1.0/3.0;

  jr=(real)j;

  return (jr*jr+jr-2.0)/(2.0*jr*(jr+1.0));
}

int RadiationRKL2Stages(real dt,real dtfe) {
  int s;
  real ratio;
  real capacity;

  if (dt <= 0.0)
    return 0;

  if (dtfe <= 0.0 || !isfinite(dtfe))
    prs_error("RadiationRKL2Stages received an invalid forward-Euler timestep.");

  if (dt <= dtfe)
    return 1;

  ratio=dt/dtfe;

  s=(int)ceil(0.5*(-1.0+sqrt(9.0+16.0*ratio)));

  if (s < 2)
    s=2;

  while (1) {
    real sr=(real)s;

    capacity=0.25*(sr*sr+sr-2.0);

    if (capacity*dtfe >= dt)
      break;

    s++;

    if (s > RADRKL2_MAX_STAGES)
      prs_error("RKL2 requires too many stages.");
  }

  while (s > 2) {
    real sm1=(real)(s-1);

    capacity=0.25*(sm1*sm1+sm1-2.0);

    if (capacity*dtfe < dt)
      break;

    s--;
  }

  return s;
}

int RadiationRKL2(real dt,real dtfe,real clight,Field* Erad,Field* Rho,Field* KappaR) {
  int j;
  int stages;

  real sreal;
  real w1;

  real bj;
  real bjm1;
  real bjm2;
  real ajm1;

  real mu;
  real nu;
  real c0;
  real mut;
  real gamt;
  real mut1;

  if (dt <= 0.0)
    return 0;

  if (dtfe <= 0.0 || !isfinite(dtfe))
    prs_error("RadiationRKL2 received an invalid forward-Euler timestep.");

  if (clight <= 0.0 || !isfinite(clight))
    prs_error("RadiationRKL2 received an invalid speed of light.");

  stages=RadiationRKL2Stages(dt,dtfe);

  /*
   * A single forward-Euler step is already stable.
   */
  if (stages == 1) {

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    FARGO_SAFE(
      RadiationFLDFaceStep(
        dt,
        clight,
        Erad,
        Rho,
        KappaR,
        EnergyradNew
      )
    );

    FARGO_SAFE(copy_field(Erad,EnergyradNew));

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    return 1;
  }

  /*
   * RKL2 stability parameter.
   */
  sreal=(real)stages;
  w1=4.0/(sreal*sreal+sreal-2.0);

  /*
   * Y0 = E_R^n.
   *
   * RadStage1 initially stores Y0 and becomes Y_{j-2}
   * during the recurrence.
   */
  FARGO_SAFE(copy_field(RadStage0,Erad));
  FARGO_SAFE(copy_field(RadStage1,Erad));

  /*
   * L0 = L(Y0).
   *
   * RadDt is no longer required after dt_FE has been found,
   * so it is reused to store L(Y0).
   */
  FARGO_SAFE(FillGhosts(PrimitiveVariables()));

  FARGO_SAFE(
    RadiationFLDOperator(
      clight,
      Erad,
      Rho,
      KappaR,
      RadDt
    )
  );

  /*
   * First RKL2 stage:
   *
   * Y1 = Y0 + b1*w1*dt*L0.
   */
  mut1=(1.0/3.0)*w1;

  FARGO_SAFE(
    RadiationRKL2Update(
      0.0,
      0.0,
      1.0,
      mut1,
      0.0,
      dt,
      RadStage0,
      RadStage0,
      RadStage0,
      RadDt,
      RadDt,
      EnergyradNew
    )
  );

  FARGO_SAFE(copy_field(Erad,EnergyradNew));

  /*
   * General RKL2 recurrence.
   */
  for (j=2;j<=stages;j++) {

    /*
     * Erad      = Y_{j-1}
     * RadStage1 = Y_{j-2}
     * RadStage0 = Y0
     * RadDt     = L(Y0)
     */

    FARGO_SAFE(FillGhosts(PrimitiveVariables()));

    /*
     * RadDiff is temporarily reused as L(Y_{j-1}).
     */
    FARGO_SAFE(
      RadiationFLDOperator(
        clight,
        Erad,
        Rho,
        KappaR,
        RadDiff
      )
    );

    bj=RadiationRKL2b(j);
    bjm1=RadiationRKL2b(j-1);
    bjm2=RadiationRKL2b(j-2);

    ajm1=1.0-bjm1;

    mu=((2.0*(real)j-1.0)/(real)j)*(bj/bjm1);

    nu=-(((real)j-1.0)/(real)j)*(bj/bjm2);

    c0=1.0-mu-nu;

    mut=mu*w1;

    gamt=-ajm1*mut;

    /*
     * Y_j -> EnergyradNew.
     */
    FARGO_SAFE(
      RadiationRKL2Update(
        mu,
        nu,
        c0,
        mut,
        gamt,
        dt,
        Erad,
        RadStage1,
        RadStage0,
        RadDiff,
        RadDt,
        EnergyradNew
      )
    );

    /*
     * Preserve Y_{j-1} as Y_{j-2} for the next recurrence.
     */
    FARGO_SAFE(copy_field(RadStage1,Erad));

    /*
     * Current radiation state becomes Y_j.
     */
    FARGO_SAFE(copy_field(Erad,EnergyradNew));
  }

  FARGO_SAFE(FillGhosts(PrimitiveVariables()));

  return stages;
}

#endif