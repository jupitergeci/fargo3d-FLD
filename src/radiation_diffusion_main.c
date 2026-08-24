#include "fargo3d.h"

#ifdef RADIATION

real RadiationDiffusionDt(real diffcoef) {
  int i,j,k;
  real dx,dy,dz,rate;
  real ratemax_local=0.0;
  real ratemax_global=0.0;

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

void RadiationDiffusion(real dt) {
  int n,nsub,gas_index=-1;
  real dtdiff,dtsub;

  if (dt <= 0.0 || RADDIFFCOEF <= 0.0) return;

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

  dtdiff=RadiationDiffusionDt(RADDIFFCOEF);

  nsub=(int)ceil(dt/dtdiff);
  if (nsub < 1) nsub=1;

  dtsub=dt/(real)nsub;

  for (n=0;n<nsub;n++) {
    FARGO_SAFE(FillGhosts(ENERGYRAD));
    FARGO_SAFE(RadiationDiffusionStep(dtsub,RADDIFFCOEF,
                                      Energyrad,EnergyradNew));
    FARGO_SAFE(copy_field(Energyrad,EnergyradNew));
  }

  FARGO_SAFE(FillGhosts(ENERGYRAD));
}

#endif