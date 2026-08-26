//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#include "radiation_opacity_table.h"
//<\INCLUDES>

#ifndef __GPU

void RadiationOpacityTableLoad(void) {
  FILE *fp;
  char line[4096];
  int nrho,ntemp,i,j,n;
  double value,axis_value;
  size_t nbytes=sizeof(real)*SEMENOV_TABLE_SIZE;

  if (SemenovOpacityLoaded) return;

  SemenovLogKappa_cpu=NULL;
#ifdef GPU
  SemenovLogKappa_gpu=NULL;
#endif

  fp=fopen(SEMENOV_OPACITY_FILE,"r");
  if (fp == NULL) {
    if (CPU_Master) fprintf(stderr,"Cannot open Semenov opacity table: %s\n",SEMENOV_OPACITY_FILE);
    prs_error("Cannot open Semenov Rosseland opacity table.");
  }

  for (n=0;n<4;n++) {
    if (fgets(line,sizeof(line),fp) == NULL) {
      fclose(fp);
      prs_error("Invalid Semenov opacity table header.");
    }
  }

  if (fscanf(fp,"%d %d",&nrho,&ntemp) != 2) {
    fclose(fp);
    prs_error("Cannot read Semenov opacity table dimensions.");
  }

  if (nrho != SEMENOV_NRHO || ntemp != SEMENOV_NTEMP) {
    if (CPU_Master) fprintf(stderr,"Semenov table: %d x %d; expected: %d x %d\n",nrho,ntemp,SEMENOV_NRHO,SEMENOV_NTEMP);
    fclose(fp);
    prs_error("Unexpected Semenov opacity table dimensions.");
  }

  for (i=0;i<SEMENOV_NRHO;i++) {
    if (fscanf(fp,"%lf",&axis_value) != 1 || !isfinite(axis_value)) {
      fclose(fp);
      prs_error("Cannot read Semenov density axis.");
    }
  }

  for (j=0;j<SEMENOV_NTEMP;j++) {
    if (fscanf(fp,"%lf",&axis_value) != 1 || !isfinite(axis_value)) {
      fclose(fp);
      prs_error("Cannot read Semenov temperature axis.");
    }
  }

  SemenovLogKappa_cpu=(real*)malloc(nbytes);
  if (SemenovLogKappa_cpu == NULL) {
    fclose(fp);
    prs_error("Cannot allocate Semenov opacity table.");
  }

  for (i=0;i<SEMENOV_NRHO;i++) {
    for (j=0;j<SEMENOV_NTEMP;j++) {
      if (fscanf(fp,"%lf",&value) != 1 || !isfinite(value)) {
        fclose(fp);
        free(SemenovLogKappa_cpu);
        SemenovLogKappa_cpu=NULL;
        prs_error("Cannot read Semenov opacity value.");
      }
      SemenovLogKappa_cpu[i*SEMENOV_NTEMP+j]=(real)value;
    }
  }

  fclose(fp);

#ifdef GPU
  if (DevMalloc((void*)&SemenovLogKappa_gpu,nbytes) != 0) {
    free(SemenovLogKappa_cpu);
    SemenovLogKappa_cpu=NULL;
    prs_error("Cannot allocate Semenov opacity table on GPU.");
  }

  if (DevMemcpyH2D(SemenovLogKappa_gpu,SemenovLogKappa_cpu,nbytes) != 0) {
    free(SemenovLogKappa_cpu);
    SemenovLogKappa_cpu=NULL;
    prs_error("Cannot copy Semenov opacity table to GPU.");
  }
#endif

  SemenovOpacityLoaded=1;

  if (CPU_Master) {
    real rho_unit=(MSTAR_CGS/MSTAR)*pow(R0/R0_CGS,3.0);
    real kappa_unit=(MSTAR_CGS/MSTAR)*pow(R0/R0_CGS,2.0);
    real temp_factor=(G_CGS*MSTAR_CGS/(R0_CGS*RAD_OPACITY_R_MU_CGS))*(R0/(G*MSTAR));

    printf("\nSEMENOV_OPACITY loaded\n");
    printf("  model                 = nrm / h / s\n");
    printf("  mean                  = Rosseland\n");
    printf("  table                 = %d x %d\n",SEMENOV_NRHO,SEMENOV_NTEMP);
    printf("  mu                    = %.8f\n",(double)RAD_OPACITY_MU);
    printf("  rho unit [g cm^-3]    = %.12e\n",(double)rho_unit);
    printf("  kappa unit [code/cgs] = %.12e\n",(double)kappa_unit);
    printf("  T factor [K]          = %.12e\n",(double)temp_factor);
    printf("  file                  = %s\n\n",SEMENOV_OPACITY_FILE);
    fflush(stdout);
  }
}

void RadiationOpacityTableFree(void) {
  if (SemenovLogKappa_cpu != NULL) {
    free(SemenovLogKappa_cpu);
    SemenovLogKappa_cpu=NULL;
  }
  SemenovOpacityLoaded=0;
}

#endif

void RadiationOpacityTable_cpu(Field* Rho,Field* Energy,Field* KappaR) {

//<USER_DEFINED>
  INPUT(Rho);
  INPUT(Energy);
  OUTPUT(KappaR);
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho=Rho->field_cpu;
  real* energy=Energy->field_cpu;
  real* kappar=KappaR->field_cpu;
  real* opacity_table=SemenovLogKappa_cpu;
  real rho_unit=(MSTAR_CGS/MSTAR)*(R0/R0_CGS)*(R0/R0_CGS)*(R0/R0_CGS);
  real kappa_unit=(MSTAR_CGS/MSTAR)*(R0/R0_CGS)*(R0/R0_CGS);
  real temp_factor=(G_CGS*MSTAR_CGS/(R0_CGS*RAD_OPACITY_R_MU_CGS))*(R0/(G*MSTAR));
  int pitch=Pitch_cpu;
  int stride=Stride_cpu;
  int size_x=Nx+2*NGHX;
  int size_y=Ny+2*NGHY;
  int size_z=Nz+2*NGHZ;
//<\EXTERNAL>

//<CONSTANT>
// real GAMMA(1);
//<\CONSTANT>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ir;
  int it;
  int idx00;
  int idx10;
  int idx01;
  int idx11;
  real rho_code;
  real e_code;
  real rho_cgs;
  real temperature;
  real rho_eval;
  real temp_eval;
  real logrho;
  real logtemp;
  real xr;
  real xt;
  real fr;
  real ft;
  real z00;
  real z10;
  real z01;
  real z11;
  real logkappa;
  real kappa_cgs;
  real kappa_code;
//<\INTERNAL>

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
        rho_code=rho[l];
        e_code=energy[l];

        if (!isfinite(rho_code) || !isfinite(e_code) || rho_code <= 0.0 || e_code <= 0.0) {
          kappar[l]=0.0;
        } else {
          rho_cgs=rho_code*rho_unit;
          temperature=(GAMMA-1.0)*(e_code/rho_code)*temp_factor;

          rho_eval=rho_cgs;
          if (rho_eval < 1.0e-18) rho_eval=1.0e-18;
          if (rho_eval > 1.0e-7) rho_eval=1.0e-7;

          temp_eval=temperature;
          if (temp_eval < SEMENOV_T_MIN) temp_eval=SEMENOV_T_MIN;
          if (temp_eval > SEMENOV_T_MAX) temp_eval=SEMENOV_T_MAX;

          logrho=log10(rho_eval);
          logtemp=log10(temp_eval);

          xr=(logrho-SEMENOV_LOGRHO_MIN)*SEMENOV_INV_DLOGRHO;
          xt=(logtemp-SEMENOV_LOGT_MIN)*SEMENOV_INV_DLOGT;

          ir=(int)floor(xr);
          it=(int)floor(xt);

          if (ir < 0) ir=0;
          if (it < 0) it=0;
          if (ir > SEMENOV_NRHO-2) ir=SEMENOV_NRHO-2;
          if (it > SEMENOV_NTEMP-2) it=SEMENOV_NTEMP-2;

          fr=xr-(real)ir;
          ft=xt-(real)it;

          if (fr < 0.0) fr=0.0;
          if (fr > 1.0) fr=1.0;
          if (ft < 0.0) ft=0.0;
          if (ft > 1.0) ft=1.0;

          idx00=ir*SEMENOV_NTEMP+it;
          idx10=(ir+1)*SEMENOV_NTEMP+it;
          idx01=ir*SEMENOV_NTEMP+it+1;
          idx11=(ir+1)*SEMENOV_NTEMP+it+1;

          z00=opacity_table[idx00];
          z10=opacity_table[idx10];
          z01=opacity_table[idx01];
          z11=opacity_table[idx11];

          logkappa=(1.0-fr)*(1.0-ft)*z00
                  +fr*(1.0-ft)*z10
                  +(1.0-fr)*ft*z01
                  +fr*ft*z11;

          kappa_cgs=pow(10.0,logkappa);
          kappa_code=kappa_cgs*kappa_unit;

          if (!isfinite(kappa_code) || kappa_code <= 0.0) kappar[l]=0.0;
          else kappar[l]=kappa_code;
        }
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