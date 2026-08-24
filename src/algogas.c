#include "fargo3d.h"

TimeProcess t_Comm;
TimeProcess t_Hydro;
TimeProcess t_Mhd;
TimeProcess t_sub1;
TimeProcess t_sub1_x;
TimeProcess t_sub1_y;
TimeProcess t_sub1_z;


void FillGhosts(int var) {

  InitSpecificTime(&t_Comm,"MPI Communications");
  FARGO_SAFE(comm(var));
  GiveSpecificTime(t_Comm);


  FARGO_SAFE(boundaries());


#if defined(Y)
  if (NY==1) CheckMuteY();
#endif
#if defined(Z)
  if (NZ==1) CheckMuteZ();
#endif
}

void Sources(real dt) {
  SetupHook1();

#ifdef ADIABATIC
  FARGO_SAFE(ComputePressureFieldAd());
#endif
#ifdef ISOTHERMAL
  FARGO_SAFE(ComputePressureFieldIso());
#endif
#ifdef POLYTROPIC
  FARGO_SAFE(ComputePressureFieldPoly());
#endif

  InitSpecificTime(&t_Hydro,"Eulerian Hydro (no transport) algorithms");

#ifdef POTENTIAL
  FARGO_SAFE(compute_potential(dt));
  if (Corotating) FARGO_SAFE(CorrectVtheta(Domega));
#endif

#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC(Vy));
#endif

#ifdef X
  FARGO_SAFE(SubStep1_x(dt));
#endif
#ifdef Y
  FARGO_SAFE(SubStep1_y(dt));
#endif
#ifdef Z
  FARGO_SAFE(SubStep1_z(dt));
#endif

#if (defined(VISCOSITY) || defined(ALPHAVISCOSITY))
  if (Fluidtype==GAS) viscosity(dt);
#endif

#ifndef NOSUBSTEP2
  FARGO_SAFE(SubStep2_a(dt));
  FARGO_SAFE(SubStep2_b(dt));
#endif

#ifdef ADIABATIC
  if (Fluidtype==GAS) {
#ifdef DISKCOOLING
    FARGO_SAFE(disk_cooling_opacity_bl(dt));
#endif
#ifdef STELLARHEATING
    FARGO_SAFE(BuildStellarOpticalDepth());
    FARGO_SAFE(stellar_heating_term(dt));
#endif
#ifdef VISCOUSHEATING
    FARGO_SAFE(Reset_field(Viscous_heat_pow));
    FARGO_SAFE(viscous_heating_term(dt));
#endif
    FARGO_SAFE(SubStep3(dt));
  }
#endif

  GiveSpecificTime(t_Hydro);

#ifdef MHD
  if (Fluidtype==GAS) {
    InitSpecificTime(&t_Mhd,"MHD algorithms");
    FARGO_SAFE(copy_velocities(VTEMP2V));
#ifndef STANDARD
    FARGO_SAFE(ComputeVmed(Vx));
    FARGO_SAFE(ChangeFrame(-1,Vx,VxMed));
    VxIsResidual=YES;
#endif

    ComputeMHD(dt);

#ifndef STANDARD
    FARGO_SAFE(ChangeFrame(+1,Vx,VxMed));
    VxIsResidual=NO;
#endif
    FARGO_SAFE(copy_velocities(V2VTEMP));
    GiveSpecificTime(t_Mhd);
  }
#endif

  InitSpecificTime(&t_Hydro,"Transport algorithms");

#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC(Vy_temp));
#endif

  FARGO_SAFE(copy_velocities(VTEMP2V));
  FARGO_SAFE(FillGhosts(PrimitiveVariables()));
  FARGO_SAFE(copy_velocities(V2VTEMP));

#ifdef MHD
  if (Fluidtype==GAS) {
    FARGO_SAFE(UpdateMagneticField(dt,1,0,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,1,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,0,1));

#if !defined(STANDARD)
    FARGO_SAFE(MHD_fargo(dt));
#endif
  }
#endif
}

void Transport(real dt) {
#ifdef X
#ifndef STANDARD
  FARGO_SAFE(ComputeVmed(Vx_temp));
#endif
#endif

  transport(dt);

  GiveSpecificTime(t_Hydro);

  if (ForwardOneStep==YES) prs_exit(EXIT_SUCCESS);

#ifdef MHD
  if (Fluidtype==GAS) {
    *(Emfx->owner)=Emfx;
    *(Emfy->owner)=Emfy;
    *(Emfz->owner)=Emfz;
  }
#endif
}