/** \file main.c

Main file of the distribution.
Manages the call to initialization
functions, then the main loop.

*/

#include "fargo3d.h"

int begin_i=0,NbRestart=0;
int InnerOutputCounter=0,StillWriteOneOutput;
real dt;
real dtemp=0.0;


int main(int argc,char *argv[]) {
  int i=0,OutputNumber=0,d;
  char sepline[]="===========================";

  sprintf(FirstCommand,"%s",argv[0]);
  sprintf(CommandLine,"%s ",argv[0]);
  while (++i<argc) {
    strncat(CommandLine,argv[i],1023);
    strncat(CommandLine," ",1023);
  }

#ifdef LONGSUMMARY
  strncpy(StickyOptions,ExtractFromExecutable(YES,"",1),1023);
  strncpy(BoundaryFile,ExtractFromExecutable(YES,"",3),4095);
#endif

  strcpy(ParameterFile,"");
  for (i=1;i<argc;i+=d) {
    d=1;

    if (*(argv[i])=='+') {
      if (strspn(argv[i],"+S#D")!=strlen(argv[i])) PrintUsage(argv[0]);

      if (strchr(argv[i],'#')) {
        d=2;
        ArrayNb=atoi(argv[i+1]);
        EarlyOutputRename=YES;
        if (ArrayNb<=0) {
          masterprint("Incorrect Array number after +# flag\n");
          PrintUsage(argv[0]);
        }
      }

      if (strchr(argv[i],'D')) {
        d=2;
        strcpy(DeviceFile,argv[i+1]);
        DeviceFileSpecified=YES;
      }

      if (strchr(argv[i],'S')) {
        d=2;
        StretchNumber=atoi(argv[i+1]);
        StretchOldOutput=YES;
      }
    }

    if (*(argv[i])=='-') {
      if (strspn(argv[i],"-tofCmkspSVBD0#")!=strlen(argv[i])) PrintUsage(argv[0]);

      if (strchr(argv[i],'t')) TimeInfo=YES;
      if (strchr(argv[i],'f')) ForwardOneStep=YES;
      if (strchr(argv[i],'0')) OnlyInit=YES;

      if (strchr(argv[i],'C')) {
        EverythingOnCPU=YES;
#ifdef GPU
        mastererr("WARNING: Forcing execution of all functions on CPU\n");
#else
        mastererr("WARNING: Flag -C meaningless for a CPU built\n");
#endif
      }

      if (strchr(argv[i],'m')) Merge=YES;
      if (strchr(argv[i],'k')) Merge=NO;

      if (strchr(argv[i],'o')) {
        RedefineOptions=YES;
        ParseRedefinedOptions(argv[i+1]);
        d=2;
      }

      if (strchr(argv[i],'s')) {
        Restart=YES;
        d=2;
        NbRestart=atoi(argv[i+1]);
        if (NbRestart<0) {
          masterprint("Incorrect restart number\n");
          PrintUsage(argv[0]);
        }
      }

      if (strchr(argv[i],'#')) {
        d=2;
        ArrayNb=atoi(argv[i+1]);
        if (ArrayNb<=0) {
          masterprint("Incorrect Array number after -# flag\n");
          PrintUsage(argv[0]);
        }
      }

      if (strchr(argv[i],'p')) PostRestart=YES;

      if (strchr(argv[i],'S')) {
        Restart_Full=YES;
        d=2;
        NbRestart=atoi(argv[i+1]);
        if (NbRestart<0) {
          masterprint("Incorrect restart number\n");
          PrintUsage(argv[0]);
        }
      }

      if (strchr(argv[i],'V')) {
        Dat2vtk=YES;
        Restart_Full=YES;
        d=2;
        NbRestart=atoi(argv[i+1]);
        if (NbRestart<0) {
          masterprint("Incorrect output number\n");
          PrintUsage(argv[0]);
        }
      }

      if (strchr(argv[i],'B')) {
        Vtk2dat=YES;
        Restart_Full=YES;
        d=2;
        NbRestart=atoi(argv[i+1]);
        if (NbRestart<0) {
          masterprint("Incorrect output number\n");
          PrintUsage(argv[0]);
        }
      }

      if (strchr(argv[i],'D')) {
        d=2;
        DeviceManualSelection=atoi(argv[i+1]);
      }
    } else strcpy(ParameterFile,argv[i]);
  }

#ifdef WRITEGHOSTS
  if (Merge==YES) {
    mastererr("Cannot merge outputs when dumping ghost values.\n");
    mastererr("'make nofulldebug' could fix this problem.\n");
    mastererr("Using the -k flag could be another solution.\n");
    prs_exit(1);
  }
#endif

#ifdef MPICUDA
  EarlyDeviceSelection();
#endif

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&CPU_Rank);
  MPI_Comm_size(MPI_COMM_WORLD,&CPU_Number);
  CPU_Master=(CPU_Rank==0 ? 1 : 0);

  if (strlen(xstr(VERSION))<2)
    sprintf(VersionString,"FARGO3D Public version 1.3");
  else
    sprintf(VersionString,"FARGO3D git version %s",xstr(VERSION));

  masterprint("\n\n%s\n%s\nSETUP: '%s'\n%s\n\n",
              sepline,VersionString,xstr(SETUPNAME),sepline);

  if ((ParameterFile[0]==0) || (argc==1)) PrintUsage(argv[0]);

#ifndef MPICUDA
  SelectDevice(CPU_Rank);
#endif

  InitVariables();
  MPI_Barrier(MPI_COMM_WORLD);
  ReadDefaultOut();
  ReadVarFile(ParameterFile);

  if (strcmp(PLANETCONFIG,"NONE")!=0) ThereArePlanets=YES;

  if (ORBITALRADIUS>1.0e-30) {
    YMIN*=ORBITALRADIUS;
    YMAX*=ORBITALRADIUS;
    DT*=sqrt(ORBITALRADIUS*ORBITALRADIUS*ORBITALRADIUS);
  }

  SubsDef(OUTPUTDIR,DefaultOut);

  if ((ArrayNb) && (EarlyOutputRename==YES)) {
    i=strlen(OUTPUTDIR);
    if (OUTPUTDIR[i-1]=='/') OUTPUTDIR[i-1]=0;
    sprintf(OUTPUTDIR,"%s%06d/",OUTPUTDIR,ArrayNb);
    printf("\n\n***\n\nNew Output Directory is %s\n\n***\n\n",OUTPUTDIR);
    MakeDir(OUTPUTDIR);
  }

  MakeDir(OUTPUTDIR);

#if !defined(X)
  NX=1;
#endif
#if !defined(Y)
  NY=1;
#endif
#if !defined(Z)
  NZ=1;
#endif

  SelectWriteMethod();

#if !defined(Y) && !defined(Z)
  if (CPU_Rank==1)
    prs_error("You cannot split a 1D mesh in x. Sequential runs only!");

  if (CPU_Number>1) {
    MPI_Finalize();
    prs_exit(EXIT_FAILURE);
  }
#endif

  ListVariables("variables.par");
  ListVariablesIDL("IDL.var");

  ChangeArch();
  split(&Gridd);
  InitSpace();
  WriteDim();
  InitSurfaces();
  LightGlobalDev();
  CreateFields();

  Sys=InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();

  if (Corotating) OMEGAFRAME=GetPsysInfo(FREQUENCY);
  OMEGAFRAME0=OMEGAFRAME;

  if (Restart==YES || Restart_Full==YES) {
    CondInit();


    MULTIFLUID(begin_i=RestartSimulation(NbRestart));


    if (ThereArePlanets) {
      PhysicalTime=GetfromPlanetFile(NbRestart,9,0);
      OMEGAFRAME=GetfromPlanetFile(NbRestart,10,0);
      RestartPlanetarySystem(NbRestart,Sys);
    }
  } else {
    if (ThereArePlanets) EmptyPlanetSystemFiles();
    CondInit();

  }

  if (StretchOldOutput==YES) StretchOutput(StretchNumber);


  MULTIFLUID(comm(ENERGY));


  if ((ArrayNb) && (EarlyOutputRename==NO)) {
    i=strlen(OUTPUTDIR);
    if (OUTPUTDIR[i-1]=='/') OUTPUTDIR[i-1]=0;

    sprintf(OUTPUTDIR,"%s%06d/",OUTPUTDIR,ArrayNb);
    printf("\n\n***\n\nNew Output Directory is %s\n\n***\n\n",OUTPUTDIR);

    MakeDir(OUTPUTDIR);
    ListVariables("variables.par");
    ListVariablesIDL("IDL.var");
    InitSpace();
    WriteDim();
  }

#if !defined(STANDARD) && !defined(RAM)
  if (*SPACING=='N') {
    mastererr("Warning: Mesh must be unifrom in x-direction to use FARGO advection.\n");
    mastererr("Warning: If SPACING is set to 'N', use -DSTANDARD or -DRAM in .opt file.\n");
    prs_exit(1);
  }
#endif

#if defined(MHD) && defined(RAM)
  mastererr("Warning: Current version of MHD module does not include RAM.\n");
  mastererr("Warning: If SPACING is set to 'N', use -DSTANDARD in .opt file.\n");
  prs_exit(1);
#endif

  GetHostsList();
  DumpToFargo3drc(argc,argv);

  SelectArchFileName();

#ifdef LONGSUMMARY
  ExtractFromExecutable(NO,ArchFile,2);
#endif


  MULTIFLUID(FillGhosts(PrimitiveVariables()));


#ifdef STOCKHOLM
  FARGO_SAFE(init_stockholm());
#endif

#ifdef GHOSTSX
  masterprint("\n\nNew version with ghost zones in X activated\n");
#else
  masterprint("Standard version with no ghost zones in X\n");
#endif

  for (i=begin_i;i<=NTOT;i++) {
    if (NINTERM*(TimeStep=(i/NINTERM))==i) {

#if defined(MHD) && defined(DEBUG)
      FARGO_SAFE(ComputeDivergence(Bx,By,Bz));
#endif

      if (ThereArePlanets)
        WritePlanetSystemFile(TimeStep,NO);

#ifndef NOOUTPUTS

      MULTIFLUID(WriteOutputs(ALL));


#ifdef MATPLOTLIB
      Display();
#endif

      if (CPU_Master)
        printf("OUTPUTS %d at date t = %f OK\n",TimeStep,PhysicalTime);
#endif

      if (TimeInfo==YES) GiveTimeInfo(TimeStep);
    }

    if (NSNAP!=0) {
      if (NSNAP*(TimeStep=(i/NSNAP))==i) {
        MULTIFLUID(WriteOutputs(SPECIFIC));

#ifdef MATPLOTLIB
        Display();
#endif
      }
    }

    if (i==NTOT) break;

    dtemp=0.0;

    while (dtemp<DT) {


#ifdef X
#ifndef STANDARD
      MULTIFLUID(ComputeVmed(Vx));
#endif
#endif


#ifdef FLOOR
      MULTIFLUID(Floor());
#endif


#ifdef MHD
#ifdef OHMICDIFFUSION
      FARGO_SAFE(OhmicDiffusion_coeff());
#endif
#ifdef AMBIPOLARDIFFUSION
      FARGO_SAFE(AmbipolarDiffusion_coeff());
#endif
#ifdef HALLEFFECT
      FARGO_SAFE(HallEffect_coeff());
#endif
#endif


      MULTIFLUID(cfl());


      CflFluidsMin();
      dt=StepTime;

      dtemp+=dt;
      if (dtemp>DT) dt=DT-(dtemp-dt);

#ifdef RADDIFFTEST

#ifdef RADIATION
      FARGO_SAFE(RadiationDiffusion(dt));
#endif

#else

      FARGO_SAFE(Reset_field(Total_Density));
      MULTIFLUID(ComputeTotalDensity());

#ifdef COLLISIONPREDICTOR
      FARGO_SAFE(Collisions(0.5*dt,0));
#endif

      MULTIFLUID(Sources(dt));

#ifdef DRAGFORCE
      FARGO_SAFE(Collisions(dt,1));
#endif

#ifdef DUSTDIFFUSION
      FARGO_SAFE(DustDiffusion_Main(dt));
#endif

      MULTIFLUID(Transport(dt));

#ifdef RADIATION
      FARGO_SAFE(RadiationDiffusion(dt));
#endif

#endif

      PhysicalTime+=dt;
      Timestepcount++;

#ifdef STOCKHOLM

      MULTIFLUID(StockholmBoundary(dt));

#endif


      MULTIFLUID(FillGhosts(PrimitiveVariables()));


      if (CPU_Master) {
        if (FullArrayComms)
          printf("%s","!");
        else {
          if (ContourComms)
            printf("%s",":");
          else
            printf("%s",".");
        }

#ifndef NOFLUSH
        fflush(stdout);
#endif
      }

      FullArrayComms=0;
      ContourComms=0;
    }

    if (CPU_Master) printf("%s","\n");

    MULTIFLUID(MonitorGlobal(MONITOR2D     |
                             MONITORY      |
                             MONITORY_RAW  |
                             MONITORSCALAR |
                             MONITORZ      |
                             MONITORZ_RAW));

    if (ThereArePlanets) {
      WritePlanetSystemFile(TimeStep,YES);
      SolveOrbits(Sys);
    }
  }

  MPI_Finalize();

  masterprint("End of the simulation!\n");
  return 0;
}