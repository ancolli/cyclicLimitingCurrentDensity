volScalarField DEff("DEff",turbulence->nu()/Sc + turbulence->nut()/Sct);

fvScalarMatrix CEqn
(
//fvm::ddt(C)
  fvm::div(phi, C_ox)
- fvm::laplacian(DEff, C_ox)
//+ (flowDir & U)*gamma
);

CEqn.relax();

fvOptions.constrain(CEqn);

CEqn.solve();

