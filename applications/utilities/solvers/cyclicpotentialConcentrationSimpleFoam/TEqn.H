volScalarField alphaEff("alphaEff",turbulence->nu()/Pr + turbulence->nut()/Prt);

fvScalarMatrix TEqn
(
//fvm::ddt(C)
  fvm::div(phi, T)
- fvm::laplacian(alphaEff, T)
//+ (flowDir & U)*gamma
);

TEqn.solve();

