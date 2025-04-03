volScalarField alphaEff("alphaEff",turbulence->nu()/Pr + turbulence->nut()/Prt);

fvScalarMatrix TEqn
(
  fvm::ddt(T)
+ fvm::div(phi, T)
- fvm::laplacian(alphaEff, T)
//+ (flowDir & U)*gamma
);

TEqn.solve();

