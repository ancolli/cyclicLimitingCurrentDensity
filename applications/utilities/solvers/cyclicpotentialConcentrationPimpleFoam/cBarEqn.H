fvScalarMatrix cBarEqn
(
	fvm::ddt(cBar)
	+ fvm::div(phi, cBar)
	==
	fvm::laplacian(DEff, cBar)
	+ fvm::SuSp(((DEff*pow(lambda,2)) + (lambda*(flowDir & (U - fvc::grad(DEff))))),cBar) //- (2*DEff*lambda*(flowDir & fvc::grad(cBar)))
);

cBarEqn.relax();

fvOptions.constrain(cBarEqn);

cBarEqn.solve();

cBar = cmax*cBar/gMax(cBar.boundaryField()[inletID]);

/*forAll(theta, celli)
{
    theta[celli] = Foam::exp(-lambda.value()*(flowDir & mesh.C()[celli]))*thetaBar[celli];
}*/

C_ox = Foam::exp(-lambda*(flowDir & mesh.C()))*cBar;//theta + Cinlet;
