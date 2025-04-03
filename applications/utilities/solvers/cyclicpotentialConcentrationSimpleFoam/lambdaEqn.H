// Calculate decay constant lambda

volScalarField DEff("DEff",turbulence->nu()/Sc + turbulence->nut()/Sct);

label inletID = mesh.boundaryMesh().findPatchID("inlet");
label wallsID = mesh.boundaryMesh().findPatchID("wall_cathode");

if(runTime.value() >= delay.value())
{
    //if(int(runTime.value()) % int(interval.value()) == 0)
    {
        
        // Calculate decay constant lambda
        dimensionedScalar lambdaNum
        (
            "lambdaNum",
            dimensionSet(0, 0, -1, 0, 1, 0, 0),
            gSum(DEff.boundaryField()[wallsID]*-C_ox.boundaryField()[wallsID].snGrad()*mesh.boundary()[wallsID].magSf())
        );

        dimensionedScalar lambdaDen
        (
            "lambdaDen",
            dimensionSet(0, 0, -1, 0, 1, 0, 0),
            gSum(cBar.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf()) +  //
            gSum(DEff.boundaryField()[inletID]*(flowDir & fvc::grad(cBar)().boundaryField()[inletID])*mesh.boundary()[inletID].magSf())
        );

        //Info << "Ux*Theta: " << gSum((thetaBar.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf())) << endl;
        //Info << "alpha*dTheta/dx: " << gSum(alphaEff.boundaryField()[inletID]*(flowDir & fvc::grad(thetaBar)().boundaryField()[inletID])*mesh.boundary()[inletID].magSf()) << endl;

        //Info << "lambda_num: " << lambdaNum.value() << endl;
        //Info << "lambda_den: " << lambdaDen.value() << endl;
        //Info << "lambda_num/den: " << lambdaNum.value()/lambdaDen.value() << endl;
	//dimensionedScalar lambda0 = lambda;
        lambda = -(1/Lc)*Foam::log(1 - (lambdaNum/lambdaDen));
        //lambda = lambda0 + 0.2*(lambda-lambda0);
    }
}
else
{
   lambda = lambda0;
}

Info << "Lambda: " << lambda.value() << endl;
