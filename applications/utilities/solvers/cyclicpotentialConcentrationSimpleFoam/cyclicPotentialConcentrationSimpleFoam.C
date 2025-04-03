/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flows.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting pressure/momentum loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"

            if(method == "wall")
            {   
                //Info << "wall" << endl;
                #include "lambdaEqn.H"
                #include "cBarEqn.H"
                //#include "thetaEqn.H"
                
                //#include "fiEqn.H"
            }
            else if(method == "Real")
            {
                Info << "Real" << endl;
                #include "CEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    const boundBox& globalBb = mesh.bounds();
    scalar Volume = mag((globalBb.max()[1]-globalBb.min()[1])*(globalBb.max()[2]-globalBb.min()[2])*(globalBb.max()[3]-globalBb.min()[3]));
    
    reduce(Volume, sumOp<scalar>());
     
    Info << "Total Volume: " << Volume << " m^3 "<< endl;
    
    //calculate mesh volume    
    scalar sumVf = 0.0;
    forAll(mesh.V(), cellI)
    {
	    {
		sumVf += mesh.V()[cellI];
	    }
    }

    // Sync sum across processors
    reduce(sumVf, sumOp<scalar>());

    Info<< "Partial Volume: " << sumVf << " m^3 "<< endl;
    
    //porosity
    
    Info<< "Porosity: " << (1.0-sumVf/Volume) <<  endl;

    
    // Calculates the specific surface area
        
    /*label patchi = mesh.boundaryMesh().findPatchID("wall_mesh");

    scalar area = gSum(mesh.magSf().boundaryField()[patchi]);
    
    // Sync sum across processors
    reduce(area, sumOp<scalar>());
    
    Info<< "Total Area = " << area << endl;
    
    scalar A_s  = area/Volume;

    Info<< "Specific surface area = " << A_s << " m^-1\n " << endl;    */
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
