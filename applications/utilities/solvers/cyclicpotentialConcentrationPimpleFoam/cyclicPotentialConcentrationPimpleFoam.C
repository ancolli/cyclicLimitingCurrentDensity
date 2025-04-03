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
    pimpleFoam

Group
    grpIncompressibleSolvers

Description
     Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();
    
    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // needed to keep track of how many time steps have happened
    int iterations = 0.0;
    
    
    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // update iterations
        iterations++;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }           
                      
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
                      
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

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
