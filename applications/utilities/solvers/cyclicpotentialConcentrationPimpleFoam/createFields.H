#include "createRDeltaT.H"

// * * * * * * * * * * * * * Properties fields * * * * * * * * * * * * * * * //

Info<< "Reading speciesProperties\n" << endl;

IOdictionary speciesProperties
(
    IOobject
    (
        "speciesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field C_ox\n" << endl;
volScalarField C_ox
(
    IOobject
    (
        "C_ox",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//From here

Info<< "Reading field fi\n" << endl;

volScalarField fi
(
    IOobject
    (
	    "fi",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Flux
(
	IOobject
	(
	    "Flux",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE//AUTO_WRITE//
	),
	mesh,
	dimensionedScalar("0", dimensionSet(0,-2,0,0,0,1,0), 0.0)
);

// * * * * * * * * * * * * * Reading Properties * * * * * * * * * * * * * * * * //

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
	IOobject
	(
	    "transportProperties",
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ_IF_MODIFIED,
	    IOobject::NO_WRITE
	)
);

#include "createPhi.H"
#include "readTransportProperties.H"

/*volScalarField keff
(
	IOobject
	(
	    "keff",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("keff", kf)
);*/

Info<< "Reading controlProperties\n" << endl;

IOdictionary controlProperties
(
    IOobject
    (
        "controlProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

/*IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    )
);*/


IOdictionary Results
(
    IOobject
    (
        "Results",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    )
);

/*volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh
    ),
    C_ox - Cinlet  
);*/

volScalarField cBar
(
    IOobject
    (
        "cBar",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
