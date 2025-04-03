/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::functionObjects::speciesWallFlux

Group

Description
    This function object...

    Example of function object specification:
    \verbatim
    speciesWallFlux1
    {
        type           speciesWallFlux;
        libs ("libspeciesWallFluxFunctionObject.so");
        ...
        region      fluid;
        patches     (".*Wall");
	CName	    speciesName;
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: speciesWallFlux        | yes        |
        patches  | list of patches to process | no         | all wall patches
        region   | region to be evaluated     | no         | default region
    \endtable

SourceFiles
    speciesWallFluxMean.C

\*---------------------------------------------------------------------------*/

#ifndef speciesWallFluxMean_H
#define speciesWallFluxMean_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

//#include "volFields.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class speciesWallFluxMean Declaration
\*---------------------------------------------------------------------------*/

class speciesWallFluxMean
:
    public fvMeshFunctionObject,
    public writeFile
{


protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;

	//- Faraday's constant
	dimensionedScalar F_;

	//- Diffusion coefficient
        dimensionedScalar D_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the flux
        void calcFlux
        (
            //const volScalarField& D_red_,
            const volScalarField& C_, //const
            volScalarField& speciesWallFluxMean
        );

private:

//Private data

//- Name of flux field
        word CName_;
	
        //volScalarField C_;

// Private Member Functions

        //- Disallow default bitwise copy construct
        speciesWallFluxMean(const speciesWallFluxMean&);

        //- Disallow default bitwise assignment
        void operator=(const speciesWallFluxMean&);


public:



    //- Runtime type information
    TypeName("speciesWallFluxMean");


    // Constructors

        //- Construct from Time and dictionary
        speciesWallFluxMean
        (
            const word& name,
            const Time& runTime,	    
            const dictionary& dict
        );


    //- Destructor
    virtual ~speciesWallFluxMean();


    // Member Functions

        //- Read the speciesWallFluxMean data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the speciesWallFluxMean
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
