/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::general

Description
    General relative velocity model

SourceFiles
    general.C

\*---------------------------------------------------------------------------*/

#ifndef Maninnen_H
#define Maninnen_H

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{

/*---------------------------------------------------------------------------*\
                         Class general Declaration
\*---------------------------------------------------------------------------*/

class Maninnen
:
    public relativeVelocityModel
{
    // Private data

        //- a coefficient
        dimensionedScalar a_;

        //- a1 coefficient
        dimensionedScalar a1_;

        //- Drift velocity
        dimensionedVector g0_;
        
         dimensionedScalar dt_;

        //- Residual phase fraction
     //   dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("Maninnen");


    // Constructors

        //- Construct from components
        Maninnen
        (
            const dictionary& dict,
            const incompressibleNinePhaseInteractingMixture& mixture
        );


    //- Destructor
    ~Maninnen();


    // Member Functions

        //- Update the diffusion velocity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativeVelocityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
