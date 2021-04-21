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
    Foam::relativeVelocityModel

Description

SourceFiles
    relativeVelocityModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativeVelocityModel_H
#define relativeVelocityModel_H

#include "fvCFD.H"
#include "dictionary.H"
#include "incompressibleNinePhaseInteractingMixture.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relativeVelocityModel Declaration
\*---------------------------------------------------------------------------*/

class relativeVelocityModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        relativeVelocityModel(const relativeVelocityModel&);

        //- Disallow default bitwise assignment
        void operator=(const relativeVelocityModel&);

protected:

    // Protected data

        //- Mixture properties
        const incompressibleNinePhaseInteractingMixture& mixture_;

        //- Name of the continuous phase
        const word continuousPhaseName_;
        
        const volVectorField& U_;
        //- Continuous phase fraction
        const volScalarField& alphac_;

        //- Dispersed phase fraction
        const volScalarField& alphad1_;
        const volScalarField& alphad2_;
        const volScalarField& alphad3_;
        const volScalarField& alphad4_;
        const volScalarField& alphad5_;
        const volScalarField& alphad6_;
        const volScalarField& alphad7_;
        const volScalarField& alphad8_;
        //- Continuous density
        const dimensionedScalar& rhoc_;

        //- Dispersed density
        const dimensionedScalar& rhod1_;
        const dimensionedScalar& rhod2_;
        const dimensionedScalar& rhod3_;
        const dimensionedScalar& rhod4_;
        const dimensionedScalar& rhod5_;
        const dimensionedScalar& rhod6_;
        const dimensionedScalar& rhod7_;
        const dimensionedScalar& rhod8_;
        
        const volScalarField& mu_; 
        
        const dimensionedScalar& dd1_;
        const dimensionedScalar& dd2_;
        const dimensionedScalar& dd3_;
        const dimensionedScalar& dd4_;
        const dimensionedScalar& dd5_;
        const dimensionedScalar& dd6_;
        const dimensionedScalar& dd7_;
 //       const dimensionedScalar& dd8_;
               
        const dimensionedVector grav_;
        const dimensionedScalar DeltaT_;
      
        mutable volVectorField MixForce_;
        mutable volVectorField LiftForce_;
        mutable volVectorField LiftForce1_;
        mutable volVectorField LiftForce2_;
        mutable volVectorField LiftForce3_;
        mutable volVectorField LiftForce4_;
        mutable volVectorField LiftForce5_;
        mutable volVectorField LiftForce6_;
       
        
        mutable volScalarField Rep1_;
        mutable volScalarField Rep2_;
        mutable volScalarField Rep3_;
        mutable volScalarField Rep4_;
        
        mutable volScalarField Rep5_;
        mutable volScalarField Rep6_;
        mutable volScalarField Rep7_;
        mutable volScalarField Rep8_;
        
        //- Dispersed diffusion velocity
        mutable volVectorField Udm1_;
        mutable volVectorField Udm2_;
        mutable volVectorField Udm3_;
        mutable volVectorField Udm4_;
         mutable volVectorField Udm5_;
        mutable volVectorField Udm6_;
        mutable volVectorField Udm7_;
        mutable volVectorField Udm8_;


public:

    //- Runtime type information
    TypeName("relativeVelocityModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relativeVelocityModel,
        dictionary,
        (const dictionary& dict,
         const incompressibleNinePhaseInteractingMixture& mixture),
        (dict, mixture)
    );


    // Constructors

        //- Construct from components
        relativeVelocityModel
        (
            const dictionary& dict,
            const incompressibleNinePhaseInteractingMixture& mixture
        );


    // Selector
    static autoPtr<relativeVelocityModel> New
    (
        const dictionary& dict,
        const incompressibleNinePhaseInteractingMixture& mixture
    );


    //- Destructor
    virtual ~relativeVelocityModel();


    // Member Functions

        //- Mixture properties
        const incompressibleNinePhaseInteractingMixture& mixture() const
        {
            return mixture_;
        }

        //- Return the mixture mean density
        tmp<volScalarField> rho() const;
        tmp<volVectorField> PrevU() const;
        tmp<volScalarField> mu() const;
        tmp<volScalarField> Rep1() const;
        tmp<volScalarField> Rep2() const;
        tmp<volScalarField> Rep3() const;
        tmp<volScalarField> Rep4() const;
         tmp<volScalarField> Rep5() const;
        tmp<volScalarField> Rep6() const;
  //      tmp<volScalarField> Rep7() const;
  //      tmp<volScalarField> Rep8() const;
        tmp<volVectorField> Force() const;
        tmp<volVectorField> LiftForce1() const;
        tmp<volVectorField> LiftForce2() const;
        tmp<volVectorField> LiftForce3() const;
        tmp<volVectorField> LiftForce4() const;
         tmp<volVectorField> LiftForce5() const;
        tmp<volVectorField> LiftForce6() const;
        tmp<volVectorField> LiftForce7() const;
        tmp<volVectorField> LiftForce8() const;
        //- Return the diffusion velocity of the dispersed phase
        
        const volVectorField& Udm1() const
        {
            return Udm1_;
        }
        
        const volVectorField& Udm2() const
        {
            return Udm2_;
        }
        
        const volVectorField& Udm3() const
        {
            return Udm3_;
        }
        
        const volVectorField& Udm4() const
        {
            return Udm4_;
        }
          const volVectorField& Udm5() const
        {
            return Udm5_;
        }
          const volVectorField& Udm6() const
        {
            return Udm6_;
        }
          const volVectorField& Udm7() const
        {
            return Udm7_;
        }
          const volVectorField& Udm8() const
        {
            return Udm8_;
        }

        //- Return the stress tensor due to the phase transport
        tmp<volSymmTensorField> tauDm() const;

        //- Update the diffusion velocity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

using Foam::relativeVelocityModel;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
