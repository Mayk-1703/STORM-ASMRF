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
    Foam::incompressibleTwoPhaseInteractingMixture

Description
    A two-phase incompressible transportModel for interacting phases
    requiring the direct evaluation of the mixture viscosity,
    e.g. activated sludge or slurry.

SourceFiles
    incompressibleTwoPhaseInteractingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleNinePhaseInteractingMixture_H
#define incompressibleNinePhaseInteractingMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "mixtureViscosityModel.H"
#include "ninePhaseMixture.H"
#include "../ninePhaseInterfaceProperties/ninePhaseInterfaceProperties.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class incompressibleTwoPhaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleNinePhaseInteractingMixture
:
    public IOdictionary,
    public transportModel,
    public ninePhaseMixture,
    public ninePhaseInterfaceProperties
{
protected:

    // Protected data

        autoPtr<mixtureViscosityModel> muModel1_;
        autoPtr<mixtureViscosityModel> muModel2_;
        autoPtr<mixtureViscosityModel> muModel3_;
        autoPtr<mixtureViscosityModel> muModel4_;
        autoPtr<mixtureViscosityModel> muModel5_;
        autoPtr<mixtureViscosityModel> muModel6_;
        autoPtr<mixtureViscosityModel> muModel7_;
        autoPtr<viscosityModel> muModel8_;
        autoPtr<viscosityModel> nucModel_;

        dimensionedScalar rhod1_;
        dimensionedScalar rhod2_;
        dimensionedScalar rhod3_;
        dimensionedScalar rhod4_;
        dimensionedScalar rhod5_;
        dimensionedScalar rhod6_;
        dimensionedScalar rhod7_;
        dimensionedScalar rhod8_;
        dimensionedScalar rhoc_;

        dimensionedScalar nuc_;
        dimensionedScalar nud_;
        //- Optional diameter of the dispersed phase particles
        dimensionedScalar dd1_;
        dimensionedScalar dd2_;
        dimensionedScalar dd3_;
        dimensionedScalar dd4_;
        dimensionedScalar dd5_;
        dimensionedScalar dd6_;
        dimensionedScalar dd7_;
        dimensionedScalar dd8_;
         dimensionedScalar dt_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        const volVectorField& U_;
      
        const surfaceScalarField& phi_;

        volScalarField mu_;
      
      //  volScalarField alphas_;


public:

    TypeName("incompressibleNinePhaseInteractingMixture");

    incompressibleNinePhaseInteractingMixture ( const volVectorField& U, const surfaceScalarField& phi);

    //- Destructor
    virtual ~incompressibleNinePhaseInteractingMixture()
    {}


    // Member Functions

        //- Return const-access to the mixture viscosityModel
        const mixtureViscosityModel& muModel1() const
        {
            return muModel1_();
        }
        const mixtureViscosityModel& muModel2() const
        {
            return muModel2_();
        }
        const mixtureViscosityModel& muModel3() const
        {
            return muModel3_();
        }
        const mixtureViscosityModel& muModel4() const
        {
            return muModel4_();
        }
         const mixtureViscosityModel& muModel5() const
        {
            return muModel5_();
        }
         const mixtureViscosityModel& muModel6() const
        {
            return muModel6_();
        }
         const mixtureViscosityModel& muModel7() const
        {
            return muModel7_();
        }
         const viscosityModel& muModel8() const
        {
            return muModel8_();
        }

        //- Return const-access to the continuous-phase viscosityModel
        const viscosityModel& nucModel() const
        {
            return nucModel_();
        }

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rhod1() const
        {
            return rhod1_;
        }
        const dimensionedScalar& rhod2() const
        {
            return rhod2_;
        }
        const dimensionedScalar& rhod3() const
        {
            return rhod3_;
        }
        const dimensionedScalar& rhod4() const
        {
            return rhod4_;
        }
        
        const dimensionedScalar& rhod5() const
        {
            return rhod5_;
        }
        
         const dimensionedScalar& rhod6() const
        {
            return rhod6_;
        }
        
         const dimensionedScalar& rhod7() const
        {
            return rhod7_;
        }
         const dimensionedScalar& rhod8() const
        {
            return rhod8_;
        }
        //- Return const-access to continuous-phase density
        const dimensionedScalar& rhoc() const
        {
            return rhoc_;
        };

        //- Return the diameter of the dispersed-phase particles
        const dimensionedScalar& dd1() const
        {
            return dd1_;
        }
        const dimensionedScalar& dd2() const
        {
            return dd2_;
        }
        const dimensionedScalar& dd3() const
        {
            return dd3_;
        }
        const dimensionedScalar& dd4() const
        {
            return dd4_;
        }
        const dimensionedScalar& dd5() const
        {
            return dd5_;
        }
        const dimensionedScalar& dd6() const
        {
            return dd6_;
        }
        const dimensionedScalar& dd7() const
        {
            return dd7_;
        }
//        const dimensionedScalar& dd8() const
//        {
//            return dd8_;
//        }
//        
        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the dynamic mixture viscosity for patch
        virtual tmp<scalargpuField> mu(const label patchi) const
        {
            return mu_.boundaryField()[patchi];
        }

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const
        {
            return alpha1_*rhod1_ + alpha2_*rhod2_ + alpha3_*rhod3_ + alpha4_*rhod4_ + alpha5_*rhod5_+alpha6_*rhod6_ + alpha7_*rhod7_ + alpha8_*rhod8_ + alpha9_*rhoc_;
        }

        //- Return the mixture density for patch
        virtual tmp<scalargpuField> rho(const label patchi) const
        {
            return
                alpha1_.boundaryField()[patchi]*rhod1_.value()
              + alpha2_.boundaryField()[patchi]*rhod2_.value()
              + alpha3_.boundaryField()[patchi]*rhod3_.value()
              + alpha4_.boundaryField()[patchi]*rhod4_.value()
              + alpha5_.boundaryField()[patchi]*rhod5_.value()
              + alpha6_.boundaryField()[patchi]*rhod6_.value()
              + alpha7_.boundaryField()[patchi]*rhod7_.value()
              + alpha8_.boundaryField()[patchi]*rhod8_.value()
              + alpha9_.boundaryField()[patchi]*rhoc_.value();
        }

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
           
            volScalarField alphas = alpha1_+alpha2_+alpha3_+alpha4_+alpha5_+alpha6_;
            
            const volScalarField limitedAlpha8
            (
                "limitedAlpha8",
                min(max(alpha8_, scalar(0)), scalar(1))
            );
            const volScalarField limitedAlphas
            (
                "limitedAlphas",
                min(max(alphas, scalar(0)), scalar(1))
            );

            // Average kinematic viscosity calculated from dynamic viscosity
            return mu()/(limitedAlpha8*rhod8_ + (scalar(1) - limitedAlpha8 )*rhoc_ + limitedAlphas*rhod2_);
            
        }

        //- Return the mixture viscosity for patch
        virtual tmp<scalargpuField> nu(const label patchi) const
        {
            return mu_.boundaryField()[patchi]/rho(patchi);
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            ninePhaseInterfaceProperties::correct();
          
            const volScalarField limitedAlpha8
            (
                min(max(alpha8_, scalar(0)), scalar(1))
            );   
            
            volScalarField alphas = alpha1_+alpha2_+alpha3_+alpha4_+alpha5_+alpha6_;
     //       volScalarField alphal = alpha1_ + alpha2_ + alpha3_+ alpha4_ + alpha5_ + alpha6_ + alpha9;
            const volScalarField limitedAlphas
            (
                min(max(alphas, scalar(0)), scalar(0.62))
            );   
            
            mu_ =  (limitedAlpha8*(rhod8_)*nud_)  + (((scalar(1) -limitedAlpha8)*(rhoc_)*(nuc_)) * (pow((1.0 - ((limitedAlphas+1e-20)/0.63)),-1.55))); // * (1.0 + 2.5*limitedAlphas + 10.05*sqr(limitedAlphas) + 0.00273*exp(16.6*limitedAlphas)
          
        }
        
        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

using Foam::incompressibleNinePhaseInteractingMixture;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
