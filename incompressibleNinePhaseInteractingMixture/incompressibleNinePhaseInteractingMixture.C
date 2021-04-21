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

\*---------------------------------------------------------------------------*/

#include "incompressibleNinePhaseInteractingMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleNinePhaseInteractingMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleNinePhaseInteractingMixture::
incompressibleNinePhaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    ninePhaseMixture(U.mesh(), *this, U),
    ninePhaseInterfaceProperties
    (
        static_cast<ninePhaseMixture&>(*this)
    ),
    muModel1_
    (
        mixtureViscosityModel::New
        (
            "mu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
     muModel2_
    (
        mixtureViscosityModel::New
        (
            "mu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
     muModel3_
    (
        mixtureViscosityModel::New
        (
            "mu3",
            subDict(phase3Name_),
            U,
            phi
        )
    ),
     muModel4_
    (
        mixtureViscosityModel::New
        (
            "mu4",
            subDict(phase4Name_),
            U,
            phi
        )
    ),
    
    muModel5_
    (
        mixtureViscosityModel::New
        (
            "mu5",
            subDict(phase5Name_),
            U,
            phi
        )
    ),
    
    muModel6_
    (
        mixtureViscosityModel::New
        (
            "mu6",
            subDict(phase6Name_),
            U,
            phi
        )
    ),
    
    muModel7_
    (
        mixtureViscosityModel::New
        (
            "mu7",
            subDict(phase7Name_),
            U,
            phi
        )
    ),

     muModel8_
    (
       viscosityModel::New
        (
            "nud",
            subDict(phase8Name_),
            U,
            phi
        )
    ),

    
    nucModel_
    (
        viscosityModel::New
        (
            "nuc",
            subDict(phase9Name_),
            U,
            phi
        )
    ),

    rhod1_("rho", dimDensity, muModel1_->viscosityProperties().lookup("rho")),
    rhod2_("rho", dimDensity, muModel2_->viscosityProperties().lookup("rho")),
    rhod3_("rho", dimDensity, muModel3_->viscosityProperties().lookup("rho")),
    rhod4_("rho", dimDensity, muModel4_->viscosityProperties().lookup("rho")),
    rhod5_("rho", dimDensity, muModel5_->viscosityProperties().lookup("rho")),
    rhod6_("rho", dimDensity, muModel6_->viscosityProperties().lookup("rho")),
    rhod7_("rho", dimDensity, muModel7_->viscosityProperties().lookup("rho")),
    rhod8_("rho", dimDensity, muModel8_->viscosityProperties().lookup("rho")),
    rhoc_("rho", dimDensity, nucModel_->viscosityProperties().lookup("rho")),
    
    nuc_("nu", dimensionSet(0, 2, -1, 0, 0), nucModel_->viscosityProperties().lookup("nu")),
    nud_("nu", dimensionSet(0, 2, -1, 0, 0), muModel8_->viscosityProperties().lookup("nu")),   
    
    dd1_
    (
        "d",
        dimLength,
        muModel1_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    dd2_
    (
        "d",
        dimLength,
        muModel2_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    dd3_
    (
        "d",
        dimLength,
        muModel3_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    dd4_
    (
        "d",
        dimLength,
        muModel4_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    dd5_
    (
        "d",
        dimLength,
        muModel5_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    dd6_
    (
        "d",
        dimLength,
        muModel6_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
    dd7_
    (
        "d",
        dimLength,
        muModel7_->viscosityProperties().lookupOrDefault("d", 0.0)
    ),
//    dd8_
//    (
//        "d",
//        dimLength,
//        muModel8_->viscosityProperties().lookupOrDefault("d", 0.0)
//    ),
    alphaMax_(muModel1_->viscosityProperties().lookupOrDefault("alphaMax", 0.63)),

    U_(U),
    phi_(phi),
    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )

{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::incompressibleNinePhaseInteractingMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            muModel1_().read(subDict(phase1Name_))
         && muModel2_().read(subDict(phase2Name_))
         && muModel3_().read(subDict(phase3Name_))
         && muModel4_().read(subDict(phase4Name_)) 
         && muModel5_().read(subDict(phase5Name_)) 
         && muModel6_().read(subDict(phase6Name_)) 
         && muModel7_().read(subDict(phase7Name_))   
         && muModel8_().read(subDict(phase8Name_)) 
         && nucModel_().read(subDict(phase9Name_))
        )
        {
            muModel1_->viscosityProperties().lookup("rho") >> rhod1_;
            muModel2_->viscosityProperties().lookup("rho") >> rhod2_;
            muModel3_->viscosityProperties().lookup("rho") >> rhod3_;
            muModel4_->viscosityProperties().lookup("rho") >> rhod4_;
            muModel5_->viscosityProperties().lookup("rho") >> rhod5_;
            muModel6_->viscosityProperties().lookup("rho") >> rhod6_;
            muModel7_->viscosityProperties().lookup("rho") >> rhod7_;
            muModel8_->viscosityProperties().lookup("rho") >> rhod8_;
            nucModel_->viscosityProperties().lookup("rho") >> rhoc_;

            dd1_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel1_->viscosityProperties().lookupOrDefault("d", 0)
            );
             dd2_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel2_->viscosityProperties().lookupOrDefault("d", 0)
            );
             dd3_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel3_->viscosityProperties().lookupOrDefault("d", 0)
            );
             dd4_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel4_->viscosityProperties().lookupOrDefault("d", 0)
            );
            
            dd5_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel5_->viscosityProperties().lookupOrDefault("d", 0)
            );
            
            dd6_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel6_->viscosityProperties().lookupOrDefault("d", 0)
            );
            
            dd7_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel7_->viscosityProperties().lookupOrDefault("d", 0)
            );
            
//            dd8_ = dimensionedScalar
//            (
//                "d",
//                dimLength,
//                muModel8_->viscosityProperties().lookupOrDefault("d", 0)
//            );
            alphaMax_ =
                muModel1_->viscosityProperties().lookupOrDefault
                (
                    "alphaMax",
                    0.63
                );

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
