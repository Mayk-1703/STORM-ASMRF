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

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const dictionary& dict,
    const incompressibleNinePhaseInteractingMixture& mixture
)
:
    mixture_(mixture),
    U_(mixture.U()),
    alphac_(mixture.alpha9()),
    alphad1_(mixture.alpha1()),
    alphad2_(mixture.alpha2()),
    alphad3_(mixture.alpha3()),
    alphad4_(mixture.alpha4()),
    alphad5_(mixture.alpha5()),
    alphad6_(mixture.alpha6()),
    alphad7_(mixture.alpha7()),
    alphad8_(mixture.alpha8()),
    rhoc_(mixture.rhoc()),
    rhod1_(mixture.rhod1()),
    rhod2_(mixture.rhod2()),
    rhod3_(mixture.rhod3()),
    rhod4_(mixture.rhod4()),
    rhod5_(mixture.rhod5()),
    rhod6_(mixture.rhod6()),
    rhod7_(mixture.rhod7()),
    rhod8_(mixture.rhod8()),
    mu_(mixture.mu()),
    dd1_(mixture.dd1()),
    dd2_(mixture.dd2()),
    dd3_(mixture.dd3()),
    dd4_(mixture.dd4()),
    dd5_(mixture.dd5()),
    dd6_(mixture.dd6()),
    dd7_(mixture.dd7()),
//    dd8_(mixture.dd8()), 
    Rep1_
    (
        IOobject
        (
            "Rep1",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep1", dimless, 0)
    ),
    
    Rep2_
    (
        IOobject
        (
            "Rep2",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep2", dimless, 0)
    ),
    
    Rep3_
    (
        IOobject
        (
            "Rep3",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep3", dimless, 0)
    ),
    
    Rep4_
    (
        IOobject
        (
            "Rep4",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep4", dimless, 0)
    ),
     Rep5_
    (
        IOobject
        (
            "Rep5",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep5", dimless, 0)
    ),
     Rep6_
    (
        IOobject
        (
            "Rep6",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep6", dimless, 0)
    ),
     Rep7_
    (
        IOobject
        (
            "Rep7",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep7", dimless, 0)
    ),
     Rep8_
    (
        IOobject
        (
            "Rep8",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphac_.mesh(),
        dimensionedScalar ("Rep8", dimless, 0)
    ),
    
    DeltaT_
    (
        dimensionedScalar("DeltaT", dimTime , 0.00001)
    ),
    
    MixForce_
    (
        IOobject
        (
            "MixForce",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("MixForce", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    LiftForce_
    (
        IOobject
        (
            "LiftForce",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("LiftForce", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
     LiftForce1_
    (
        IOobject
        (
            "LiftForce1",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("LiftForce1", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
     LiftForce2_
    (
        IOobject
        (
            "LiftForce2",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("LiftForce2", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
     LiftForce3_
    (
        IOobject
        (
            "LiftForce3",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("LiftForce3", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
     LiftForce4_
    (
        IOobject
        (
            "LiftForce4",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("LiftForce4", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
     LiftForce5_
    (
        IOobject
        (
            "LiftForce5",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("LiftForce5", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
     LiftForce6_
    (
        IOobject
        (
            "LiftForce6",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("LiftForce6", dimAcceleration, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    Udm1_
    (
        IOobject
        (
            "Udm1",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm1", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    ),
    Udm2_
    (
        IOobject
        (
            "Udm2",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm2", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    Udm3_
    (
        IOobject
        (
            "Udm3",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm3", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    Udm4_
    (
        IOobject
        (
            "Udm4",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm4", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    Udm5_
    (
        IOobject
        (
            "Udm5",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm5", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    Udm6_
    (
        IOobject
        (
            "Udm6",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm6", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    Udm7_
    (
        IOobject
        (
            "Udm7",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm7", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    ),
    
    Udm8_
    (
        IOobject
        (
            "Udm8",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("Udm8", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const incompressibleNinePhaseInteractingMixture& mixture
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relativeVelocityModel::New"
            "("
                "const dictionary&"
            ")"
        )   << "Unknown time scale model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid time scale model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return
        autoPtr<relativeVelocityModel>
        (
            cstrIter()
            (
                dict.subDict(modelType + "Coeffs"),
                mixture
            )
        );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<volScalarField> Foam::relativeVelocityModel::rho() const
{
    return alphac_ * rhoc_ + alphad1_ * rhod1_ + alphad2_ * rhod2_ + alphad3_ * rhod3_ + alphad4_ * rhod4_ + alphad5_ * rhod5_ + alphad6_ * rhod6_ + alphad7_ * rhod7_ + alphad8_ * rhod8_ ;
}
tmp<volScalarField> Foam::relativeVelocityModel::mu() const
{
    
    return mu_;//(scalar(3.8) * pow((scalar(1) - (alphad_)/scalar(0.62)),-1.55));
}

tmp<volScalarField> Foam::relativeVelocityModel::Rep1() const
{
    Rep1_ = (rhod1_ * dd1_ * mag(Udm1_) / mu());
    return Rep1_;
}
tmp<volScalarField> Foam::relativeVelocityModel::Rep2() const
{
    Rep2_ = (rhod2_ * dd2_ * mag(Udm2_) / mu());
    return Rep2_;
}
tmp<volScalarField> Foam::relativeVelocityModel::Rep3() const
{
    Rep3_ = (rhod3_ * dd3_ * mag(Udm3_) / mu());
    return Rep3_;
}
tmp<volScalarField> Foam::relativeVelocityModel::Rep4() const
{
    Rep4_ = (rhod4_ * dd4_ * mag(Udm4_) / mu());
    return Rep4_;
}
tmp<volScalarField> Foam::relativeVelocityModel::Rep5() const
{
    Rep5_ = (rhod5_ * dd5_ * mag(Udm5_) / mu());
    return Rep5_;
}
tmp<volScalarField> Foam::relativeVelocityModel::Rep6() const
{
    Rep6_ = (rhod6_ * dd6_ * mag(Udm6_) / mu());
    return Rep6_;
}
//tmp<volScalarField> Foam::relativeVelocityModel::Rep7() const
//{
//    Rep7_ = (rhod7_ * dd7_ * mag(Udm7_) / mu());
//    return Rep7_;
//}
//tmp<volScalarField> Foam::relativeVelocityModel::Rep8() const
//{
//    Rep8_ = (rhod8_ * dd8_ * mag(Udm8_) / mu());
//    return Rep8_;
//}


tmp<volVectorField> Foam::relativeVelocityModel::PrevU() const
{
    return ((U_-U_.oldTime()));
}

tmp<volVectorField> Foam::relativeVelocityModel::LiftForce1() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce1_ = 0.75 * (rhoc_/(rhod1_-rho())) * (omega ^ Udm1_);
    return LiftForce1_ ; 
}
tmp<volVectorField> Foam::relativeVelocityModel::LiftForce2() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce2_ = 0.75 * (rhoc_/(rhod2_-rho())) * (omega ^ Udm2_);
    return LiftForce2_ ; 
}
tmp<volVectorField> Foam::relativeVelocityModel::LiftForce3() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce3_ = 0.75 * (rhoc_/(rhod3_-rho())) * (omega ^ Udm3_);
    return LiftForce3_ ; 
}
tmp<volVectorField> Foam::relativeVelocityModel::LiftForce4() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce4_ = 0.75 * (rhoc_/(rhod4_-rho())) * (omega ^ Udm4_);
    return LiftForce4_ ; 
}
tmp<volVectorField> Foam::relativeVelocityModel::LiftForce5() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce5_ = 0.75 * (rhoc_/(rhod5_-rho())) * (omega ^ Udm5_);
    return LiftForce5_ ; 
}
tmp<volVectorField> Foam::relativeVelocityModel::LiftForce6() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce6_ = 0.75 * (rhoc_/(rhod6_-rho())) * (omega ^ Udm6_);
    return LiftForce6_ ; 
}
tmp<volVectorField> Foam::relativeVelocityModel::LiftForce7() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce_ = 0.75 * (rhoc_/(rhod7_-rho())) * (omega ^ Udm7_);
    return LiftForce_ ; 
}

tmp<volVectorField> Foam::relativeVelocityModel::LiftForce8() const
{
    tmp<volVectorField> omega;
    omega = fvc::curl(U_);
    LiftForce_ = 0.75 * (rhoc_/(rhod8_-rho())) * (omega ^ Udm8_);
    return LiftForce_ ; 
}

tmp<volVectorField> Foam::relativeVelocityModel::Force() const
{
    tmp<volTensorField> GradU;
 //   dimensionedScalar dT ("dT",dimensionSet(0 0 -1 0 0),1));
//    dT = 1/DeltaT_;
    //tmp<volVectorField> MixGradient;
    //dimensionedScalar DeltaT = (runTime.deltaTValue());
    GradU = fvc::grad(U_);
    MixForce_ = (PrevU() / DeltaT_);
    MixForce_ = MixForce_ +  (GradU & U_);   
    return (MixForce_);
}


tmp<volSymmTensorField> Foam::relativeVelocityModel::tauDm() const
{
    volScalarField betac(alphac_*rhoc_ + alphad8_*rhod8_);
    volScalarField betad1(alphad1_*rhod1_);// + alphad2_*rhod2_ + alphad3_ * rhod3_ + _alphad4_*rhod4 );
    volScalarField betad2(alphad2_*rhod2_);
    volScalarField betad3(alphad3_*rhod3_);
    volScalarField betad4(alphad4_*rhod4_);
    volScalarField betad5(alphad5_*rhod5_);
    volScalarField betad6(alphad6_*rhod6_);
    
    volVectorField Ucm((betad1*(Udm1_)+betad2*(Udm2_)+betad3*(Udm3_)+betad4*(Udm4_)+betad5*(Udm5_)+betad6*(Udm6_))/(rho()*(alphac_+alphad8_)));

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            "tauDm",
            betad1*sqr(Udm1_) + betad2*sqr(Udm2_) + betad3*sqr(Udm3_) + betad4*sqr(Udm4_) + betad5*sqr(Udm5_) + betad6*sqr(Udm6_)   + betac*sqr(Ucm)
        )
    );
}
// ************************************************************************* //