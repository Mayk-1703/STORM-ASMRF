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

#include "Maninnen.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(Maninnen, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, Maninnen, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::Maninnen::Maninnen
(
    const dictionary& dict,
    const incompressibleNinePhaseInteractingMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    a_("a", dimless, dict.lookup("a")),
    a1_("a1", dimless, dict.lookup("a1")),
    g0_("g0", dimAcceleration, dict.lookup("g0")),
    dt_("dt", dimTime, 0.00001)
  //  residualAlpha_(dict.lookup("residualAlpha"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::Maninnen::~Maninnen()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::Maninnen::correct()
{
     volScalarField alphas = 1 - (alphad1_+alphad2_+alphad3_+alphad4_+alphad5_+alphad6_);
     const volScalarField limitedAlphas
     (
                "limitedAlphas",
                min(max(alphas, scalar(1e-10)), scalar(1))
     );
     const volScalarField limitedAlpha1
     (
                "limitedAlpha1",
                min(max(alphad1_, scalar(1e-10)), scalar(1))
     );
     const volScalarField limitedAlpha2
     (
                "limitedAlpha2",
                min(max(alphad2_, scalar(1e-10)), scalar(1))
     );
     const volScalarField limitedAlpha3
     (
                "limitedAlpha3",
                min(max(alphad3_, scalar(1e-10)), scalar(1))
     );
     const volScalarField limitedAlpha4
     (
                "limitedAlpha4",
                min(max(alphad4_, scalar(1e-10)), scalar(1))
     );
     const volScalarField limitedAlpha5
     (
                "limitedAlpha5",
                min(max(alphad5_, scalar(1e-10)), scalar(1))
     );
     const volScalarField limitedAlpha6
     (
                "limitedAlpha6",
                min(max(alphad6_, scalar(1e-10)), scalar(1))
     );
     
     const volScalarField SumSij1
     (
                "SumSij1",
                (((-3.50-(1.10*(dd2_/dd1_))-(1.02*(pow((dd2_/dd1_),2)))-(0.002*pow((dd2_/dd1_),3))) - a1_) * limitedAlpha2) +  (((-3.50-(1.10*(dd3_/dd1_))-(1.02*(pow((dd3_/dd1_),2)))-(0.002*pow((dd3_/dd1_),3))) - a1_) * limitedAlpha3)+ (((-3.50-(1.10*(dd4_/dd1_))-(1.02*(pow((dd4_/dd1_),2)))-(0.002*pow((dd4_/dd1_),3))) - a1_)*limitedAlpha4)+ (((-3.50-(1.10*(dd5_/dd1_))-(1.02*(pow((dd5_/dd1_),2)))-(0.002*pow((dd5_/dd1_),3))) - a1_) * limitedAlpha5)+ (((-3.50-(1.10*(dd6_/dd1_))-(1.02*(pow((dd6_/dd1_),2)))-(0.002*pow((dd6_/dd1_),3))) - a1_) * limitedAlpha6)
     );
     const volScalarField SumSij2
     (
                "SumSij2",
                (((-3.50-(1.10*(dd1_/dd2_))-(1.02*(pow((dd1_/dd2_),2)))-(0.002*pow((dd1_/dd2_),3))) - a1_) * limitedAlpha1) +  (((-3.50-(1.10*(dd3_/dd2_))-(1.02*(pow((dd3_/dd2_),2)))-(0.002*pow((dd3_/dd2_),3))) - a1_) * limitedAlpha3)+ (((-3.50-(1.10*(dd4_/dd2_))-(1.02*(pow((dd4_/dd2_),2)))-(0.002*pow((dd4_/dd2_),3))) - a1_) * limitedAlpha4)+ (((-3.50-(1.10*(dd5_/dd2_))-(1.02*(pow((dd5_/dd2_),2)))-(0.002*pow((dd5_/dd2_),3))) - a1_) * limitedAlpha5)+ (((-3.50-(1.10*(dd6_/dd2_))-(1.02*(pow((dd6_/dd2_),2)))-(0.002*pow((dd6_/dd2_),3))) - a1_) * limitedAlpha6)
     );
     
     const volScalarField SumSij3
     (
                "SumSij3",
                (((-3.50-(1.10*(dd1_/dd3_))-(1.02*(pow((dd1_/dd3_),2)))-(0.002*pow((dd1_/dd3_),3))) - a1_) * limitedAlpha1) +  (((-3.50-(1.10*(dd2_/dd3_))-(1.02*(pow((dd2_/dd3_),2)))-(0.002*pow((dd2_/dd3_),3))) - a1_) * limitedAlpha2)+ (((-3.50-(1.10*(dd4_/dd3_))-(1.02*(pow((dd4_/dd3_),2)))-(0.002*pow((dd4_/dd3_),3))) - a1_) * limitedAlpha4)+ (((-3.50-(1.10*(dd5_/dd3_))-(1.02*(pow((dd5_/dd3_),2)))-(0.002*pow((dd5_/dd3_),3))) - a1_) * limitedAlpha5)+ (((-3.50-(1.10*(dd6_/dd3_))-(1.02*(pow((dd6_/dd3_),2)))-(0.002*pow((dd6_/dd3_),3))) - a1_) * limitedAlpha6)
     );
     const volScalarField SumSij4
     (
                "SumSij4",
                (((-3.50-(1.10*(dd1_/dd4_))-(1.02*(pow((dd1_/dd4_),2)))-(0.002*pow((dd1_/dd4_),3))) - a1_) * limitedAlpha1) +  (((-3.50-(1.10*(dd2_/dd4_))-(1.02*(pow((dd2_/dd4_),2)))-(0.002*pow((dd2_/dd4_),3))) - a1_) * limitedAlpha2)+ (((-3.50-(1.10*(dd3_/dd4_))-(1.02*(pow((dd3_/dd4_),2)))-(0.002*pow((dd3_/dd4_),3))) - a1_) * limitedAlpha3)+ (((-3.50-(1.10*(dd5_/dd4_))-(1.02*(pow((dd5_/dd4_),2)))-(0.002*pow((dd5_/dd4_),3))) - a1_) * limitedAlpha5)+ (((-3.50-(1.10*(dd6_/dd4_))-(1.02*(pow((dd6_/dd4_),2)))-(0.002*pow((dd6_/dd4_),3))) - a1_) * limitedAlpha6)
     );
     const volScalarField SumSij5
     (
                "SumSij5",
                (((-3.50-(1.10*(dd1_/dd5_))-(1.02*(pow((dd1_/dd5_),2)))-(0.002*pow((dd1_/dd5_),3))) - a1_) * limitedAlpha1) +  (((-3.50-(1.10*(dd2_/dd5_))-(1.02*(pow((dd2_/dd5_),2)))-(0.002*pow((dd2_/dd5_),3))) - a1_) * limitedAlpha2)+ (((-3.50-(1.10*(dd4_/dd5_))-(1.02*(pow((dd4_/dd5_),2)))-(0.002*pow((dd4_/dd5_),3))) - a1_) * limitedAlpha4)+ (((-3.50-(1.10*(dd3_/dd5_))-(1.02*(pow((dd3_/dd5_),2)))-(0.002*pow((dd3_/dd5_),3))) - a1_) * limitedAlpha3)+ (((-3.50-(1.10*(dd6_/dd5_))-(1.02*(pow((dd6_/dd5_),2)))-(0.002*pow((dd6_/dd5_),3))) - a1_) * limitedAlpha6)
     );
     const volScalarField SumSij6
     (
                "SumSij6",
                (((-3.50-(1.10*(dd1_/dd6_))-(1.02*(pow((dd1_/dd6_),2)))-(0.002*pow((dd1_/dd6_),3))) - a1_) * limitedAlpha1) +  (((-3.50-(1.10*(dd2_/dd6_))-(1.02*(pow((dd2_/dd6_),2)))-(0.002*pow((dd2_/dd6_),3))) - a1_) * limitedAlpha2)+ (((-3.50-(1.10*(dd4_/dd6_))-(1.02*(pow((dd4_/dd6_),2)))-(0.002*pow((dd4_/dd6_),3))) - a1_) * limitedAlpha4)+ (((-3.50-(1.10*(dd3_/dd6_))-(1.02*(pow((dd3_/dd6_),2)))-(0.002*pow((dd3_/dd6_),3))) - a1_) * limitedAlpha3)+ (((-3.50-(1.10*(dd5_/dd6_))-(1.02*(pow((dd5_/dd6_),2)))-(0.002*pow((dd5_/dd6_),3))) - a1_) * limitedAlpha5)
     );
     
     Udm1_ = (dd1_* dd1_*(rhod1_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep1(),0.687))*(pow((limitedAlphas),a1_)*(1/(1+SumSij1))) )))) * (g0_ - Force() + a_*LiftForce1());
     Udm2_ = (dd2_* dd2_*(rhod2_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep2(),0.687))*(pow((limitedAlphas),a1_)*(1/(1+SumSij2))) )))) * (g0_ - Force() + a_*LiftForce2());
     Udm3_ = (dd3_* dd3_*(rhod3_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep3(),0.687))*(pow((limitedAlphas),a1_)*(1/(1+SumSij3))) )))) * (g0_ - Force() + a_*LiftForce3());
     Udm4_ = (dd4_* dd4_*(rhod4_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep4(),0.687))*(pow((limitedAlphas),a1_)*(1/(1+SumSij4))) )))) * (g0_ - Force() + a_*LiftForce4());
     Udm5_ = (dd5_* dd5_*(rhod5_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep5(),0.687))*(pow((limitedAlphas),a1_)*(1/(1+SumSij5))) )))) * (g0_ - Force() + a_*LiftForce5());
     Udm6_ = (dd6_* dd6_*(rhod6_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep6(),0.687))*(pow((limitedAlphas),a1_)*(1/(1+SumSij6))) )))) * (g0_ - Force() + a_*LiftForce6());
 //    Udm7_ = (dd7_* dd7_*(rhod7_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep7(),0.687)) )))) * (g0_ - Force() + LiftForce7());
 //    Udm8_ = (dd8_* dd8_*(rhod8_ - rho())/((scalar(18) * mu() * ((scalar(1)+scalar(0.15) * pow(Rep8(),0.687)) )))) * (g0_ - Force() + LiftForce8());
    
     Udm1_ = Udm1_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_;// - alphad7_ * Udm7_ - alphad8_ * Udm8_;
     Udm2_ = Udm2_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_;// - alphad7_ * Udm7_ - alphad8_ * Udm8_;
     Udm3_ = Udm3_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_;// - alphad7_ * Udm7_ - alphad8_ * Udm8_;
     Udm4_ = Udm4_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_;// - alphad7_ * Udm7_ - alphad8_ * Udm8_;
     Udm5_ = Udm5_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_;// - alphad7_ * Udm7_ - alphad8_ * Udm8_;
     Udm6_ = Udm6_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_;// - alphad7_ * Udm7_ - alphad8_ * Udm8_;
 //    Udm7_ = Udm7_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_ - alphad7_ * Udm7_ - alphad8_ * Udm8_;
 //    Udm8_ = Udm8_ - alphad1_*Udm1_ - alphad2_*Udm2_ - alphad3_ *Udm3_ - alphad4_ * Udm4_ - alphad5_ * Udm5_ - alphad6_ * Udm6_ - alphad7_ * Udm7_ - alphad8_ * Udm8_; 
     
     Info<<"Correct velocity"<<endl;
}


// ************************************************************************* //
