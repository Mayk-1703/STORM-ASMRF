/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "ninePhaseInterfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::ninePhaseInterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
namespace Foam
{
	struct interfacePropertiesCorrectContactAngleFunctor : public std::binary_function<scalar,scalar,scalar>{
		__HOST____DEVICE__
		scalar operator ()(const scalar& a12, const scalar& theta){
			return cos(acos(a12) - theta);
		}
	};
}


void Foam::ninePhaseInterfaceProperties::correctContactAngle
(
    surfaceVectorField::GeometricBoundaryField& nHatb
) const
{
    const volScalarField::GeometricBoundaryField& alpha8 =
        mixture_.alpha8().boundaryField();
    const volScalarField::GeometricBoundaryField& alpha9 =
        mixture_.alpha9().boundaryField();
    const volScalarField::GeometricBoundaryField& alpha3 =
        mixture_.alpha3().boundaryField();
    const volVectorField::GeometricBoundaryField& U =
        mixture_.U().boundaryField();

    const fvMesh& mesh = mixture_.U().mesh();
    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha8[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& a2cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha8[patchi]);

           

            scalargpuField twoPhaseAlpha2(max(a2cap, scalar(0)));
          
//            scalargpuField sumTwoPhaseAlpha
//            (
//                twoPhaseAlpha2  + SMALL
//            );

            twoPhaseAlpha2 /= 1;//sumTwoPhaseAlpha;
          
            fvsPatchVectorField& nHatp = nHatb[patchi];

            scalargpuField theta
            (
                convertToRad
              * (
                   twoPhaseAlpha2*(180 - a2cap.theta(U[patchi], nHatp))
                )
            );

            vectorgpuField nf(boundary[patchi].nf());

            // Reset nHatPatch to correspond to the contact angle

            scalargpuField a12(nHatp & nf);

            scalargpuField b1(cos(theta));

            scalargpuField b2(nHatp.size());

            thrust::transform(a12.begin(),a12.end(),theta.begin(),b2.begin(),
                              interfacePropertiesCorrectContactAngleFunctor());

            scalargpuField det(1.0 - a12*a12);

            scalargpuField a((b1 - a12*b2)/det);
            scalargpuField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;

            nHatp /= (mag(nHatp) + deltaN_.value());
        }
    }
}


void Foam::ninePhaseInterfaceProperties::calculateK()
{
    const volScalarField& alpha8 = mixture_.alpha8();

    const fvMesh& mesh = alpha8.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(alpha8));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(nHatfv.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    // volVectorField nHat = gradAlpha/(mag(gradAlpha) + deltaN_);
    // nHat.boundaryField() = nHatfv.boundaryField();
    // K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ninePhaseInterfaceProperties::ninePhaseInterfaceProperties
(
    const ninePhaseMixture& mixture
)
:
    mixture_(mixture),
    cAlpha_
    (
        readScalar
        (
            mixture.U().mesh().solverDict
            (
                mixture_.alpha8().name()
            ).lookup("cAlpha")
        )
    ),
    sigma89_
    (
      dimensionedScalar("sigma", dimensionSet(1, 0, -2, 0, 0) , 0.072)
    ),
//      sigma89_("sigma89", dimensionSet(1, 0, -2, 0, 0), mixture),
//    sigma13_("sigma13", dimensionSet(1, 0, -2, 0, 0), mixture),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mixture.U().mesh().V()), 1.0/3.0)
    ),
    nHatf_
    (
        IOobject
        (
            "nHatf",
            mixture.alpha8().time().timeName(),
            mixture.alpha8().mesh()
        ),
        mixture.alpha8().mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),
    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            mixture.alpha8().time().timeName(),
            mixture.alpha8().mesh()
        ),
        mixture.alpha8().mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::ninePhaseInterfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(mixture_.alpha8());
}


Foam::tmp<Foam::volScalarField>
Foam::ninePhaseInterfaceProperties::nearInterface() const
{
    return max
    (
        pos(mixture_.alpha8() - 0.01)*pos(0.99 - mixture_.alpha8()),
        pos(mixture_.alpha9() - 0.01)*pos(0.99 - mixture_.alpha9())
    );
}


// ************************************************************************* //
