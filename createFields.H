    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
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

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
  
    incompressibleNinePhaseInteractingMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());
    volScalarField& alpha3(mixture.alpha3());
    volScalarField& alpha4(mixture.alpha4());
    volScalarField& alpha5(mixture.alpha5());
    volScalarField& alpha6(mixture.alpha6());
    volScalarField& alpha7(mixture.alpha7());
    volScalarField& alpha8(mixture.alpha8());
    volScalarField& alpha9(mixture.alpha9());
    
    const dimensionedScalar& rho1 = mixture.rhod1();
    const dimensionedScalar& rho2 = mixture.rhod2();
    const dimensionedScalar& rho3 = mixture.rhod3();
    const dimensionedScalar& rho4 = mixture.rhod4();
    const dimensionedScalar& rho5 = mixture.rhod5();
    const dimensionedScalar& rho6 = mixture.rhod6();
    const dimensionedScalar& rho7 = mixture.rhod7();
    const dimensionedScalar& rho8 = mixture.rhod8();
    const dimensionedScalar& rho9 = mixture.rhoc();

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

    // Mixture density
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.rho()
    );
    
     volScalarField alphas
    (
        IOobject
        (
            "alphas",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7 + alpha8
    );

    // Mass flux
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );
     surfaceScalarField rhoPhi1
    (
        IOobject
        (
            "rhoPhi1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha1*rho1)*phi
    );
      surfaceScalarField rhoPhi2
    (
        IOobject
        (
            "rhoPhi2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha2*rho1)*phi
    );
      surfaceScalarField rhoPhi3
    (
        IOobject
        (
            "rhoPhi3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha3*rho1)*phi
    );
      surfaceScalarField rhoPhi4
    (
        IOobject
        (
            "rhoPhi4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha4*rho1)*phi
    );
      surfaceScalarField rhoPhi5
    (
        IOobject
        (
            "rhoPhi5",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha5*rho1)*phi
    );
      surfaceScalarField rhoPhi6
    (
        IOobject
        (
            "rhoPhi6",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha6*rho1)*phi
    );
    // Relative Velocity
    autoPtr<relativeVelocityModel> UdmModelPtr
    (
        relativeVelocityModel::New
        (
            transportProperties,
            mixture
        )
    );

    relativeVelocityModel& UdmModel(UdmModelPtr());


//     Construct compressible turbulence model
    autoPtr
    <
        CompressibleTurbulenceModel<incompressibleNinePhaseInteractingMixture>
    > turbulence
    (
        CompressibleTurbulenceModel<incompressibleNinePhaseInteractingMixture>
        ::New(rho, U, rhoPhi, mixture)
    );
//    
//      autoPtr<incompressible::turbulenceModel> turbulence
//    (
//        incompressible::turbulenceModel::New(U, phi, mixture)
//    );

    
    Info<< "Calculating field (g.h)f\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("gh", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );
    
   
    
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }


    // MULES Correction
    tmp<surfaceScalarField> tphiAlphaCorr0;
 