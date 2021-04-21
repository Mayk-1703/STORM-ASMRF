{
     surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    surfaceScalarField phiAlpha1
    (
        IOobject
        (
            "phiAlpha1",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    
    surfaceScalarField phiAlpha2
    (
        IOobject
        (
            "phiAlpha2",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    
    surfaceScalarField phiAlpha3
    (
        IOobject
        (
            "phiAlpha3",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    
    surfaceScalarField phiAlpha4
    (
        IOobject
        (
            "phiAlpha4",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    
      surfaceScalarField phiAlpha5
    (
        IOobject
        (
            "phiAlpha5",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    
      surfaceScalarField phiAlpha6
    (
        IOobject
        (
            "phiAlpha6",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    
      surfaceScalarField phiAlpha7
    (
        IOobject
        (
            "phiAlpha7",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );
    
      surfaceScalarField phiAlpha8
    (
        IOobject
        (
            "phiAlpha8",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );

    surfaceScalarField phir1
    (
        mesh.Sf() & fvc::interpolate(UdmModel.Udm1())
    );
    
     surfaceScalarField phir2
    (
        mesh.Sf() & fvc::interpolate(UdmModel.Udm2())
    );
     surfaceScalarField phir3
    (
        mesh.Sf() & fvc::interpolate(UdmModel.Udm3())
    );
     surfaceScalarField phir4
    (
        mesh.Sf() & fvc::interpolate(UdmModel.Udm4())
    );
     surfaceScalarField phir5
    (
        mesh.Sf() & fvc::interpolate(UdmModel.Udm5())
    );
     surfaceScalarField phir6
    (
        mesh.Sf() & fvc::interpolate(UdmModel.Udm6())
    );
//     surfaceScalarField phir7
//    (
//        mesh.Sf() & fvc::interpolate(UdmModel.Udm7())
//    );
//     surfaceScalarField phir8
//    (
//        mesh.Sf() & fvc::interpolate(UdmModel.Udm8())
//    );
    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField phiAlpha1Sum
        (
            IOobject
            (
                "phiAlpha1Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
         surfaceScalarField phiAlpha2Sum
        (
            IOobject
            (
                "phiAlpha2Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
         surfaceScalarField phiAlpha3Sum
        (
            IOobject
            (
                "phiAlpha3Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
         surfaceScalarField phiAlpha4Sum
        (
            IOobject
            (
                "phiAlpha4Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
        
         surfaceScalarField phiAlpha5Sum
        (
            IOobject
            (
                "phiAlpha5Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
        
         surfaceScalarField phiAlpha6Sum
        (
            IOobject
            (
                "phiAlpha6Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
        
         surfaceScalarField phiAlpha7Sum
        (
            IOobject
            (
                "phiAlpha7Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
        
         surfaceScalarField phiAlpha8Sum
        (
            IOobject
            (
                "phiAlpha8Sum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha8, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
           
            #include "alpha8Eqn.H"
            phiAlpha8Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha8;
            #include "alpha1Eqn.H"
            phiAlpha1Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha1;
            #include "alpha2Eqn.H"
            phiAlpha2Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha2;
            #include "alpha3Eqn.H"
            phiAlpha3Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha3;
            #include "alpha4Eqn.H"
            phiAlpha4Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha4;
            #include "alpha5Eqn.H"
            phiAlpha5Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha5;
            #include "alpha6Eqn.H"
            phiAlpha6Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha6;
//            #include "alpha7Eqn.H"
//            phiAlpha7Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha7;
          
        }
        phiAlpha1 = phiAlpha1Sum;
        phiAlpha2 = phiAlpha2Sum;
        phiAlpha3 = phiAlpha3Sum;
        phiAlpha4 = phiAlpha4Sum;
        phiAlpha5 = phiAlpha5Sum;
        phiAlpha6 = phiAlpha6Sum;
//        phiAlpha7 = phiAlpha7Sum;
        phiAlpha8 = phiAlpha8Sum;
        
    }
    else
    {
         #include "alpha8Eqn.H"
        #include "alpha1Eqn.H"
        #include "alpha2Eqn.H"
        #include "alpha3Eqn.H"
        #include "alpha4Eqn.H"
        #include "alpha5Eqn.H"
        #include "alpha6Eqn.H"
//        #include "alpha7Eqn.H"
      
    }

    // Apply the diffusion term separately to allow implicit solution
    // and boundedness of the explicit advection
    {
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1) - fvc::ddt(alpha1)
         - fvm::laplacian(turbulence->nut(), alpha1)
        );

        alpha1Eqn.solve(mesh.solver("alpha1Diffusion"));

        phiAlpha1 += alpha1Eqn.flux();
        
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
        
        fvScalarMatrix alpha2Eqn
        (
            fvm::ddt(alpha2) - fvc::ddt(alpha2)
         - fvm::laplacian(turbulence->nut(), alpha2)
        );

        alpha2Eqn.solve(mesh.solver("alpha2Diffusion"));

        phiAlpha2 += alpha2Eqn.flux();
        
        Info<< "Phase-2 volume fraction = "
            << alpha2.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha2) = " << min(alpha2).value()
            << "  Max(alpha2) = " << max(alpha2).value()
            << endl;
            
        fvScalarMatrix alpha3Eqn
        (
            fvm::ddt(alpha3) - fvc::ddt(alpha3)
         - fvm::laplacian(turbulence->nut(), alpha3)
        );

        alpha3Eqn.solve(mesh.solver("alpha3Diffusion"));

        phiAlpha3 += alpha3Eqn.flux();
        
        Info<< "Phase-3 volume fraction = "
            << alpha3.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha3) = " << min(alpha3).value()
            << "  Max(alpha3) = " << max(alpha3).value()
            << endl;
        
        fvScalarMatrix alpha4Eqn
        (
            fvm::ddt(alpha4) - fvc::ddt(alpha4)
         - fvm::laplacian(turbulence->nut(), alpha4)
        );

        alpha4Eqn.solve(mesh.solver("alpha4Diffusion"));

        phiAlpha4 += alpha4Eqn.flux();
        
        Info<< "Phase-4 volume fraction = "
            << alpha4.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha4) = " << min(alpha4).value()
            << "  Max(alpha4) = " << max(alpha4).value()
            << endl;
        
         fvScalarMatrix alpha5Eqn
        (
            fvm::ddt(alpha5) - fvc::ddt(alpha5)
         - fvm::laplacian(turbulence->nut(), alpha5)
        );

        alpha5Eqn.solve(mesh.solver("alpha5Diffusion"));

        phiAlpha5 += alpha5Eqn.flux();
        
        Info<< "Phase-5 volume fraction = "
            << alpha5.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha5) = " << min(alpha5).value()
            << "  Max(alpha5) = " << max(alpha5).value()
            << endl;
        
         fvScalarMatrix alpha6Eqn
        (
            fvm::ddt(alpha6) - fvc::ddt(alpha6)
         - fvm::laplacian(turbulence->nut(), alpha6)
        );

        alpha6Eqn.solve(mesh.solver("alpha4Diffusion"));

        phiAlpha6 += alpha6Eqn.flux();
        
        Info<< "Phase-6 volume fraction = "
            << alpha6.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha6) = " << min(alpha6).value()
            << "  Max(alpha6) = " << max(alpha6).value()
            << endl;
        
  
        phiAlpha = phiAlpha1 + phiAlpha2 + phiAlpha3 + phiAlpha4 + phiAlpha5 + phiAlpha6 ;
        
        rhoPhi1 == phiAlpha1*rho1;
        rhoPhi2 == phiAlpha2*rho1;
        rhoPhi3 == phiAlpha3*rho1;
        rhoPhi4 == phiAlpha4*rho1;
        rhoPhi5 == phiAlpha5*rho1;
        rhoPhi6 == phiAlpha6*rho1;
  //      divTauDm = rho1*phiAlpha1*fvc::interpolate(UdmModel.Udm1()) + rho1*phiAlpha2*fvc::interpolate(UdmModel.Udm2()) +rho1*phiAlpha3*fvc::interpolate(UdmModel.Udm3()) + rho1*phiAlpha4*fvc::interpolate(UdmModel.Udm4()) + rho1*phiAlpha5*fvc::interpolate(UdmModel.Udm5()) + rho1*phiAlpha6*fvc::interpolate(UdmModel.Udm6());
        alphas =    alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha8;// + alpha7 + alpha8;
        alpha9 = (1.0 - alphas);
    }
        
    rhoPhi = phiAlpha*(rho1) + (phiAlpha8*rho8) + (phi-phiAlpha - phiAlpha8)*rho9;
    rho == mixture.rho();
}
