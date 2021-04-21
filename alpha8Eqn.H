{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    
      surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phic.boundaryField()[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> tphiAlpha;

    
    if (MULESCorr)
    {
        fvScalarMatrix alpha8Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha8)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha8)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha8)
        );

        solve(alpha8Eqn);
        
        Info<< "Phase-8 volume fraction = "
            << alpha8.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha8) = " << min(alpha8).value()
            << "  Max(alpha8) = " << max(alpha8).value()
            << endl;
        
        
        
        tmp<surfaceScalarField> tphiAlphaUD(alpha8Eqn.flux());
        phiAlpha8 = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha8,
                phiAlpha8,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha8,
                phiAlpha8,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #endif

            phiAlpha8 += tphiAlphaCorr0();
        }

        // Cache the upwind-flux
        tphiAlphaCorr0 = tphiAlphaUD;
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir8(phic*mixture.nHatf());
        
        tmp<surfaceScalarField> tphiAlphaUn
        (
            fvc::flux
            (
                phi,
                alpha8,
                alphaScheme
            )
            + fvc::flux
            (
               -fvc::flux(-phir8, alpha9, alpharScheme),
                alpha8,
                alpharScheme
            )
            
//            + fvc::flux
//            (
//               -fvc::flux(-phir8, alpha1, alpharScheme),
//                alpha8,
//                alpharScheme
//            )
//            + fvc::flux
//            (
//               -fvc::flux(-phir8, alpha2, alpharScheme),
//                alpha8,
//                alpharScheme
//            )
//            + fvc::flux
//            (
//               -fvc::flux(-phir8, alpha3, alpharScheme),
//                alpha8,
//                alpharScheme
//            )
//            + fvc::flux
//            (
//               -fvc::flux(-phir8, alpha4, alpharScheme),
//                alpha8,
//                alpharScheme
//            )
//            + fvc::flux
//            (
//               -fvc::flux(-phir8, alpha5, alpharScheme),
//                alpha8,
//                alpharScheme
//            )
//            + fvc::flux
//            (
//               -fvc::flux(-phir8, alpha6, alpharScheme),
//                alpha8,
//                alpharScheme
//            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha8);
            volScalarField alpha10(alpha8);

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha8,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha8,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha8 += tphiAlphaCorr();
            }
            else
            {
                alpha8 = 0.5*alpha8 + 0.5*alpha10;
                phiAlpha8 += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha8 = tphiAlphaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve
            (
                alpha8,
                phi,
                phiAlpha8,
                mixture.alphaMax(),
                0
            );
            #else
            MULES::explicitSolve
            (
                alpha8,
                phi,
                phiAlpha8,
                mixture.alphaMax(),
                0
            );
            #endif
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha8 - tphiAlphaCorr0;
    }

   

    Info<< "Phase-8 volume fraction = "
        << alpha8.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha8) = " << min(alpha8).value()
        << "  Max(alpha8) = " << max(alpha8).value()
        << endl;
}
