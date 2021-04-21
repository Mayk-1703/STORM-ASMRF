{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alpha4Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha4)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha4)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha4)
        );

        solve(alpha4Eqn);
        
        Info<< "Phase-4 volume fraction = "
            << alpha4.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha4) = " << min(alpha4).value()
            << "  Max(alpha4) = " << max(alpha4).value()
            << endl;
        
        
        
        tmp<surfaceScalarField> tphiAlphaUD(alpha4Eqn.flux());
        phiAlpha4 = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha4,
                phiAlpha4,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha4,
                phiAlpha4,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #endif

            phiAlpha4 += tphiAlphaCorr0();
        }

        // Cache the upwind-flux
        tphiAlphaCorr0 = tphiAlphaUD;
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> tphiAlphaUn
        (
            fvc::flux
            (
                phi,
                alpha4,
                alphaScheme
            )
          + fvc::flux
            (
                phir4,
                alpha4,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha4);
            volScalarField alpha10(alpha4);

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha4,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha4,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha4 += tphiAlphaCorr();
            }
            else
            {
                alpha4 = 0.5*alpha4 + 0.5*alpha10;
                phiAlpha4 += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha4 = tphiAlphaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve
            (
                alpha4,
                phi,
                phiAlpha4,
                mixture.alphaMax(),
                0
            );
            #else
            MULES::explicitSolve
            (
                alpha4,
                phi,
                phiAlpha4,
                mixture.alphaMax(),
                0
            );
            #endif
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha4 - tphiAlphaCorr0;
    }

   

    Info<< "Phase-4 volume fraction = "
        << alpha4.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha3) = " << min(alpha4).value()
        << "  Max(alpha3) = " << max(alpha4).value()
        << endl;
}
