{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alpha7Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha7)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha7)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha7)
        );

        solve(alpha7Eqn);
        
        Info<< "Phase-7 volume fraction = "
            << alpha7.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha7) = " << min(alpha7).value()
            << "  Max(alpha7) = " << max(alpha7).value()
            << endl;
        
        
        
        tmp<surfaceScalarField> tphiAlphaUD(alpha7Eqn.flux());
        phiAlpha7 = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha7,
                phiAlpha7,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha7,
                phiAlpha7,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #endif

            phiAlpha7 += tphiAlphaCorr0();
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
                alpha7,
                alphaScheme
            )
          + fvc::flux
            (
                phir7,
                alpha7,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha7);
            volScalarField alpha10(alpha7);

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha7,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha7,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha7 += tphiAlphaCorr();
            }
            else
            {
                alpha7 = 0.5*alpha7 + 0.5*alpha10;
                phiAlpha7 += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha7 = tphiAlphaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve
            (
                alpha7,
                phi,
                phiAlpha7,
                mixture.alphaMax(),
                0
            );
            #else
            MULES::explicitSolve
            (
                alpha7,
                phi,
                phiAlpha7,
                mixture.alphaMax(),
                0
            );
            #endif
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha7 - tphiAlphaCorr0;
    }

   

    Info<< "Phase-7 volume fraction = "
        << alpha7.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha7) = " << min(alpha7).value()
        << "  Max(alpha7) = " << max(alpha7).value()
        << endl;
}
