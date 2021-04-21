{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alpha2Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha2)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha2)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha2)
        );

        solve(alpha2Eqn);
        
        Info<< "Phase-2 volume fraction = "
            << alpha2.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha1) = " << min(alpha2).value()
            << "  Max(alpha1) = " << max(alpha2).value()
            << endl;
        
        
        
        tmp<surfaceScalarField> tphiAlphaUD(alpha2Eqn.flux());
        phiAlpha2 = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha2,
                phiAlpha2,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha2,
                phiAlpha2,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #endif

            phiAlpha2 += tphiAlphaCorr0();
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
                alpha2,
                alphaScheme
            )
          + fvc::flux
            (
                phir2,
                alpha2,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha2);
            volScalarField alpha10(alpha2);

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha2,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha2,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha2 += tphiAlphaCorr();
            }
            else
            {
                alpha2 = 0.5*alpha2 + 0.5*alpha10;
                phiAlpha2 += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha2 = tphiAlphaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve
            (
                alpha2,
                phi,
                phiAlpha2,
                mixture.alphaMax(),
                0
            );
            #else
            MULES::explicitSolve
            (
                alpha2,
                phi,
                phiAlpha2,
                mixture.alphaMax(),
                0
            );
            #endif
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha2 - tphiAlphaCorr0;
    }

   

    Info<< "Phase-2 volume fraction = "
        << alpha2.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha2) = " << min(alpha2).value()
        << "  Max(alpha2) = " << max(alpha2).value()
        << endl;
}
