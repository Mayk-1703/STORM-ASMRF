{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alpha3Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha3)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha3)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha3)
        );

        solve(alpha3Eqn);
        
        Info<< "Phase-3 volume fraction = "
            << alpha3.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha3) = " << min(alpha3).value()
            << "  Max(alpha3) = " << max(alpha3).value()
            << endl;
        
        
        
        tmp<surfaceScalarField> tphiAlphaUD(alpha3Eqn.flux());
        phiAlpha3 = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha3,
                phiAlpha3,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha3,
                phiAlpha3,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #endif

            phiAlpha3 += tphiAlphaCorr0();
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
                alpha3,
                alphaScheme
            )
          + fvc::flux
            (
                phir3,
                alpha3,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha3);
            volScalarField alpha10(alpha3);

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha3,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha3,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha3 += tphiAlphaCorr();
            }
            else
            {
                alpha3 = 0.5*alpha3 + 0.5*alpha10;
                phiAlpha3 += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha3 = tphiAlphaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve
            (
                alpha3,
                phi,
                phiAlpha3,
                mixture.alphaMax(),
                0
            );
            #else
            MULES::explicitSolve
            (
                alpha3,
                phi,
                phiAlpha3,
                mixture.alphaMax(),
                0
            );
            #endif
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha3 - tphiAlphaCorr0;
    }

   

    Info<< "Phase-3 volume fraction = "
        << alpha3.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha3) = " << min(alpha3).value()
        << "  Max(alpha3) = " << max(alpha3).value()
        << endl;
}
