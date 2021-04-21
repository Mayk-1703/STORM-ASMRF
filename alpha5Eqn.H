{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alpha5Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha5)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha5)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha5)
        );

        solve(alpha5Eqn);
        
        Info<< "Phase-5 volume fraction = "
            << alpha5.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha5) = " << min(alpha5).value()
            << "  Max(alpha5) = " << max(alpha5).value()
            << endl;
        
        
        
        tmp<surfaceScalarField> tphiAlphaUD(alpha5Eqn.flux());
        phiAlpha5 = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha5,
                phiAlpha5,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha5,
                phiAlpha5,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #endif

            phiAlpha5 += tphiAlphaCorr0();
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
                alpha5,
                alphaScheme
            )
          + fvc::flux
            (
                phir5,
                alpha5,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha5);
            volScalarField alpha10(alpha5);

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha5,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha5,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha5 += tphiAlphaCorr();
            }
            else
            {
                alpha5 = 0.5*alpha5 + 0.5*alpha10;
                phiAlpha5 += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha5 = tphiAlphaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve
            (
                alpha5,
                phi,
                phiAlpha5,
                mixture.alphaMax(),
                0
            );
            #else
            MULES::explicitSolve
            (
                alpha5,
                phi,
                phiAlpha5,
                mixture.alphaMax(),
                0
            );
            #endif
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha5 - tphiAlphaCorr0;
    }

   

    Info<< "Phase-5 volume fraction = "
        << alpha5.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha5) = " << min(alpha5).value()
        << "  Max(alpha5) = " << max(alpha5).value()
        << endl;
}
