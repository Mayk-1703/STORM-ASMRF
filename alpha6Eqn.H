{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alpha6Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha6)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha6)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha6)
        );

        solve(alpha6Eqn);
        
        Info<< "Phase-6 volume fraction = "
            << alpha6.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha6) = " << min(alpha6).value()
            << "  Max(alpha6) = " << max(alpha6).value()
            << endl;
        
        
        
        tmp<surfaceScalarField> tphiAlphaUD(alpha6Eqn.flux());
        phiAlpha6 = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha6,
                phiAlpha6,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha6,
                phiAlpha6,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );
            #endif

            phiAlpha6 += tphiAlphaCorr0();
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
                alpha6,
                alphaScheme
            )
          + fvc::flux
            (
                phir6,
                alpha6,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha6);
            volScalarField alpha10(alpha6);

            #ifdef LTSSOLVE
            MULES::LTScorrect
            (
                alpha6,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #else
            MULES::correct
            (
                alpha6,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );
            #endif

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha6 += tphiAlphaCorr();
            }
            else
            {
                alpha6 = 0.5*alpha6 + 0.5*alpha10;
                phiAlpha6 += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha6 = tphiAlphaUn;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve
            (
                alpha6,
                phi,
                phiAlpha6,
                mixture.alphaMax(),
                0
            );
            #else
            MULES::explicitSolve
            (
                alpha6,
                phi,
                phiAlpha6,
                mixture.alphaMax(),
                0
            );
            #endif
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha6 - tphiAlphaCorr0;
    }

   

    Info<< "Phase-6 volume fraction = "
        << alpha6.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha6) = " << min(alpha6).value()
        << "  Max(alpha6) = " << max(alpha6).value()
        << endl;
}
