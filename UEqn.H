    // Solve the Momentum equation
     Info<< "In U eqn  "<<endl;
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(U) 
      - fvc::div(rhoPhi1,UdmModel.Udm1()) 
      - fvc::div(rhoPhi2,UdmModel.Udm2())
      - fvc::div(rhoPhi3,UdmModel.Udm3())
      - fvc::div(rhoPhi4,UdmModel.Udm4())
      - fvc::div(rhoPhi5,UdmModel.Udm5())
      - fvc::div(rhoPhi6,UdmModel.Udm6())
       ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                   mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
