
#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "incompressibleNinePhaseInteractingMixture.H"
#include "relativeVelocityModel.H"
#include "turbulenceModel.H"
#include "CompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    pimpleControl pimple(mesh);
       #include "createPrghCorrTypes.H"
    #include "correctPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
            UdmModel.correct();
            #include "alphaEqnSubCycle.H"
            mixture.correct();
            #include "UEqn.H"
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << "s"
            << "  ClockTime = " << runTime.elapsedClockTime() << "s"
            << nl << endl;
    }
    Info<< "End\n" << endl;
    return 0;
}
// ************************************************************************* //
