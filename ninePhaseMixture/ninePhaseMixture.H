/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::twoPhaseMixture

Description
    A two-phase mixture model

SourceFiles
    twoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef ninePhaseMixture_H
#define ninePhaseMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class ninePhaseMixture

{
protected:

    // Protected data

        word phase1Name_;
        word phase2Name_;
        word phase3Name_;
        word phase4Name_;
        word phase5Name_;
        word phase6Name_;
        word phase7Name_;
        word phase8Name_;
        word phase9Name_;
        
        volScalarField alpha1_;
        volScalarField alpha2_;
        volScalarField alpha3_;
        volScalarField alpha4_;
        volScalarField alpha5_;
        volScalarField alpha6_;
        volScalarField alpha7_;
        volScalarField alpha8_;
        volScalarField alpha9_;
        const volVectorField& U_;

public:

    // Constructors

        //- Construct from components
        ninePhaseMixture
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const volVectorField& U
        );


    //- Destructor
    ~ninePhaseMixture()
    {}


    // Member Functions
        
        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }
        
         const word& phase3Name() const
        {
            return phase3Name_;
        }
        
         const word& phase4Name() const
        {
            return phase4Name_;
        }
        
         const word& phase5Name() const
        {
            return phase5Name_;
        }
         const word& phase6Name() const
        {
            return phase6Name_;
        }
         const word& phase7Name() const
        {
            return phase7Name_;
        }
         const word& phase8Name() const
        {
            return phase8Name_;
        }
         const word& phase9Name() const
        {
            return phase9Name_;
        }

        //- Return the phase-fraction of phase 1
       const volVectorField& U() const
        {
            return U_;
        }
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 2
        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha2()
        {
            return alpha2_;
        }
        
        const volScalarField& alpha3() const
        {
            return alpha3_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha3()
        {
            return alpha3_;
        }
        
        const volScalarField& alpha4() const
        {
            return alpha4_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha4()
        {
            return alpha4_;
        }
        
        const volScalarField& alpha5() const
        {
            return alpha5_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha5()
        {
            return alpha5_;
        }
        
         const volScalarField& alpha6() const
        {
            return alpha6_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha6()
        {
            return alpha6_;
        }
        
         const volScalarField& alpha7() const
        {
            return alpha7_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha7()
        {
            return alpha7_;
        }
        
         const volScalarField& alpha8() const
        {
            return alpha8_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha8()
        {
            return alpha8_;
        }
        
         const volScalarField& alpha9() const
        {
            return alpha9_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha9()
        {
            return alpha9_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
