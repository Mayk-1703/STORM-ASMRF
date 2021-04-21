/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "ninePhaseMixture.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ninePhaseMixture::ninePhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict,
    const volVectorField& U
)
:
    phase1Name_(wordList(dict.lookup("phases"))[0]),
    phase2Name_(wordList(dict.lookup("phases"))[1]),
    phase3Name_(wordList(dict.lookup("phases"))[2]),
    phase4Name_(wordList(dict.lookup("phases"))[3]),
    phase5Name_(wordList(dict.lookup("phases"))[4]),
    phase6Name_(wordList(dict.lookup("phases"))[5]),
    phase7Name_(wordList(dict.lookup("phases"))[6]),
    phase8Name_(wordList(dict.lookup("phases"))[7]),
    phase9Name_(wordList(dict.lookup("phases"))[8]),
   
    U_(U),
   
    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha3_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha4_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase4Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha5_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase5Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha6_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase6Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha7_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase7Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    alpha8_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase8Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    alpha9_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase9Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - alpha1_ - alpha2_ - alpha3_ - alpha4_ - alpha5_ - alpha6_ - alpha7_ - alpha8_
    )
{}


// ************************************************************************* //
