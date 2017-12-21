Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

bool O2AxialDiffusion = transportProperties.lookupOrDefault("O2AxialDiffusion", true);

tensor diffusionTensor(1,0,0,0,1,0,0,0,1);
if (!O2AxialDiffusion)
{
    diffusionTensor.xx() = 0.0;
    Info<< "Disabling oxygen diffusion in the axial direction" << endl;
}

// Read the oxygen consumption rate. 
expressionSource<scalar> O2ConsumptionRate
(
    transportProperties.subDict("O2ConsumptionRate"),
    mesh
);

dimensionedScalar PO2Crit
(
    transportProperties.lookup("PO2Crit")
);

dimensionedScalar alphaRBC
(
    transportProperties.lookup("alphaRBC")
);

dimensionedScalar alphaPlasma
(
    transportProperties.lookup("alphaPlasma")
);

dimensionedScalar kappaO2RBC
(
    transportProperties.lookup("kappaO2RBC")
);

dimensionedScalar kappaHb
(
    transportProperties.lookup("kappaHb")
);

dimensionedScalar PO2RBCInlet
(
    transportProperties.lookup("PO2RBCInlet")
);

dimensionedScalar PO2PlasmaInlet
(
    transportProperties.lookup("PO2PlasmaInlet")
);

dimensionedScalar dissociationRate
(
    transportProperties.lookup("dissociationRate")
);

dimensionedScalar VMolO2
(
    transportProperties.lookup("VMolO2")
);

dimensionedScalar NHb
(
    transportProperties.lookup("NHb")
);

// Variables defined on the Eulerian mesh

    Info<< "Creating field alphaNoRBC\n" << endl;
    volScalarField alphaNoRBC
    (
        IOobject
        (
            "alphaNoRBC",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Creating field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaNoRBC
    );

    Info<< "Creating field kappaO2NoRBC\n" << endl;
    volScalarField kappaO2NoRBC
    (
        IOobject
        (
            "kappaO2NoRBC",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Creating field kappaO2\n" << endl;
    volScalarField kappaO2
    (
        IOobject
        (
            "kappaO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        kappaO2NoRBC
    );

    Info<< "Creating field DalphaSurf\n" << endl;
    surfaceScalarField DalphaSurf
    (
        IOobject
        (
            "DalphaSurf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(kappaO2*alpha)
    );

    Info<< "Reading field PO2\n" << endl;
    volScalarField PO2
    (
        IOobject
        (
            "PO2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field PO2Adv\n" << endl;
    volScalarField PO2Adv
    (
        IOobject
        (
            "PO2Adv",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field PO2Phi\n" << endl;
    surfaceScalarField PO2Phi
    (
        IOobject
        (
            "PO2Phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        inv(dimLength)
    );

    Info<< "Reading field C\n" << endl;
    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field HbEuler\n" << endl;
    volScalarField HbEuler
    (
        IOobject
        (
            "HbEuler",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("HbEuler", dimless, 0.0)
    );

    Info<< "Creating field inRBCEuler\n" << endl;
    volScalarField inRBCEuler
    (
        IOobject
        (
            "inRBCEuler",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("inRBCEuler", dimless, 0.0)
    );

    Info<< "Creating field oldInRBCEuler\n" << endl;
    volScalarField oldInRBCEuler
    (
        IOobject
        (
            "oldInRBCEuler",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("oldInRBCEuler", dimless, 0.0)
    );

    Info<< "Creating field inRBCMask\n" << endl;
    volScalarField inRBCMask
    (
        IOobject
        (
            "inRBCMask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("inRBCMask", dimless, 0.0)
    );

    Info<< "Creating field inTissue\n" << endl;
    volScalarField inTissue
    (
        IOobject
        (
            "inTissue",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Creating field O2ConsumptionRateOutput\n" << endl;
    volScalarField O2ConsumptionRateOutput
    (
        IOobject
        (
            "O2ConsumptionRateOutput",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        inv(dimVolume*dimTime)
    );

    Info<< "Creating field xiRBC\n" << endl;
    volScalarField xiRBC
    (
        IOobject
        (
            "xiRBC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("xiRBC", dimless, 0.0)
    );

    Info<< "Creating phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimVolume/dimTime
    );

