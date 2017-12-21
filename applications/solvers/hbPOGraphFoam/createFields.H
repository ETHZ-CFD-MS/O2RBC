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

dimensionedScalar dissociationRate
(
    transportProperties.lookup("dissociationRate")
);

dimensionedScalar kappaHb
(
    transportProperties.lookup("kappaHb")
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

    Info<< "Creating field alpha\n" << endl;
    regionDependentField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        transportProperties.subDict("alpha"),
        vesselRegions
    );

    Info<< "Creating field kappa_O2\n" << endl;
    regionDependentField kappa_O2
    (
        IOobject
        (
            "kappa_O2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        transportProperties.subDict("kappa"),
        vesselRegions
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
        dimensionSet(0,-1,0,0,0,0,0)
    );

    Info<< "Creating field upwindCells\n" << endl;
    labelIOField upwindCells
    (
        IOobject
        (
            "upwindCells",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.nCells()
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

    Info<< "Creating field Hb_euler\n" << endl;
    volScalarField Hb_euler
    (
        IOobject
        (
            "Hb_euler",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Hb_euler", dimless, 0.0)
    );

    Info<< "Creating field in_RBC_euler\n" << endl;
    volScalarField in_RBC_euler
    (
        IOobject
        (
            "in_RBC_euler",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("in_RBC_euler", dimless, 0.0)
    );

    Info<< "Creating field in_RBC_mask\n" << endl;
    volScalarField in_RBC_mask
    (
        IOobject
        (
            "in_RBC_mask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("in_RBC_mask", dimless, 0.0)
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
        dimensionSet(0,-3,-1,0,0,0,0)
    );

