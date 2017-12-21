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

dimensionedScalar kappa_O2_RBC
(
    transportProperties.lookup("kappa_O2_RBC")
);

dimensionedScalar PO2_crit
(
    transportProperties.lookup("PO2_crit")
);

dimensionedScalar RBC_velocity
(
    transportProperties.lookup("RBC_velocity")
);

dimensionedScalar O2_consumption_rate
(
    transportProperties.lookup("O2_consumption_rate")
);

dimensionedScalar PO2_RBC_inlet
(
    transportProperties.lookup("PO2_RBC_inlet")
);

dimensionedScalar PO2_plasma
(
    transportProperties.lookup("PO2_plasma")
);

dimensionedVector U_plasma
(
    "U_plasma",
    dimensionSet(0,1,-1,0,0,0,0),
    vector(RBC_velocity.value(), 0, 0)
);

dimensionedScalar dissociation_rate
(
    transportProperties.lookup("dissociation_rate")
);

dimensionedScalar P_50
(
    transportProperties.lookup("P_50")
);

dimensionedScalar hill_exponent
(
    transportProperties.lookup("hill_exponent")
);

dimensionedScalar kappa_Hb
(
    transportProperties.lookup("kappa_Hb")
);

dimensionedScalar alpha_RBC
(
    transportProperties.lookup("alpha_RBC")
);

dimensionedScalar vol_mol_O2
(
    transportProperties.lookup("vol_mol_O2")
);

dimensionedScalar N_Hb
(
    transportProperties.lookup("N_Hb")
);

// Variables defined on the Eulerian mesh

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
#include "createPhi.H"

    // Variable for the update of alpha
    Info<< "Reading field alpha_no_RBC\n" << endl;
    volScalarField alpha_no_RBC
    (
        IOobject
        (
            "alpha_no_RBC",
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
        alpha_no_RBC
    );

    // Variable for the update of kappa_O2
    Info<< "Reading field kappa_O2_no_RBC\n" << endl;
    volScalarField kappa_O2_no_RBC
    (
        IOobject
        (
            "kappa_O2_no_RBC",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Creating field kappa_O2\n" << endl;
    volScalarField kappa_O2
    (
        IOobject
        (
            "kappa_O2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        kappa_O2_no_RBC
    );

    Info<< "Reading field in_tissue\n" << endl;
    volScalarField in_tissue
    (
        IOobject
        (
            "in_tissue",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
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
        dimensionSet(0,0,0,0,0,0,0)
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
        0.0*PO2
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
        dimensionSet(0,0,0,0,0,0,0)
    );

    // Derivative of consumption term w.r.t. PO2
    Info<< "Creating field Mp\n" << endl;
    volScalarField Mp
    (
        IOobject
        (
            "Mp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0,0,0)
    );
    // Constant part of consumption term
    Info<< "Creating field Mc\n" << endl;
    volScalarField Mc
    (
        IOobject
        (
            "Mc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0,0,0)
    );
    // Derivative of reaction term w.r.t. PO2
    Info<< "Creating field Rp\n" << endl;
    volScalarField Rp
    (
        IOobject
        (
            "Rp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0,0,0)
    );
    // Constant part of reaction term
    Info<< "Creating field Rc\n" << endl;
    volScalarField Rc
    (
        IOobject
        (
            "Rc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0,0,0)
    );


#   include "createPhi.H"

