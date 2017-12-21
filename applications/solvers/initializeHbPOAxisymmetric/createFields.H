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

Info<< "Reading geometricProperties\n" << endl;

IOdictionary geometricProperties
(
    IOobject
    (
        "geometricProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar RBCVolume
(
    geometricProperties.lookup("RBCVolume")
);


// Variables defined on the Eulerian mesh

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
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    PO2
);


