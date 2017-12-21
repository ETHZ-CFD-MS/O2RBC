// create field on the RBC mesh
PtrList<volScalarField> Hb(RBCRegions.size());
PtrList<volScalarField> PO2_RBC(RBCRegions.size());

// Populate fluid field pointer lists
forAll(RBCRegions, i)
{
    Info<< "*** Reading/creating fields for RBC mesh "
        << RBCRegions[i].name() << nl << endl;

    Info<< "    Adding to Hb\n" << endl;
    Hb.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Hb",
                runTime.timeName(),
                RBCRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            RBCRegions[i]
        )
    );

    Info<< "    Adding to PO2_RBC\n" << endl;
    PO2_RBC.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "PO2_RBC",
                runTime.timeName(),
                RBCRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            RBCRegions[i],
            dimensionSet(0,0,0,0,0,0,0)
        )
    );
}

