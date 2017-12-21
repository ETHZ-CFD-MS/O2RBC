{
    volScalarField& Hb     = RBCs.Hb();
    volScalarField& PO2RBC = RBCs.PO2();

    fvScalarMatrix HbEqn
    (
        fvm::ddt(Hb)
      - fvm::laplacian(kappaHb, Hb)
    );

    volScalarField RpHb = -(DC.ddS(PO2RBC, Hb));
    volScalarField RcHb = -(DC.reactionTerm(PO2RBC, Hb) - Hb*DC.ddS(PO2RBC, Hb));

    HbEqn -= dissociationRate * (fvm::Sp(RpHb,Hb) + RcHb);
    HbEqn.solve();

    bool hbVerbose = catchSmallOrNanValues(RBCs.Hb());

    if (hbVerbose)
    {
        Pout<< "Hb boundary Field = " << Hb.boundaryField() << endl;
    }
}

