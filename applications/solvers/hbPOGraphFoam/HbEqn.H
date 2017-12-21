{
    volScalarField& Hb      = RBCs.Hb();
    volScalarField& PO2_RBC = RBCs.PO2();

    // Info<< "PO2 RBC   :  min = " << min(PO2_RBC) << ", "
                    // << " max = " << max(PO2_RBC) << endl;
    // Info<< "Hb        :  min = " << min(Hb) << ", "
                    // << " max = " << max(Hb) << endl;

    // solve hemoglobin equation
    fvScalarMatrix HbEqn
    (
        fvm::ddt(Hb)
      - fvm::laplacian(kappaHb, Hb)
    );

    // linearize reaction terms
    volScalarField Rp_Hb = -(DC.ddS(PO2_RBC, Hb));
    volScalarField Rc_Hb = -(DC.reactionTerm(PO2_RBC, Hb) - Hb*DC.ddS(PO2_RBC, Hb));

    HbEqn -= dissociationRate * (fvm::Sp(Rp_Hb,Hb) + Rc_Hb); // implicit version

    HbEqn.solve();

    bool hbVerbose = catchSmallOrNanValues(RBCs.Hb());

    if (hbVerbose)
    {
        Pout<< "Hb boundary Field = " << Hb.boundaryField() << endl;
        // Pout<< "Hb0 boundary Field = " << Hb.oldTime().boundaryField() << endl;
        // Pout<< "HbEqn.A(): " << HbEqn.A() << endl;
        // Pout<< "HbEqn.D(): " << HbEqn.D() << endl;
        // Pout<< "HbEqn.H(): " << HbEqn.H() << endl;
        // Pout<< "HbEqn.source(): " << HbEqn.source() << endl;
    }
}


