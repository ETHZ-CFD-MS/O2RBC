// set values at time 0 for the following fields
// Hb_euler, in_RBC_euler, kappa_O2
//
// Those values are defined using the following fields, which are read from
// files:
// PO2, Hb, kappa_O2_no_RBC

Hb_euler = 0.0;
PO2Phi = fvc::snGrad(PO2);

// interpolate * Hb      -> Hb_euler
//             * in_RBC  -> in_RBC_euler
interpolator.interpolateRBCToEuler("Hb", Hb_euler);
interpolator.interpolateRBCToEuler("in_RBC", in_RBC_euler);

// compute alpha and kappa_O2
alpha    = in_RBC_euler*alpha_RBC    + (1-in_RBC_euler)*alpha_no_RBC;
kappa_O2 = in_RBC_euler*kappa_O2_RBC + (1-in_RBC_euler)*kappa_O2_no_RBC;

// compute oxygen concentration
C = alpha*PO2;

// write fields at time 0
Hb_euler.write();
in_RBC_euler.write();
kappa_O2.write();
alpha.write();
C.write();

forAll(RBCs, i)
{
    RBCs[i].in_RBC().write();
    RBCs[i].PO2   ().write();
}


