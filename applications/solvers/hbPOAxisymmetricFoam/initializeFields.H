// Set values of dependent fields as a function of base fields.
// Write these fields to disk.

HbEuler = 0.0;
PO2Phi = fvc::snGrad(PO2);

// interpolate * Hb      -> Hb_euler
//             * in_RBC  -> in_RBC_euler
interpolator.interpolateRBCToEuler("Hb", HbEuler);
interpolator.interpolateRBCToEuler("in_RBC", inRBCEuler);

// compute alpha and kappa_O2
alpha   = inRBCEuler*alphaRBC   + (1 - inRBCEuler)*alphaNoRBC;
kappaO2 = inRBCEuler*kappaO2RBC + (1 - inRBCEuler)*kappaO2NoRBC;
DalphaSurf = fvc::interpolate(kappaO2*alpha);

// compute oxygen concentration
C = alpha*PO2;

// write fields at time 0
HbEuler.write();
inRBCEuler.write();
kappaO2.write();
alpha.write();
C.write();
RBCs.in_RBC().write();
RBCs.PO2().write();

