// The advection-diffusion-reaction diffusion for PO2 is solved using Godunov
// splitting (see Leveque, Finite Volume Methods for Hyperbolic Problems, Chap.
// 17).
// First, the advection equation is solved.
// Then, the diffusion-reaction equation is solved.
//
// This file contains is for the first step.

// Solve advection equation for the oxygen concentration.
// Since oxygen concentration is the conserved variable, we use this quantity.
// However, the corresponding equation with alpha*P (where alpha takes values
// *after* RBC displacement) yields almost the same results (< 0.3% difference). 
// For now, we stick to the formulation with C since it is more rigorous and 
// would (probably) work with other time stepping schemes.

{
    C = alphaPrevious*PO2;

    Info<< "Min/max of C before PO2Adv equation: " << min(C) << ", " << max(C) << endl;
    fvScalarMatrix CAdvEqn
    (
        fvm::ddt(C)
      + fvm::div(phi, C)
    );

    CAdvEqn.solve(mesh.solver("CAdv"));

    Info<< "Min/max of C after PO2Adv equation: " << min(C) << ", " << max(C) << endl;

    catchSmallOrNanValues(C);

    PO2 = C/alpha.field();
    PO2Adv = PO2;

    PO2.oldTime() = PO2;
}

