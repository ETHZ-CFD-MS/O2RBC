// The advection-diffusion-reaction diffusion for PO2 is solved using Godunov
// splitting (see Leveque, Finite Volume Methods for Hyperbolic Problems, Chap.
// 17).
// First, the advection equation is solved (with some diffusion to stabilize it,
// optionally). 
// Then, the diffusion-reaction equation is solved.
//
// This file is for the second step, which is in a loop for solving iteratively.

{
    // diffusion term
    fvScalarMatrix diffusionTerm(PO2, dimless/dimTime);
    if (O2AxialDiffusion)
    {
        diffusionTerm = fvm::laplacian(kappaO2*alpha, PO2);
    }
    else
    {
        volTensorField DO2
        (
            IOobject
            (
                "DO2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            kappaO2*diffusionTensor
        );
        diffusionTerm = fvm::laplacian(DO2*alpha, PO2);
    }

    // PO2 equation without consumption and reaction term
    fvScalarMatrix PO2Eqn
    (
        alpha*fvm::ddt(PO2)
      - diffusionTerm
    );

    // add linearized consumption term
    volScalarField Mp = PO2Crit/sqr(PO2 + PO2Crit);
    volScalarField Mc = sqr(PO2/(PO2 + PO2Crit));
    PO2Eqn -= - inTissue*O2ConsumptionRate()
              * (fvm::Sp(Mp, PO2) + Mc);

    // add linearized reaction term
    if (sum(inRBCMask).value() > 0)
    {
        volScalarField Rp = DC.ddP(PO2, HbEuler, inRBCEuler);
        volScalarField Rc = DC.reactionTerm(PO2, HbEuler, inRBCEuler)
           - PO2*DC.ddP(PO2, HbEuler, inRBCEuler);
        PO2Eqn -= inRBCMask*O2ReactionFactor
                    * (fvm::Sp(Rp, PO2) + Rc);
    }
                
    // solve equation
    PO2Eqn.solve(mesh.solver("PO2"));
    catchSmallOrNanValues(PO2);

    PO2Phi = fvc::snGrad(PO2);
}

