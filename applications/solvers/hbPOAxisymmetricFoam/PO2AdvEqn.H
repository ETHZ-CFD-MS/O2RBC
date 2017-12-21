// The advection-diffusion-reaction diffusion for PO2 is solved using Godunov
// splitting (see Leveque, Finite Volume Methods for Hyperbolic Problems, Chap.
// 17).
// First, the advection equation is solved.
// Then, the diffusion-reaction equation is solved.
//
// This file contains deals with the advection step. The numerical method
// employed here is based on Section 2.2.4 of the PhD dissertation of Adrien
// Luecker available at
// 
// https://www.research-collection.ethz.ch/handle/20.500.11850/181551
//
// It is also published in the journal publication "The relative influence of
// hematocrit and red blood velocity on oxygen transport from capillaries to
// tissue" by Luecker, Secomb, Weber and Jenny. Microcirculation (2016),
// available at
//
// https://dx.doi.org/10.1111/micc.12337

// cell-centered interpolation weights for the oxygen concentration
volScalarField CWeights
(
    IOobject
    (
        "CWeights",
        mesh.time().timeName(),
        mesh, 
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CWeight", dimless, 0.0)
);

// interpolation coefficients associated to the face owners
surfaceScalarField CLambdas
(
    IOobject
    (
        "CLambdas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CLambda", dimless, 0.0)
);

// interpolation coefficients associated to the face neighbours
surfaceScalarField CYs
(
    IOobject
    (
        "CYs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CY", dimless, 0.0)
);

vector UVector(1, 0, 0);
scalar motionTol = 1e-6;

forAll(mesh.C(), cI)
{
    scalar deltaX = (mesh.C()[1] - mesh.C()[0]).x();
    scalar cfl = moverPtr->getRBCVelocity()*runTime.deltaT().value()/deltaX;
    if (mag(inRBCEuler[cI] - oldInRBCEuler[cI]) > motionTol) // RBC membrane moves in cell
    {
        scalar gammaRight;
        scalar alphaRight;
        scalar alphaLeft;
        if (inRBCEuler[cI] > oldInRBCEuler[cI]) // "right" RBC border moves
        {
            gammaRight = 1 - oldInRBCEuler[cI]/xiRBC[cI];
            alphaRight = alphaPlasma.value();
            alphaLeft = alphaRBC.value();
        }
        else // "left" RBC border
        {
            gammaRight = oldInRBCEuler[cI]/xiRBC[cI];
            alphaRight = alphaRBC.value();
            alphaLeft = alphaPlasma.value();
        }
        CWeights[cI] = (
                         xiRBC[cI]
                         * (  alphaRight*min(1, gammaRight/cfl) 
                            + alphaLeft*max(0, 1 - gammaRight/cfl)
                           )
                         + (1 - xiRBC[cI])*alphaPlasma.value()
                       ) / alphaPrevious[cI];
    }
    else // no region change
    {
        CWeights[cI] = 1.0;
    }
}

// set values of CWeights on the boundary
forAll(CWeights.boundaryField(), patchI)
{
    const fvPatch& patch = mesh.boundary()[patchI];
    fvPatchField<scalar>& pCw = CWeights.boundaryFieldRef()[patchI];
    // For boundary faces with outflux, use the cell centered weights.
    // For boundary faces with influx, set it to zero
    pCw = pos(phi.boundaryField()[patchI])*patch.patchInternalField(CWeights);
}

// set CLambdas and CYs on internal faces
forAll(CLambdas, faceI)
{
    if (phi[faceI] > 0)
    {
        CLambdas[faceI] = CWeights[mesh.owner()[faceI]];
        CYs     [faceI] = 0;
    }
    else if (phi[faceI] < 0)
    {
        CLambdas[faceI] = 0;
        CYs     [faceI] = CWeights[mesh.neighbour()[faceI]];
    }
}

// exchange processor patch fields for those whose patch neighbour field
// will be needed
forAll(CWeights.boundaryField(), patchI)
{
    if (isA<processorFvPatch>(mesh.boundary()[patchI]))
    {
        CWeights.boundaryFieldRef()[patchI].initEvaluate();
        C.boundaryFieldRef()[patchI].initEvaluate();
    }
}
forAll(CWeights.boundaryField(), patchI)
{
    if (isA<processorFvPatch>(mesh.boundary()[patchI]))
    {
        CWeights.boundaryFieldRef()[patchI].evaluate();
        C.boundaryFieldRef()[patchI].evaluate();
    }
}

// compute the surface-centered interpolation weights
forAll(CLambdas.boundaryField(), patchI)
{
    fvsPatchField<scalar>& pCLambdas = CLambdas.boundaryFieldRef()[patchI];
    fvsPatchField<scalar>& pCYs      = CYs     .boundaryFieldRef()[patchI];
    if (pCLambdas.coupled())
    {
        scalarField neighbourCWeights = CWeights.boundaryField()[patchI].patchNeighbourField();
        scalarField internalCWeights = CWeights.boundaryField()[patchI].patchInternalField();
        forAll(pCLambdas, faceI)
        {
            // inflow boundary
            if (phi.boundaryField()[patchI][faceI] < -VSMALL)
            {
                pCLambdas[faceI] = 0.0;
                pCYs     [faceI] = neighbourCWeights[faceI];
            }
            // outflow boundary
            else if (phi.boundaryField()[patchI][faceI] > VSMALL)
            {
                pCLambdas[faceI] = internalCWeights[faceI];
                pCYs     [faceI] = 0.0;
            }
        }
    }
    else
    {
        pCLambdas = CWeights.boundaryField()[patchI];
        pCYs      = CWeights.boundaryField()[patchI];
    }
}

// compute explicitly the correct inflow boundary field for C
forAll(C.boundaryField(), patchI)
{
    const fvPatch& patch = mesh.boundary()[patchI];
    if (!patch.coupled())
    {
        const labelList internalCellIndices(patch.patchInternalField(identity(mesh.nPoints())));
        fvPatchField<scalar>& pC = C.boundaryFieldRef()[patchI];
        // for boundary faces with outflux, use the previously computed weights
        pC = pos(phi.boundaryField()[patchI])*CLambdas.boundaryField()[patchI]
             * patch.patchInternalField(C);
        // for boundary faces with influx, compute by hand the boundary value of C
        vectorField faceNormals = mesh.boundary()[patchI].nf();
        forAll(pC, faceI)
        {
            label cI = internalCellIndices[faceI];
            if (mag(UVector/mag(UVector) & faceNormals[faceI]) > motionTol
                && phi.boundaryField()[patchI][faceI] < -VSMALL)
            {
                scalar deltaX = (mesh.C()[1] - mesh.C()[0]).x();
                scalar cfl = moverPtr->getRBCVelocity()*runTime.deltaT().value()/deltaX;
                // compute C at the left and right of the RBC interface in the
                // ghost cell, and the RBC volume fraction in the ghost cell
                bool membraneCrossing = false;
                scalar gammaRightGhost;
                scalar CRight;
                scalar CLeft;

                if (xiRBC[cI] > 0.0)
                {
                    // "right" RBC border crosses cell face
                    if (inRBCEuler[cI]/xiRBC[cI] > SMALL && oldInRBCEuler[cI] < SMALL) 
                    {
                        // Info<< "PO2AdvEqn.H: Right RBC border crossing, cI = " << cI << endl;
                        membraneCrossing = true;
                        gammaRightGhost = cfl - inRBCEuler[cI]/xiRBC[cI];
                        CRight = alphaPlasma.value()*PO2PlasmaInlet.value();
                        CLeft = alphaRBC.value()*PO2RBCInlet.value();
                    }
                    // "left" RBC border crosses cell face
                    else if (inRBCEuler[cI] < oldInRBCEuler[cI] - motionTol
                          && oldInRBCEuler[cI]/xiRBC[cI] > 1. - motionTol)
                    {
                        // Info<< "PO2AdvEqn.H: Left RBC border crossing, cI = " << cI << endl;
                        membraneCrossing = true;
                        gammaRightGhost = 1 + inRBCEuler[cI]/xiRBC[cI] - cfl;
                        CRight = alphaRBC.value()*PO2RBCInlet.value();
                        CLeft = alphaPlasma.value()*PO2PlasmaInlet.value();
                    }
                }
                if (membraneCrossing)
                {
                    pC[faceI] = xiRBC[cI]
                              * (  CRight*min(1, gammaRightGhost/cfl) 
                                 + CLeft*max(0, 1 - gammaRightGhost/cfl)
                                )
                              + (1 - xiRBC[cI])*alphaPlasma.value()*PO2PlasmaInlet.value();
                }
                else
                {
                    // if the face is in the RBC
                    if (xiRBC[cI] > 0.0
                     && (inRBCEuler[cI]/xiRBC[cI] > 1. - motionTol
                         || inRBCEuler[cI] > oldInRBCEuler[cI] + motionTol))
                    {
                        pC[faceI] = xiRBC[cI]*alphaRBC.value()*PO2RBCInlet.value() 
                                  + (1 - xiRBC[cI])*alphaPlasma.value()*PO2PlasmaInlet.value();
                    }
                    // if the face is in the plasma
                    else
                    {
                        pC[faceI] = alphaPlasma.value()*PO2PlasmaInlet.value();
                    }
                }
            }
        }
    }
}

// read the computation scheme for the divergence operator
Istream& schemeData(mesh.divScheme("div(phi,C)"));
word schemeName(schemeData); // read first word from stream
tmp<surfaceInterpolationScheme<scalar> > tDivScheme
(
    surfaceInterpolationScheme<scalar>::New
    (
        mesh, 
        phi, 
        schemeData
    )
);

// interpolate C to cell faces
tmp<surfaceScalarField> tCs
(
    tDivScheme().interpolate(C, CLambdas, CYs)
);

// solve for C using the interpolated surface values
fvScalarMatrix CAdvEqn
(
    fvm::ddt(C)
  + fvc::surfaceIntegrate(phi*tCs())
);

CAdvEqn.solve(mesh.solver("CAdv"));

volScalarField CAdv("CAdv", C);

PO2 = C/alpha;
PO2Adv = PO2;

if (min(PO2).value() < 0)
{
    WarningIn("hbPOAxisymmetricFoam::PO2AdvEqn.H")
        << "    PO2 reached a negative minimum value (" << min(PO2).value() << ")."
        << nl;
}
PO2 = pos(PO2)*PO2;

// overwrite the oldTime
PO2.oldTime() = PO2;
