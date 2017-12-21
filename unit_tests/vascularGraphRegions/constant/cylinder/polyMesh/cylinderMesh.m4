/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])


    meshGenApp blockMesh;
    convertToMeters 1e-0;

    define(D, 4)  // cylinder diameter
    define(L, 20) // cylinder length
    define(PI, 3.14159265)
    
    define(RO,  calc(0.5*D)) // cylinder outer radius
    define(RI,  calc(0.8*RO)) // cylinder inner radius
    define(CW,  calc(0.5*RI)) // width of middle square section
    
    define(CYI, calc(RI*cos((PI/180)*45)))
    define(CZI, calc(RI*sin((PI/180)*45)))
    define(CYO, calc(RO*cos((PI/180)*45)))
    define(CZO, calc(RO*sin((PI/180)*45)))
    
    define(NS,   5) // number of cells in the square section
    define(NPI,  1) // number of cells from square section to inner perimeter
    define(NPO,  1) // number of cells from inner to outer perimeter
    define(NX, 60) // number of cells from left to right

vertices
(
    (0.0  CW  CW) vlabel(oneoclocksql)
    (0.0 -CW  CW) vlabel(fiveoclocksql)
    (0.0 -CW -CW) vlabel(sevenoclocksql)
    (0.0  CW -CW) vlabel(elevenoclocksql)
   
    (0.0  CYI  CZI) vlabel(oneoclockcil)
    (0.0 -CYI  CZI) vlabel(fiveoclockcil)
    (0.0 -CYI -CZI) vlabel(sevenoclockcil)
    (0.0  CYI -CZI) vlabel(elevenoclockcil)

    (0.0  CYO  CZO) vlabel(oneoclockcol)
    (0.0 -CYO  CZO) vlabel(fiveoclockcol)
    (0.0 -CYO -CZO) vlabel(sevenoclockcol)
    (0.0  CYO -CZO) vlabel(elevenoclockcol)

    (L  CW  CW) vlabel(oneoclocksqr)
    (L -CW  CW) vlabel(fiveoclocksqr)
    (L -CW -CW) vlabel(sevenoclocksqr)
    (L  CW -CW) vlabel(elevenoclocksqr)
   
    (L  CYI  CZI) vlabel(oneoclockcir)
    (L -CYI  CZI) vlabel(fiveoclockcir)
    (L -CYI -CZI) vlabel(sevenoclockcir)
    (L  CYI -CZI) vlabel(elevenoclockcir)

    (L  CYO  CZO) vlabel(oneoclockcor)
    (L -CYO  CZO) vlabel(fiveoclockcor)
    (L -CYO -CZO) vlabel(sevenoclockcor)
    (L  CYO -CZO) vlabel(elevenoclockcor)
);				

blocks
(
    // square block
    hex (
        sevenoclocksql sevenoclocksqr elevenoclocksqr elevenoclocksql
        fiveoclocksql  fiveoclocksqr  oneoclocksqr    oneoclocksql
    )
    (NX NS NS)
    simpleGrading (1 1 1)

    // three o'clock inner slice
    hex (
        fiveoclocksql fiveoclocksqr oneoclocksqr oneoclocksql
        fiveoclockcil fiveoclockcir oneoclockcir oneoclockcil
    )
    (NX NS NPI)
    simpleGrading (1 1 1)

    // three o'clock outer slice
    hex (
        fiveoclockcil fiveoclockcir oneoclockcir oneoclockcil
        fiveoclockcol fiveoclockcor oneoclockcor oneoclockcol
    )
    (NX NS NPO)
    simpleGrading (1 1 1)

    // six o'clock inner slice
    hex (
        sevenoclockcil sevenoclockcir sevenoclocksqr sevenoclocksql
        fiveoclockcil  fiveoclockcir  fiveoclocksqr  fiveoclocksql
    )
    (NX NPI NS)
    simpleGrading (1 1 1)

    // six o'clock outer slice
    hex (
        sevenoclockcol sevenoclockcor sevenoclockcir sevenoclockcil 
        fiveoclockcol  fiveoclockcor  fiveoclockcir  fiveoclockcil 
    )
    (NX NPO NS)
    simpleGrading (1 1 1)

    // nine o'clock inner slice
    hex (
        sevenoclockcil sevenoclockcir elevenoclockcir elevenoclockcil
        sevenoclocksql sevenoclocksqr elevenoclocksqr elevenoclocksql
    )
    (NX NS NPI)
    simpleGrading (1 1 1)

    // nine o'clock outer slice
    hex (
        sevenoclockcol sevenoclockcor elevenoclockcor elevenoclockcol
        sevenoclockcil sevenoclockcir elevenoclockcir elevenoclockcil
    )
    (NX NS NPO)
    simpleGrading (1 1 1)

    // twelve o'clock inner slice
    hex (
        elevenoclocksql elevenoclocksqr elevenoclockcir elevenoclockcil
        oneoclocksql    oneoclocksqr    oneoclockcir    oneoclockcil
       )
    (NX NPI  NS)
    simpleGrading (1 1 1)

    // twelve o'clock outer slice
    hex (
        elevenoclockcil elevenoclockcir elevenoclockcor elevenoclockcol
        oneoclockcil    oneoclockcir    oneoclockcor    oneoclockcol
       )
    (NX NPO  NS)
    simpleGrading (1 1 1)

);


// create the quarter circles
edges
(
    // left, inner arcs
    arc oneoclockcil    fiveoclockcil   (0.0  0.0   RI)
    arc fiveoclockcil   sevenoclockcil  (0.0  -RI  0.0)
    arc sevenoclockcil  elevenoclockcil (0.0  0.0  -RI)
    arc elevenoclockcil oneoclockcil    (0.0   RI  0.0)

    // left, outer arcs
    arc oneoclockcol    fiveoclockcol   (0.0  0.0   RO)
    arc fiveoclockcol   sevenoclockcol  (0.0  -RO  0.0)
    arc sevenoclockcol  elevenoclockcol (0.0  0.0  -RO)
    arc elevenoclockcol oneoclockcol    (0.0   RO  0.0)

    // right, inner arcs
    arc oneoclockcir    fiveoclockcir   (L  0.0   RI)
    arc fiveoclockcir   sevenoclockcir  (L  -RI  0.0)
    arc sevenoclockcir  elevenoclockcir (L  0.0  -RI)
    arc elevenoclockcir oneoclockcir    (L   RI  0.0)

    // right, outer arcs
    arc oneoclockcor    fiveoclockcor   (L  0.0   RO)
    arc fiveoclockcor   sevenoclockcor  (L  -RO  0.0)
    arc sevenoclockcor  elevenoclockcor (L  0.0  -RO)
    arc elevenoclockcor oneoclockcor    (L   RO  0.0)
);

patches
(
    patch left
    (
        (fiveoclocksql oneoclocksql elevenoclocksql sevenoclocksql)

        (fiveoclocksql fiveoclockcil oneoclockcil oneoclocksql)
        (fiveoclockcil fiveoclockcol oneoclockcol oneoclockcil)

        (sevenoclocksql sevenoclockcil fiveoclockcil fiveoclocksql)
        (sevenoclockcil sevenoclockcol fiveoclockcol fiveoclockcil)

        (elevenoclocksql elevenoclockcil sevenoclockcil sevenoclocksql)
        (elevenoclockcil elevenoclockcol sevenoclockcol sevenoclockcil)

        (oneoclocksql oneoclockcil elevenoclockcil elevenoclocksql)
        (oneoclockcil oneoclockcol elevenoclockcol elevenoclockcil)
    )

    patch right
    (
        (oneoclocksqr fiveoclocksqr sevenoclocksqr elevenoclocksqr)

        (oneoclocksqr oneoclockcir fiveoclockcir fiveoclocksqr)
        (oneoclockcir oneoclockcor fiveoclockcor fiveoclockcir)

        (fiveoclocksqr fiveoclockcir sevenoclockcir sevenoclocksqr)
        (fiveoclockcir fiveoclockcor sevenoclockcor sevenoclockcir)

        (sevenoclocksqr sevenoclockcir elevenoclockcir elevenoclocksqr)
        (sevenoclockcir sevenoclockcor elevenoclockcor elevenoclockcir)

        (elevenoclocksqr elevenoclockcir oneoclockcir oneoclocksqr)
        (elevenoclockcir elevenoclockcor oneoclockcor oneoclockcir)
    )

    wall walls
    (
        (fiveoclockcol   fiveoclockcor   oneoclockcor    oneoclockcol)
        (oneoclockcol    oneoclockcor    elevenoclockcor elevenoclockcol)
        (elevenoclockcol elevenoclockcor sevenoclockcor  sevenoclockcol)
        (sevenoclockcol  sevenoclockcor  fiveoclockcor   fiveoclockcol)
    )

);
