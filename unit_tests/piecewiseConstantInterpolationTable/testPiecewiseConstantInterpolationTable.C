/*---------------------------------------------------------------------------*\
 
Application
    testPiecewiseConstantInterpolationTable

Description
    Test the class piecewiseConstantInterpolationTable<Type>

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "piecewiseConstantInterpolationTable.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    
    {
        Info<< "Testing constructor from dictionary entry and key:" << endl;
        piecewiseConstantInterpolationTable<scalarList> table
        (
           IOdictionary
           (
               IOobject
               (
                   "dataDict",
                   runTime.constant(),
                   runTime,
                   IOobject::MUST_READ,
                   IOobject::AUTO_WRITE
               )
           ),
           "dataKey"
        );
        Info<< "Created object: " << table << endl;
    }


    Info<< "Testing constructor from file name:" << endl;
    piecewiseConstantInterpolationTable<scalarList> table
    (
        "interpolationTable"
    );
    table.outOfBounds(piecewiseConstantInterpolationTable<scalarList>::CLAMP);


    Info<< "First value: " << table[0] << endl;
    Info<< "Second value: " << table[1] << endl;

    Info<< "Value at t = 0: " << table(0) << endl;
    Info<< "Value at t = 0.5: " << table(0.5) << endl;
    Info<< "Value at t = 1: " << table(1) << endl;
    Info<< "Value at t = 1.5: " << table(1.5) << endl;
    Info<< "Value at t = 2: " << table(2) << endl;
    Info<< "Value at t = 2.5: " << table(2.5) << endl;
    Info<< "Value at t = -1: " << table(-1) << endl;

    Info<< "Change between t = 0 and t = 1: " << table.valueChanged(0, 1) << endl;
    Info<< "Change between t = 1 and t = 0: " << table.valueChanged(1, 0) << endl;
    Info<< "Change between t = 1.5 and t = 1.9: " << table.valueChanged(1.5, 1.9) << endl;
    Info<< "Change between t = 2.0 and t = 2.5: " << table.valueChanged(2.0, 2.5) << endl;

    return 0;
}




