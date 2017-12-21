The postProcessing library                    {#postProcessingpage}
-------------------------

This library contains custom function objects for the postprocessing
of RBC data.

The function object [RBCProbe](@ref Foam::RBCProbe) is the equivalent
of OpenFOAM's probe function object for RBCs passing through a graph.
For the graph version of the code, the function object 
[sampleRBCField](@ref Foam::sampleRBCField) samples the values of a 
field in the RBCs that are currently in the computational domain. 
The corresponding code for axisymmetric simulations is implemented
in [sampleRBCFieldAxisymmetric](@ref Foam::sampleRBCFieldAxisymmetric).


