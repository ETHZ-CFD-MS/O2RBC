The topoChangerFvMesh library                    {#topoChangerFvMeshpage}
-----------------------------

This library contains classes derived from the OpenFOAM class
topoChangerFvMesh. Their goal is to provide a convenient interface for
a single mesh object with multiple disconnected regions (that correspond
to individual RBCs).

The class [regionAddRemoveFvMesh](@ref Foam::regionAddRemoveFvMesh) implements
the topological addition of a mesh to the current mesh, as well as the removal
of regions defined by cell zones. The derived class 
[disconnectedZoneMesh](@ref Foam::disconnectedZoneMesh) supports the motion of
disconnected cell zones, as well as field setting and getting on individual
zones.


