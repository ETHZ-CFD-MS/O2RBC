The meshTools library                   {#meshToolspage}
---------------------

The meshTools library contains classes related to mesh-to-mesh interpolation and
mesh-related geometric information.

The most important class in this library is 
[meshToMeshMoving](@ref Foam::meshToMeshMoving) which is an adapted version of
the OpenFOAM class meshToMesh to moving meshes. Within the class folder, a
custom version of the intersection-volume-based interpolation scheme is implemented
([cellRelativeVolumeWeightMethod](@ref Foam::cellRelativeVolumeWeightMethod)).
The class [meshToSubMesh](@ref Foam::meshToSubMesh) supports the copy of fields
between a mesh and a submesh composed of cells of the original mesh.

RBC motion requires a description of the RBC position and shape. This is
performed using the abstract class 
[deformableBodyGeometricState](@ref Foam::deformableBodyGeometricState) and its
concrete implementation 
[axisymmetricBodyGeometricState](@ref Foam::axisymmetricBodyGeometricState) for
axisymmetric bodies.

