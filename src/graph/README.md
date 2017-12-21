The graph library                   {#graphpage}
--------------------

The graph library contains all graph-related classes, from a hierarchy of
classes for graphs to position interpolation and definition of vascular regions
based on graph topology. This library is used for simulations of oxygen
transport in capillary networks. The simulations in axisymmetric domains do not
require this library.

A class hierarchy with various levels of abstractions describes the vascular
graph structure. The hierarchy is best visualized in the documentation of
[pointGraph](@ref Foam::pointGraph). The class [IOgraph](@ref Foam::IOgraph)
provides input/output functionality for graphs defined by an adjacency list. 
The graph structure itself is held by [primitiveGraph](@ref Foam::primitiveGraph).
The library [boost](http://www.boost.org/) is used for all graph queries. 
This class defines wrappers for
boost function, for instance for adjacent vertices/edges and degree computation.
The subclass [pointGraph](@ref Foam::pointGraph) describes a graph embedded in
three-dimensional space. The next subclass, 
[geometricEdgeGraph](@ref Foam::geometricEdgeGraph), enables computations related 
to the edge geometry described by the abstract class 
[geometricEdge](@ref Foam::geometricEdge). It has concrete implementations for
straight edges ([straightEdge](@ref Foam::straightEdge)) and polygonal edges
([polygonalEdge](@ref Foam::polygonalEdge)).
The last graph class, [circularTubeGraph](@ref Foam::circularTubeGraph), is
actually used by the solver [hbPOGraphFoam](@ref hbPOGraphFoam.C). A instance of
[circularTube](@ref Foam::circularTube) corresponds to each edge. The geometric
tube information is used to construct capillary networks and define the
associated regions (plasma, capillary wall, tissue).

The computation of these regions is done in the class
[vascularGraphRegions](@ref Foam::vascularGraphRegions). Similarly, the class
[tissueGraphRegions](@ref Foam::tissueGraphRegions) partitions the tissue
based on the topological distance to the nearest edge and can compute field
averages on these regions.

Interpolation of time-dependent series of positions on a graph is implemented in
[graphCoordinateInterpolation](@ref Foam::graphCoordinateInterpolation). The
creation of cell sets based on the distance to the nearest edge is done in
[geometricEdgeGraphToCell](@ref Foam::geometricEdgeGraphToCell). This is used to
define a mesh region where the mesh should be refined.


