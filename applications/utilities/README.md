The utilities folder                    {#utilitiespage}
------------------

Various utilities used for pre- or postprocessing of oxygen transport
simulations. The source files are individually documented.

Some utilities are called by each Allrun script in the graph-based simulations.
The utility [createVascularGraphRegions](@ref createVascularGraphRegions.C) 
creates fields that are indicators for the various vascular regions (capillary
lumen, capillary wall and tissue). The RBC sample mesh is scaled by
[adjustRBCMeshVolume](@ref adjustRBCMeshVolume.C) to have exactly the required
volume. Finally, [mySplitMeshRegions](@ref mySplitMeshRegions.C) corrects a bug
in splitMeshRegions that occurs in parallel simulations.

The tissue volume closest to each vessel is computed by 
[computeTissueGraphRegions](@ref computeTissueGraphRegions.C). Its advantage is
that the computation exactly takes into account the volume occupied by the blood
vessels themselves.

The utilities [concatenateEdgeVelocities](@ref concatenateEdgeVelocities.C) and
[concatenateRBCPaths](@ref concatenateRBCPaths.C) were written to concatenate
RBC paths files when it was not possible to realize this in Python due to 
size of the pickle files usually employed.



