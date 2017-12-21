The solvers folder                    {#solverspage}
------------------

The O2RBC library contains two solvers:
  - [hbPOAxisymmetricFoam](@ref hbPOAxisymmetricFoam.C) for simulations in an axisymmetric
domain,
  - [hbPOGraphFoam](@ref hbPOGraphFoam.C) for simulations in capillary networks. 

For both solvers, there is a related executable that deals with the initial simulation
([initializeHbPOAxisymmetric](@ref initializeHbPOAxisymmetric.C) and
 [initializeHbPOGraph](@ref initializeHbPOGraph.C), respectively). Their execution is
needed to run the respective solvers.

