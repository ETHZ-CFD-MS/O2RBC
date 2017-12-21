Tutorials                                         {#tutorialspage}
=========

All tutorials are equipped with a script named `Allrun` that runs preprocessing,
the solver itself, as well as postprocessing and plotting in some cases.
Thus, they can be run in a terminal by navigating to the desired directory and
running

    ./Allrun

Structure
---------

The tutorials are split in four folders with different types of simulations.
    - axisymmetric: two-dimensional axisymmetric domain (cone
            or cylinder)
    - parallelCapillaries: arrays of parallel straight capillaries
    - artificialNetworks: artificial networks such as the "loop" topology or a
    converging bifurcation
    - realistic: microvascular networks acquired from animals

Additionally, some tutorials in the folder "axisymmetric" setup and run
parameter studies by using the capabilities of the Python library for simulation
pre- and postprocessing.

The tutorials in the folder "realistic" contain a script `postProcessCase` that
performs extensive postprocessing and plotting. To dive into simulation
postprocessing, it is recommended to use a top-down approach by starting from
this stript.

Running on a cluster with LSF
-----------------------------

O2RBC has been used on a cluster using the IBM LSF (load sharing facility).
The [introduction](https://scicomp.ethz.ch/wiki/Getting_started_with_clusters#Using_the_batch_system) 
on the ETH Zurich website is useful, as well as the
[LSF submission line advisor](https://scicomp.ethz.ch/lsf_submission_line_advisor).
The [documentation](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_welcome.html)
is available on the IBM website.

For computationally expensive cases, it is useful to use job chaining. 
In tutorials equipped with a script `restartSimulation`, a simulation can be restarted 
from the latest written time step after hitting the job time limit as follows.

The first job in the chain should be submitted for example as follows:

    bsub -n 24 -R "rusage[mem=8192]" -W 4:00 -J <jobName> -o lsf.out < Allrun

and the following jobs as 

    bsub -n 24 -R "rusage[mem=8192]" -W 4:00 -J <jobName> -w "exit('<jobName>', 140)" -o lsf.out < restartSimulation

The jobs submitted with the latter syntax will start only if the previous job
with name <jobName> has been killed by LSF upon reaching its time limit, which
makes it exit with code 140.

