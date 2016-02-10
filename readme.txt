The cpp-version of AtomicClusters code with integrated task stealing scheduling strategy.

Task stealing is implemented via MPI passive-target remote memory access.

Two models of passive target access organisation ore offered:

1. Use external asyncronious progress engine. (Recommended for better performance, library-dependant)
- If this model is chosen, the definition of POKE_P_E in Variable.h must be commented;MPI library tuning is necessary to activate the model.
- Corresponding MPI library runtime enviroment variables should be set, e.g. MPICH_ASYNC_PROGRESS=1 MV2_ENABLE_AFFINITY=0 MV2_USE_BLOCKING=0 for mvapich2.0 or MPICH_NEMESIS_ASYNC_PROGRESS=1 MPICH_MAX_THREAD_SAFETY=multiple MPICH_GNI_USE_UNASSIGNED_CPUS=enabled for cray mpich or MPICH_ASYNC_PROGRESS=1 for IntelMPI.
For more info see corresponding MPI library documentation.

2. Use target-initiated metronomic progress engine "pokes". (No real passive target).
- Activated by default. To turn off target-initiated metronomic progress engine "pokes", comment the definition of POKE_P_E in Variable.h;


The recomended MPI libraries are: MPICH and others based on it, such as IntelMPI or MVAPICH.
!The latest implementation of Open MPI does not support the asynchronous progress!
