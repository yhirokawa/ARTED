#! /bin/bash
#PJM -g xg17i042
#PJM -L rscgrp=regular-flat,node=32
#PJM -L elapse=00:30:00
#PJM --omp thread=256
#PJM --mpi proc=32
#PJM -j
#PJM -o job.log

export FORT_BUFFERED=1
export KMP_AFFINITY=balanced,granularity=fine
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0,1,68,69,136,137,204,205
export I_MPI_HBW_POLICY=hbw_bind,hbw_bind
export I_MPI_DEBUG=5

mpiexec.hydra ../../bin/ARTED_ms.mic < input_ms_Si.32.inp
