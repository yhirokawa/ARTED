# Requirement

We use Intel compiler, MPI, MKL library and CMake build tool

    $ module list
    Currently Loaded Modulefiles:
      1) impi/2017.3.196    2) intel/2017.4.196
    $ which cmake
    /usr/bin/cmake
    $ cmake --version
    cmake version 2.8.11


# Build

    $ cd ${ARTED_DOWNLOAD_DIRECTORY}
    $ mkdir build-temp
    $ cd build-temp
    $ ../configure.py -t ms --arch=intel-knl --disable-swp
    $ make -j4
    $ ls ../bin
    ARTED_ms.mic


# Execution

test/ directory has a 8--32 node execution test jobs.

    $ cd ${ARTED_DOWNLOAD_DIRECTORY}
    $ cd test
    $ ls
    n08 n16 n32
    $ cd n16
    $ pjsub run.sh


# Get timer log

All process writes timer log to file.

hpsi means Hamiltonian which is dominant computation in time development part.

    $ cat timer_YYYYMMDD_hhmmss_p?????.log | grep hpsi


# Cordinate computation time

At line 31 in input file (input_ms_Si.??.inp), ``Nt'' means # of time development part iteration counts.

A default value is 2000.
