# traczilla
FLEXPART fork that is adapted to massive calculations of high altitude trajectories

This code has forked from FLEXPART V5 and contains a large number of modifications, some of
them which have been implemented also in more recent versions of FLEXPART

Non exhaustive list
- Execution driven by a COMMAND and a RELEASE file which are both namelists. The COMMAND file reads
a variable that determines the plan to be executed. There is a namelist file for each type of plan which
is read from RELEASE.
- TRACZILLA is organised in functional blocks (not as a file for each routine like FLEXPART)
- TRACZILLA is fully fortran 90
- TRACZILLA does not interpolate to terrain following coordinates but directly from hybrid levels
- TRACZILLA reads uses modern reanalysis: ERA-Interim, ERA5, MERRA and JRA-55
- TRACZILLA is able to use spherical harmonics coefficints as input for winds
- TRACZILLA has a large number of plans to release parcels from grids, balloon or aircraft trajectories, cloud top, ...
It can also start from externally defined file defining the position and the launching time of parcels.
- TRACZILLA is parallelized with OMP
