#! /bin/bash -l
#$ -l h_rt=1:00:00
#$ -N matrixrun-32
#$ -j y
#$ -P paralg
#$ -pe omp 16

# Load the correct modules
module load gcc/5.3.0  # compiler
module load mpich/3.2  # consistent mpi compile

# Immediately form fused output/error file, besides the one with the default name.
exec >  ${SGE_O_WORKDIR}/${JOB_NAME}-${JOB_ID}.scc.out 2>&1

# Run the executable
for i in 1 2 4 8 16 ; do
	./Fluid32 $i
done 

exit
