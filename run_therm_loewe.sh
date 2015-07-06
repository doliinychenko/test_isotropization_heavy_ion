#!/bin/sh
#
# sbatch ./bin/loewe_run.sh
#
# number of jobs (mulitple of 24):
# single job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#
# job output:
#SBATCH --output=/scratch/hyihp/oliiny/therm_project/therm.o%j
#SBATCH --error=/scratch/hyihp/oliiny/therm_project/therm.e%j
#
# email
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliiny@fias.uni-frankfurt.de
#
# job name:
#SBATCH --job-name=therm_analysis
#
# mpi stuff:
#SBATCH --partition=parallel
#
# mem allocation (only 200m default)
#SBATCH --mem-per-cpu=4000
#
# default time 10min, max 8 days?:
#SBATCH --time=6-23:00:00

# output
Energy=$1
Centrality=$2
Cell_size='0.6'
Sigma='0.8'
Urqmd_input="/scratch/hyihp/oliiny/therm_project/Nev_dep/urqmd_E${Energy}_b${Centrality}_n*/urqmd-3.4/test.f14"
path="/scratch/hyihp/oliiny/therm_project/ev10000/E${Energy}_b${Centrality}"
src_files="/home/hyihp/oliiny/therm_project"

# output info
echo "Running on ${SLURM_NNODES} nodes."
echo "Running with ${SLURM_NTASKS} number of tasks."
echo "Executed from: ${SLURM_SUBMIT_DIR}"
echo "List of nodes: ${SLURM_JOB_NODELIST}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Saving data: $path"

# Run stuff in scratch space
mkdir -p "${path}/vtk"
mkdir -p "${path}/plots"
cp "${src_files}/thermalization" "${path}/"
cp -r "${src_files}/obj" "${path}/"
cd "${path}"


export OMP_NUM_THREADS=1

# start programm
my_n_cells=`echo 15/${Cell_size} | bc`
srun ./thermalization -load_from_saved F \
                      -urqmd_input "${Urqmd_input}" \
                      -save_load_file "Tmn_saved_E${Energy}_b${Centrality}_n${N}.bin" \
                      -nx ${my_n_cells} -nz ${my_n_cells} -nt 75 \
                      -dx ${Cell_size} -dz ${Cell_size} \
                      -sigma ${Sigma}
