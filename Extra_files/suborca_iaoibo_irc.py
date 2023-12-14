#!/usr/bin/env python

import argparse
from pathlib import Path
import subprocess
import os

ORCA_version_default = '5.0.3'  # ORCA default version used for the calculation.
cpus_default = '10'  # Default number of threads (CPUs) used for the calculation.
queue_default = 'SP2'  # Default queue used for the calculation.
email = 'lucas.welington.lima@usp.br'  # email address to inform the user when the job ended.

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Create and submit ORCA jobs in SLURM queue system.
            Example of use:
            > suborca INPUT.inp -v 5.0.3 -p 10 -q TEST""")

    # Filename
    parser.add_argument('fname', help='XYZ filename.')

    # Selection of the Environment Module used
    parser.add_argument('--orca_version', '-v', choices=['5.0.4', '5.0.3', '5.0.1', '5.0.0', '4.2.1'],
                        default=ORCA_version_default,
                        help='ORCA version used for the calculation.')
    parser.add_argument('--ntasks_per_node', '-p', default=cpus_default,
                        help='Number of CPUs used by SLURM for the job.')
    parser.add_argument('--partition', '-q', choices=['SP2', 'TEST'], default=queue_default,
                        help='Number of CPUs used by SLURM for the job.')
    parser.add_argument('--MO_alpha', default=True, help='Last occupied alpha MO')
    parser.add_argument('--MO_beta', default=False, help='Last occupied beta MO')


    args = parser.parse_args()

    # Job submission directory
    directory = os.getcwd()

    # Checking if the input file exists
    inpfile = Path(args.fname)

    if inpfile.is_file() == False:
        print(
            'The file \"{}\" does not exist or there is a typo in the input name.\nPlease, give a valid input for job submission.'.format(
                inpfile))
    else:
        fname_splitted = args.fname.split('.')
        basename = fname_splitted[0]

        # ORCA version definition
        if args.orca_version == '5.0.4':
            orca_module = 'orca/5.0.4-openmpi411'
            orca_dir = '/temporario/apps/gnu/orca_5_0_4_linux_x86-64_shared_openmpi411'
        elif args.orca_version == '5.0.3':
            orca_module = 'orca/5.0.3-openmpi411'
            orca_dir = '/temporario/apps/gnu/orca_5_0_3_linux_x86-64_shared_openmpi411'
        elif args.orca_version == '5.0.1':
            orca_module = 'orca/5.0.1-openmpi411'
            orca_dir = '/scratch/apps/gnu/orca_5_0_1_linux_x86-64_shared_openmpi41a'
        elif args.orca_version == '4.2.1':
            orca_module = 'orca/4.2.1-openmpi314'
            orca_dir = '/scratch/apps/gnu/orca_4_2_1_linux_x86-64_openmpi314'

            # Number of CPUs used for the job definition
        if args.ntasks_per_node != cpus_default:
            number_cpus = args.ntasks_per_node
        else:
            number_cpus = cpus_default

            # Queue and total run time for the job definition
        if args.partition == queue_default:
            total_run_time = '192:00:00'
        elif args.partition == 'TEST':
            total_run_time = '1:00:00'

            # Arguments used for the job
        SLURM_settings1 = F"""#!/bin/bash 
#
#SBATCH --partition={args.partition}
#SBATCH -J {basename} 
#SBATCH --time={total_run_time}        
#SBATCH --mem-per-cpu=24042
#SBATCH --nodes=1
#SBATCH--ntasks={number_cpus} --ntasks-per-node={number_cpus}
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user={email}

#OpenMP settings:
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo $SLURM_JOB_ID              #ID of job allocation
echo $SLURM_SUBMIT_DIR          #Directory job where was submitted
echo $SLURM_JOB_NODELIST        #File containing allocated hostnames
echo $SLURM_NTASKS              #Total number of cores for job

export ORCA_DIR={orca_dir}
export PATH=$PATH:$ORCA_DIR

echo $PATH
echo $LD_LIBRARY_PATH

module unload openmpi
module load {orca_module}
cd {directory}
"""
        if args.MO_beta == False:
            SLURM_settings2 = F"""{orca_dir}/orca {basename}.inp > {basename}.out
{orca_dir}/orca_loc {basename}_loc.inp > {basename}_loc.out            
"""
        else:
            SLURM_settings2 = F"""{orca_dir}/orca {basename}.inp > {basename}.out
{orca_dir}/orca_loc {basename}_loc_alpha.inp > {basename}_loc_alpha.out            
{orca_dir}/orca_loc {basename}_loc_beta.inp > {basename}_loc_beta.out
"""
        SLURM_settings = SLURM_settings1 + SLURM_settings2

        # Writing submission script:
        submission_filename = basename + '.slurm'
        with open(submission_filename, 'w+') as f:
            f.write(SLURM_settings)
        print('{} submission script written.'.format(submission_filename))
        print('Job {} submitted. Job number:'.format(basename))
        process = subprocess.Popen(['sbatch {}'.format(submission_filename)], shell=True)
        process.wait()
