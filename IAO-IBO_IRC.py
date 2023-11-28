import os
import argparse
import shutil
import glob
import sys
import subprocess
from datetime import datetime

# ORCA binaries path:
orca_dir = ""
orca_loc_dir = ""

def sp_inp(maxcore, pal, chrg, multip, XYZ, filename):
    if multip == '1':
        method = ['! RKS TPSS0 D3BJ def2-TZVPP def2/J RIJCOSX TightSCF\n']
    else:
        method = ['! UKS TPSS0 D3BJ def2-TZVPP def2/J RIJCOSX TightSCF\n']

    input_settings = ['! NormalPrint PrintBasis PrintMOs\n\n',
                          '%scf\n',
                          'MaxIter  500\n',
                          'end\n\n']
    memory = list(F'%maxcore {maxcore}\n\n')
    procs = F'nprocs {pal}\n'
    nprocs = ['%pal\n', procs, 'end\n\n']
    cpcm = ['%cpcm\n', 'epsilon  2.82  # DMC Dieletric constant\n', 'end\n\n']
    charge_multiplicity = list(F'* xyz {chrg} {multip}\n')
    xyz_coordinates = XYZ
    tail = ['*']
    text = method + input_settings + memory + nprocs + cpcm + charge_multiplicity + xyz_coordinates + tail

    with open(filename, "w+") as f:
        for line in text:
            f.write(line)
def loc_inp(basename, iteration, last_MO, operator, filename):
    loc_inp_text = F"""{basename}_IRC_{iteration}.gbw          # input orbitals
{basename}_IRC_{iteration}.loc.gbw      # output orbitals
0                      # orbital window: first orbital to be localized e.g. first active
{last_MO}                     # orbital window: last orbital to be localized e.g. last active
3                      # localization method: # 1=PIPEK-MEZEY,2=FOSTER-BOYS,3=IAO-IBO,4=IAO-BOYS,5=NEW-BOYS,6=AHFB
{operator}                      # operator: 0 for alpha, 1 for beta
128                    # maximum number of iterations
1e-6                   # convergence tolerance of the localization functional value
0.0                    # relative convergence tolerance of the localization functional value
0.95                   # printing thresh to call an orbital strongly localized
0.85                   # printing thresh to call an orbital bond-like
2                      # printlevel
1                      # use Cholesky Decomposition (0=false, 1=true)
1                      # randomize seed for localization (0=false, 1=true)
"""
    with open(filename, "w+") as f:
        f.write(loc_inp_text)

# Parser creation
parser = argparse.ArgumentParser(description="Calculates the IAO-IBO for structures in a IRC path.")

parser.add_argument('fname', help='XYZ file')
parser.add_argument('--multip', default='1', help="Molecule's multiplicity")
parser.add_argument('--chrg', default='0', help="Molecule's charge")
parser.add_argument('-m', '--memory', default='6000', help='Maximum memory required per processing core')
parser.add_argument('-p', '--processors', default='12', help='Processing cores')
parser.add_argument('-n', '--inpname', default=None, help='Chosen input name')
parser.add_argument('--last_MO_alpha', default='None', help='Last occupied alpha MO')
parser.add_argument('--last_MO_beta', default='None', help='Last occupied beta MO')

args = parser.parse_args()

parent_dir = os.getcwd()
file_path = os.path.join(parent_dir, args.fname)
filenames = glob.glob(file_path)

# Checking whether the input file(s) exists
if len(filenames)==0:
    print(F'The file \"{args.fname}\" does not exist or there is a typo in the input name.\nPlease, give a valid input for job submission.')

else:
    for filename in filenames:
        with open(filename, 'r') as irc_traj:
            irc = irc_traj.readlines()
            indices = [i for i, x in enumerate(irc) if x == irc[0]]

            for i in range(0, len(indices)):
                struc_dir = os.path.join(parent_dir, F"{args.inpname}_IRC_{i}")
                os.makedirs(struc_dir)

                if indices[i] < indices[-1]:
                    struc_coordinates = irc[indices[i]:indices[i + 1]]
                elif indices[i] == indices[-1]:
                    struc_coordinates = irc[indices[-1]:]

                xyz = struc_coordinates[2:]
                basename = F"{args.inpname}_IRC_{i}"
                inp_name = basename + ".inp"
                inp_path = os.path.join(parent_dir, inp_name)
                sp_inp(args.memory, args.processors, args.chrg, args.multip, xyz, inp_path)
                shutil.move(inp_path, struc_dir)
                os.chdir(struc_dir)
                subprocess.run(F"{orca_dir} {basename}.inp > {basename}.out", shell=True)
                # Escrever função para parsing das Eel e colocar código aqui
                if args.last_MO_beta == "None":
                    inp_loc_name = basename + "_loc.inp"
                    inp_loc_path = os.path.join(struc_dir, inp_loc_name)
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path)
                    subprocess.run(F"{orca_loc_dir} {basename}_loc.inp > {basename}_loc.out", shell=True)
                    # add parser here
                else:
                    inp_loc_name_alpha = basename + "_loc_alpha.inp"
                    inp_loc_path_alpha = os.path.join(struc_dir, inp_loc_name_alpha)
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path_alpha)
                    inp_loc_name_beta = basename + "_loc_beta.inp"
                    inp_loc_path_beta = os.path.join(struc_dir, inp_loc_name_beta)
                    loc_inp(basename, i, args.last_MO_beta, 1, inp_loc_path_beta)
                    subprocess.run(F"{orca_loc_dir} {basename}_loc_alpha.inp > {basename}_loc_alpha.out", shell=True)
                    subprocess.run(F"{orca_loc_dir} {basename}_loc_beta.inp > {basename}_loc_beta.out", shell=True)
                    # add parser here
                os.chdir(parent_dir)
                print(F"Now I am in {os.getcwd()}")
