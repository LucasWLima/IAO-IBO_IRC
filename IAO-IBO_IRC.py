import os
import argparse
import shutil
import glob
import sys
import subprocess
from datetime import datetime

def sp_temp(maxcore, pal, chrg, multip, XYZ, GBW):
    if multip == '1':
        method = ['! RKS TPSS0 D3BJ def2-TZVPP def2/J RIJCOSX TightSCF\n']
    else:
        method = ['! UKS TPSS0 D3BJ def2-TZVPP def2/J RIJCOSX TightSCF\n']

    if GBW == 'None':
        input_settings = ['! NormalPrint PrintBasis PrintMOs\n\n',
                          '%scf\n',
                          'MaxIter  500\n',
                          'end\n\n']
    else:
        gbw_name = F'MOInp {GBW}\n'
        input_settings = ['! NormalPrint PrintBasis PrintMOs\n\n',
                          '%scf\n',
                          'Guess MORead\n',
                          gbw_name,
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
    return text

# Parser creation
parser = argparse.ArgumentParser(description="Calculates the IAO-IBO for structures in a IRC path.")

parser.add_argument('fname', help='XYZ file')
parser.add_argument('--multip', default='1', help="Molecule's multiplicity")
parser.add_argument('--chrg', default='0', help="Molecule's charge")
parser.add_argument('-m', '--memory', default='6000', help='Maximum memory required per processing core')
parser.add_argument('-p', '--processors', default='12', help='Processing cores')
parser.add_argument('-n', '--inpname', default=None, help='Chosen input name')
parser.add_argument('--gbw', default='None', help='GBW file')

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
                inp = sp_temp(args.memory, args.processors, args.chrg, args.multip, xyz, args.gbw)
                inp_name = F"{args.inpname}_IRC_{i}.inp"
                inp_path = os.path.join(parent_dir, inp_name)
                with open(inp_path, "w+") as f:
                    for line in inp:
                        f.write(line)
                print(F"{inp_name} created")
                shutil.move(inp_path, struc_dir)
                os.chdir(struc_dir)
                print(F"I am in {os.getcwd()}")
                subprocess.run(F"vim {inp_name}", shell=True)
                os.chdir(parent_dir)
                print(F"Now I am in {os.getcwd()}")