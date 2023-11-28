import os
import argparse
import shutil
import glob
#import sys
import subprocess
import re
#from datetime import datetime

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
    loc_inp_text = F"""{basename}.gbw          # input orbitals
{basename}.loc.gbw      # output orbitals
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

def Eel_parser(file, parent_dir):
    with open(file, "r") as f:
        data = f.read()
        searchbox = re.findall(r"FINAL SINGLE POINT ENERGY(.*)\n", data)
        Electronic_energy = searchbox[0].strip()
        basename = os.path.basename(file).split(".")
        basename = basename[0]
        # Writing results into single_point_energies.dat file
        sp_energies_file = os.path.join(parent_dir, "single_point_energies.dat")
        with open(sp_energies_file, "a+") as f:
            f.write(F"{basename} {Electronic_energy}\n")
        print(F"Single point energy for {basename} written at single_point_energy.dat file.")

def IAO_chrg_parser(loc_output, iterator, parent_dir, IAO_charges_output_name):
    with open(loc_output, 'r') as f:
        data = f.read()
        searchbox = re.findall(r"^Warning!!!(.*)^Initial", data, re.MULTILINE | re.DOTALL)
        searchbox = searchbox[0].split("\n")
        charge_lines = searchbox[1:-3]
        sum_charges = searchbox[-3].split(" ")
        sum_charges = sum_charges[-1]
        splitted_lines = [line.split(":") for line in charge_lines]
        atoms = []
        charges = []
        for lst in range(0, len(splitted_lines)):
            charge = splitted_lines[lst][1]
            charge = charge.strip()
            charges.append(charge)
        for lst in range(0, len(splitted_lines)):
            atom = splitted_lines[lst][0]
            atom = atom.strip()
            atom = atom.replace(" ", "-")
            atoms.append(atom)
        charges_line = ""
        for i in range(0, len(charges)):
            charges_line += F"{charges[i]} "
        atoms_line = ""
        for i in range(0, len(atoms)):
            atoms_line += F"{atoms[i]} "
        header = atoms_line + "Sum_charges\n"
        line = charges_line + F"{sum_charges}\n"
        results_file = os.path.join(parent_dir, IAO_charges_output_name)
        # Writing the results
        if iterator == 0:
            with open(results_file, "a+") as f:
                f.write(header)
            with open(results_file, "a+") as f:
                f.write(line)
            print(F"Results for {os.path.basename(loc_output)} written at {os.path.basename(results_file)}.")
        elif iterator > 0:
            with open(results_file, "a+") as f:
                f.write(line)
            print(F"Results for {os.path.basename(loc_output)} written at {os.path.basename(results_file)}.")

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

# Parent directory definition
parent_dir = os.getcwd()
file_path = os.path.join(parent_dir, args.fname)
filenames = glob.glob(file_path)
print(F"Parent directory for calculations: {parent_dir}")

# Checking whether the input file(s) exists
if len(filenames) == 0:
    print(F"""The file \"{args.fname}\" does not exist or there is a typo in the input name.
    Please, give a valid XYZ file for job submission.""")

else:
    for filename in filenames:
        with open(filename, 'r') as irc_traj:
            irc = irc_traj.readlines()
            indices = [i for i, x in enumerate(irc) if x == irc[0]]
            print(F"{args.fname} opened.")

            # Creating directories according to each XYZ structure from IRC_traj file
            for i in range(0, len(indices)):
                struc_dir = os.path.join(parent_dir, F"{args.inpname}_IRC_{i}")
                os.makedirs(struc_dir)
                print(F"Directory {args.inpname}_IRC_{i} created.")

                # From first up to N-1 IRC_traj structures parsing
                if indices[i] < indices[-1]:
                    struc_coordinates = irc[indices[i]:indices[i + 1]]
                # Last IRC_traj structure parsing
                elif indices[i] == indices[-1]:
                    struc_coordinates = irc[indices[-1]:]

                # input creation for the parsed structure
                xyz = struc_coordinates[2:]
                basename = F"{args.inpname}_IRC_{i}"
                inp_name = basename + ".inp"
                out_name = basename + ".out"
                inp_path = os.path.join(parent_dir, inp_name)
                sp_inp(args.memory, args.processors, args.chrg, args.multip, xyz, inp_path)
                print(F"{args.inpname}_IRC_{i}.inp input file created.")
                shutil.move(inp_path, struc_dir)
                os.chdir(struc_dir)
                print(F"Moving to {struc_dir} directory.\nStarting {args.inpname}_IRC_{i}.inp calculation.")
                # Running the calculation
                subprocess.run(F"{orca_dir} {inp_name} > {out_name}", shell=True, check=True)
                print(F"{args.inpname}_IRC_{i}.inp calculation done successfully.")
                out_path = os.path.join(struc_dir, out_name)
                Eel_parser(out_path, parent_dir)

                # Localization calculation
                if args.last_MO_beta == "None":
                    inp_loc_name = basename + "_loc.inp"
                    out_loc_name = basename + "_loc.out"
                    inp_loc_path = os.path.join(struc_dir, inp_loc_name)
                    out_loc_path = os.path.join(struc_dir, out_loc_name)
                    # Localization input creation
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path)
                    # Running localization calculation
                    subprocess.run(F"{orca_loc_dir} {inp_loc_name} > {out_loc_name}",
                                   shell=True, check=True)
                    # Parsing the charges from localization calculation
                    IAO_charges_output = F"{args.inpname}_IAO_charges.dat"
                    IAO_chrg_parser(out_loc_path, i, parent_dir, IAO_charges_output)
                else:
                    inp_loc_name_alpha = basename + "_loc_alpha.inp"
                    out_loc_name_alpha = basename + "_loc_alpha.out"
                    inp_loc_path_alpha = os.path.join(struc_dir, inp_loc_name_alpha)
                    out_loc_path_alpha = os.path.join(struc_dir, out_loc_name_alpha)
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path_alpha)
                    inp_loc_name_beta = basename + "_loc_beta.inp"
                    out_loc_name_beta = basename + "_loc_beta.out"
                    inp_loc_path_beta = os.path.join(struc_dir, inp_loc_name_beta)
                    out_loc_path_beta = os.path.join(struc_dir, out_loc_name_beta)
                    loc_inp(basename, i, args.last_MO_beta, 1, inp_loc_path_beta)
                    subprocess.run(F"{orca_loc_dir} {basename}_loc_alpha.inp > {basename}_loc_alpha.out",
                                   shell=True, check=True)
                    subprocess.run(F"{orca_loc_dir} {basename}_loc_beta.inp > {basename}_loc_beta.out",
                                   shell=True, check=True)
                    # Parsing the charges from localization calculation
                    IAO_charges_output_alpha = F"{args.inpname}_IAO_charges_alpha.dat"
                    IAO_charges_output_beta = F"{args.inpname}_IAO_charges_beta.dat"
                    IAO_chrg_parser(out_loc_path_alpha, i, parent_dir, IAO_charges_output_alpha)
                    IAO_chrg_parser(out_loc_path_beta, i, parent_dir, IAO_charges_output_beta)
                print(F"------- End of step {i+1} from {len(indices)} -------\n\n")
                os.chdir(parent_dir)
            print(F"Done!")
