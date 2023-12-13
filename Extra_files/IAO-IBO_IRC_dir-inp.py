##########################################################
#                     IAO-IBO_IRC.py                     #
##########################################################
# Script written by Lucas W. de Lima.                    #
# This script automates the calculation of single points #
# and localizes the occupied molecular orbitals from     #
# these calculations through the IAO-IBO method for each #
# structure of an IRC path. Then parses the electronic   #
# energy and IAO-OBO charges from the ORCA output files. #
# Run the calculation as:                                #
# $ python3 IAO-IBO_IRC.py fname -p [No. procs]          #
#     --multip [multiplicity] -m [max. memory]           #
#     -n [basename] --last_MO_alpha [No. last alpha MO]  #
#     --last_MO_beta [No. last beta MO]                  #
#                                                        #
#          Last version written on Dec 08 2023.          #
##########################################################
import os
import argparse
import shutil
import glob
import subprocess
import re
import datetime
#import numpy as np

# ORCA binaries path:
orca_dir = "/temporario/apps/gnu/orca_5_0_3_linux_x86-64_shared_openmpi411/orca"
orca_loc_dir = "/temporario/apps/gnu/orca_5_0_3_linux_x86-64_shared_openmpi411/orca_loc"

def sp_inp(maxcore, pal, chrg, multip, XYZ, filename):
    """
    Writes single point input calculations
    :param maxcore: Maximum memory used (per core) in the calculation (given in MB)
    :param pal: Number of processor used in the parallel calculation
    :param chrg: Molecular charge
    :param multip: Molecular multiplicity
    :param XYZ: Molecular cartesian coordinates
    :param filename: Input basename
    """
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
    """
    Writes the input for the orbital localization calculation
    :param basename: Basename for the GBW file
    :param iteration: Iteration of the loop in which the function is called
    :param last_MO: Number of the last occupied MO
    :param operator: Operator of the MO: 0=Alpha, 1=Beta
    :param filename: Name of the final input file
    """
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
    """
    Parses the electronic energy of the single point calculation and writes the results
    in the single_point_energies.dat file located in the parent directory
    :param file: Name of the file to be parsed
    :param parent_dir: Name of the parent directory
    """
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
        print(F"Single point energy for {basename} written at single_point_energy.dat file.\n")

def IAO_chrg_parser(loc_output, iterator, parent_dir, IAO_charges_output_name):
    """
    Parses the IAO-IBO charges from the output of the localization calculation,
    and writes the results in the proper charges output in the parent directory
    :param loc_output: Output from the localization calculation
    :param iterator: Iteration of the loop in which the function is called
    :param parent_dir: Name of the parent directory
    :param IAO_charges_output_name: Output of the IAO-IBO charges
    """
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
            print(F"Results for {os.path.basename(loc_output)} written at {os.path.basename(results_file)}.\n")
        elif iterator > 0:
            with open(results_file, "a+") as f:
                f.write(line)
            print(F"Results for {os.path.basename(loc_output)} written at {os.path.basename(results_file)}.\n")

def square_dif(array, flip=False):
    """
    Calculates the square of the difference between the charge in the Nth point of the
    IRC path in relation to the first point
    :param array: Array with the charges
    :param flip: Reverses the order in which the differences are calculated
    :return: Array with square of the differences
    """
    if flip==False:
        return [(array[i] - array[0])**2 for i in range(0, len(array))]
    else:
        return [(array[i] - array[-1])**2 for i in range(len(array)-1, -1, -1)]

def rmsd(array):
    """
    Calculates the RMSD between the charges
    :param array: Array with square of the differences
    :return: RMSD value between the charges
    """
    sum_x = np.sum(array)
    n = len(array)
    return np.sqrt((sum_x/n))

def max_square_dif(array):
    """
    Shows the maximum value of the square of the difference between the charge in
    the Nth point and the first point of the IRC path
    :param array: Array with square of the differences
    :return: Maximum value of the square of the differences
    """
    return np.max(array)

def atoms(file):
    """
    Atoms of the molecule
    :param file: Name of the output with the parsed charges for each atom
    :return: List of numbered atoms
    """
    with open(file, "r") as f:
        header = f.readline()
        splitted_line = header.split(" ")
        atoms = splitted_line[0:-1]
    return atoms

def results(file, data):
    """
    Writes the results of the maximum value of the square of the differences
    and RMSD of the charges
    :param file: Output with the charges
    :param data: Arrays with the charges
    """
    for atom in range(0, len(atoms(file))):
        Square_dif = square_dif(data[:,atom], args.flip)
        Max_Square_dif = max_square_dif(Square_dif)
        rmsd_value = rmsd(Square_dif)
        if len(atoms(file)[atom]) == 4:
            print(F"{atoms(file)[atom]}           {Max_Square_dif:.4f}          {rmsd_value:.4f}")
        else:
            print(F"{atoms(file)[atom]}            {Max_Square_dif:.4f}          {rmsd_value:.4f}")

# Parser creation
parser = argparse.ArgumentParser(description="Calculates the IAO-IBO for structures in an IRC path.")

parser.add_argument('fname', help='XYZ file')
parser.add_argument('--multip', default='1', help="Molecule's multiplicity")
parser.add_argument('--chrg', default='0', help="Molecule's charge")
parser.add_argument('-m', '--memory', default='6000', help='Maximum memory required per processing core')
parser.add_argument('-p', '--processors', default='12', help='Processing cores')
parser.add_argument('-n', '--inpname', default=None, help='Chosen input name')
parser.add_argument('--last_MO_alpha', default='None', help='Last occupied alpha MO')
parser.add_argument('--last_MO_beta', default='None', help='Last occupied beta MO')
parser.add_argument('--flip', default=False,
                    help='Reverses the order in which the differences are calculated')

args = parser.parse_args()

# Parent directory definition
parent_dir = os.getcwd()
file_path = os.path.join(parent_dir, args.fname)
filenames = glob.glob(file_path)

# Checking whether the input file(s) exists
if len(filenames) == 0:
    print(F"""The file \"{args.fname}\" does not exist or there is a typo in the input name.
    Please, give a valid XYZ file for job submission.""")

else:
    print(F"""================= Calculation info =================
IRC trajectory file:       {args.fname}
Input basename:            {args.inpname}
Molecule's charge:         {args.chrg}
Molecule's multiplicity:   {args.multip}
Max memory per core:       {args.memory}
Number of cores:           {args.processors}
Last Alpha MO:             {args.last_MO_alpha}
Last Beta MO:              {args.last_MO_beta}
Calculation starting date: {datetime.date.today()}
====================================================\n""")
    print(F"Parent directory for calculations: {parent_dir}\n")
    for filename in filenames:
        with open(filename, 'r') as irc_traj:
            irc = irc_traj.readlines()
            indices = [i for i, x in enumerate(irc) if x == irc[0]]
            print(F"{args.fname} opened.\n")
            print("**** Starting the calculation ****\n\n")

            # Creating directories according to each XYZ structure from IRC_traj file
            for i in range(0, len(indices)):
                struc_dir = os.path.join(parent_dir, F"{args.inpname}_IRC_{i}")
                os.makedirs(struc_dir)
                print(F"Directory {args.inpname}_IRC_{i} created.\n")

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
                #out_name = basename + ".out"
                inp_path = os.path.join(parent_dir, inp_name)
                sp_inp(args.memory, args.processors, args.chrg, args.multip, xyz, inp_path)
                print(F"{args.inpname}_IRC_{i}.inp input file created.\n")
                shutil.move(inp_path, struc_dir)
                os.chdir(struc_dir)
                print(F"Moving to {struc_dir} directory.\nStarting {args.inpname}_IRC_{i}.inp calculation.\n")
                # Running the calculation
                #subprocess.run(F"{orca_dir} {inp_name} > {out_name}", shell=True, check=True)
                #print(F"{args.inpname}_IRC_{i}.inp calculation done successfully.\n")
                #out_path = os.path.join(struc_dir, out_name)
                #Eel_parser(out_path, parent_dir)

                # Localization calculation
                if args.last_MO_beta == "None":
                    inp_loc_name = basename + "_loc.inp"
                    #out_loc_name = basename + "_loc.out"
                    inp_loc_path = os.path.join(struc_dir, inp_loc_name)
                    #out_loc_path = os.path.join(struc_dir, out_loc_name)
                    # Localization input creation
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path)
                    # Running localization calculation
                    #subprocess.run(F"{orca_loc_dir} {inp_loc_name} > {out_loc_name}",
                    #               shell=True, check=True)
                    # Parsing the charges from localization calculation
                    #IAO_charges_output = F"{args.inpname}_IAO_charges.dat"
                    #IAO_chrg_parser(out_loc_path, i, parent_dir, IAO_charges_output)
                else:
                    inp_loc_name_alpha = basename + "_loc_alpha.inp"
                    #out_loc_name_alpha = basename + "_loc_alpha.out"
                    inp_loc_path_alpha = os.path.join(struc_dir, inp_loc_name_alpha)
                    #out_loc_path_alpha = os.path.join(struc_dir, out_loc_name_alpha)
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path_alpha)
                    inp_loc_name_beta = basename + "_loc_beta.inp"
                    #out_loc_name_beta = basename + "_loc_beta.out"
                    inp_loc_path_beta = os.path.join(struc_dir, inp_loc_name_beta)
                    #out_loc_path_beta = os.path.join(struc_dir, out_loc_name_beta)
                    loc_inp(basename, i, args.last_MO_beta, 1, inp_loc_path_beta)
                    #subprocess.run(F"{orca_loc_dir} {basename}_loc_alpha.inp > {basename}_loc_alpha.out",
                    #               shell=True, check=True)
                    #subprocess.run(F"{orca_loc_dir} {basename}_loc_beta.inp > {basename}_loc_beta.out",
                    #               shell=True, check=True)
                    # Parsing the charges from localization calculation
                    #IAO_charges_output_alpha = F"{args.inpname}_IAO_charges_alpha.dat"
                    #IAO_charges_output_beta = F"{args.inpname}_IAO_charges_beta.dat"
                    #IAO_chrg_parser(out_loc_path_alpha, i, parent_dir, IAO_charges_output_alpha)
                    #IAO_chrg_parser(out_loc_path_beta, i, parent_dir, IAO_charges_output_beta)
                print(F"------- End of step {i+1} from {len(indices)} -------\n\n")
                os.chdir(parent_dir)

            # Charges (charge difference)² and RMSD calculation
            #if args.last_MO_beta == "None":
            #    print("Atom   max(charge difference)²   RMSD")
            #    results_file = os.path.join(parent_dir, F"{args.inpname}_IAO_charges.dat")
            #    results_data = np.loadtxt(results_file, skiprows=1)
            #    results(results_file, results_data)

            #else:
            #    results_alpha_file = os.path.join(parent_dir, F"{args.inpname}_IAO_charges_alpha.dat")
            #    results_beta_file = os.path.join(parent_dir, F"{args.inpname}_IAO_charges_beta.dat")
            #    results_alpha_data = np.loadtxt(results_alpha_file, skiprows=1)
            #    results_beta_data = np.loadtxt(results_beta_file, skiprows=1)

            #    print("       *** Alpha MO results ***\nAtom   max(charge difference)²   RMSD")
            #    results(results_alpha_file, results_alpha_data)
            #    print("\n       *** Beta MO results ***\nAtom   max(charge difference)²   RMSD")
            #    results(results_beta_file, results_beta_data)

    print("""\n               ***** HURRAY!!! *****
Directories and input files created successfully!""")
