##########################################################
#                     IAO-IBO_IRC.py                     #
##########################################################
# Script written by Lucas W. de Lima.                    #
# This script creates directories and input files for    #
# single point calculations and localizes the occupied   #
# molecular orbitals from these calculations through     #
# the IAO-IBO method for each structure of an IRC path.  #
# Run the calculation as:                                #
# $ python3 IAO-IBO_IRC_dir-inp.py fname -p [No. procs]  #
#     --multip [multiplicity] -m [max. memory]           #
#     -n [basename] --last_MO_alpha [No. last alpha MO]  #
#     --last_MO_beta [No. last beta MO]                  #
#                                                        #
#          Last version written on Dec 13 2023.          #
##########################################################
import os
import argparse
import shutil
import glob
import datetime

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
                inp_path = os.path.join(parent_dir, inp_name)
                sp_inp(args.memory, args.processors, args.chrg, args.multip, xyz, inp_path)
                print(F"{args.inpname}_IRC_{i}.inp input file created.\n")
                shutil.move(inp_path, struc_dir)
                os.chdir(struc_dir)
                print(F"Moving to {struc_dir} directory.\nStarting {args.inpname}_IRC_{i}.inp calculation.\n")

                # Localization calculation
                if args.last_MO_beta == "None":
                    inp_loc_name = basename + "_loc.inp"
                    inp_loc_path = os.path.join(struc_dir, inp_loc_name)
                    # Localization input creation
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path)
                else:
                    inp_loc_name_alpha = basename + "_loc_alpha.inp"
                    inp_loc_path_alpha = os.path.join(struc_dir, inp_loc_name_alpha)
                    loc_inp(basename, i, args.last_MO_alpha, 0, inp_loc_path_alpha)
                    inp_loc_name_beta = basename + "_loc_beta.inp"
                    inp_loc_path_beta = os.path.join(struc_dir, inp_loc_name_beta)
                    loc_inp(basename, i, args.last_MO_beta, 1, inp_loc_path_beta)
                print(F"------- End of step {i+1} from {len(indices)} -------\n\n")
                os.chdir(parent_dir)

    print("""\n             ***** HURRAY!!! *****
Directories and input files created successfully!""")
