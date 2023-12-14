# IAO-IBO_IRC
Python script to calculate the IAO-IBO orbital localization for each IRC structure from an ORCA calculation.

_Positional arguments:_

  * fname                 XYZ file


_Optional arguments:_

  * -h, --help                   Show a help message and exit
  
  * --multip (integer)           Molecule's multiplicity
  
  * --chrg (integer)             Molecule's charge
  
  * -m, --memory (integer)       Maximum memory required per processing core (in MB)
  
  * -p, --processors (integer)   Processing cores
  
  * -n, --inpname (string)       Chosen input name
  
  * --last_MO_alpha (integer)    Last occupied alpha MO
  
  * --last_MO_beta (integer)     Last occupied beta MO

  * --MO_alpha (Boolean)         Asks for localization of alpha orbitals (in suborca_iaoibo_irc.py job submissions only)
  
  * --MO_beta (Boolean)          Asks for localization of beta orbitals (in suborca_iaoibo_irc.py job submissions only)
  
  * --flip (Boolean)             Reverses the order in which the differences are calculated

_Examples of use:_ 

  * Calculation of the IAO-IBO charges of foo0_IRC_Full_trj.xyz IRC trajectory file:

    > python3 IAO-IBO_IRC.py foo0_IRC_Full_trj.xyz -p 2 -m 2000 --chrg 1 -n foo0 --last_MO_alpha 10

    In this case, the the molecule's multiplicity is defined as 1 (singlet), since no other integer is given. Consequently, no --last_MO_beta (integer) is required.

  * Calculation of the IAO-IBO charges of foo1_IRC_Full_trj.xyz IRC trajectory file selecting the multiplicity, last beta occupied MO, and the flip option:

    > python3 IAO-IBO_IRC.py foo1_IRC_Full_trj.xyz -p 2 -m 2000 --chrg 1 --multip 2 -n foo1 --last_MO_alpha 10 --last_MO_beta 9 flip True

  The flip option was used in this case because the order in which the structures were printed in the foo1_IRC_Full_trj.xyz is the inverse that expected (i.e, the structures in the product side of the IRC were printed firstly), hence the square of the differences between charges is done in reverse order.
    
  * Use of Extra_files/IAO-IBO_IRC_dir-inp.py and Extra_files/suborca_iaoibo_irc.py scripts to create the directories and submit IAO-IBO_IRC jobs in a Slurm job scheduling system (cluster aguia4 - STI-USP):

    Firstly, create the directories and input files (single point and localization calculations):
    
    > python3 IAO-IBO_IRC_dir-inp.py foo2_IRC_Full_trj.xyz -p 10 -m 5000 --chrg 1 --multip 2 -n foo2 --last_MO_alpha 10 --last_MO_beta

    Then, secondly, run the calculations. Note that if an UKS calculation has to be performed, the flags --MO_alpha (Boolean) and --MO_beta (Boolean) should be specified. Otherwise, only one localization will be performed (using alpha orbitals as default):
    
    > python3 suborca_iaoibo_irc.py foo2_IRC_1.inp -p 10 --MO_alpha True --MO_beta True

  PS.: Note that the path for _orca_ and _orca_plot_ binaries must set both in IAO-IBO_IRC.py and suborca_iaoibo_irc.py.


