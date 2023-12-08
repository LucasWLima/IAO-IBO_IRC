# IAO-IBO_IRC
Python script to calculate the IAO-IBO orbital localization for each IRC structure from an ORCA calculation.

Positional arguments:
  fname                 XYZ file

Optional arguments:
  -h, --help                   Show a help message and exit
  --multip (integer)           Molecule's multiplicity
  --chrg (integer)             Molecule's charge
  -m, --memory (integer)       Maximum memory required per processing core (in MB)
  -p, --processors (integer)   Processing cores
  -n, --inpname (string)       Chosen input name
  --last_MO_alpha (integer)    Last occupied alpha MO
  --last_MO_beta (integer)     Last occupied beta MO
  --flip (Boolean)             Reverses the order in which the differences are calculated


