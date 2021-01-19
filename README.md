# PEP_VLIB_GENERATOR
Automated generator of virtual library of short peptides.
The workflow of the python/bash script is represented by the following steps:
1) Generation of all the peptide sequences from the combination of the 20 aminoacids. The number of sequences wil thus be 20^n where the n is the residue number of the peptide
2) Generation of a peptide structure to be used as a template for modeller
3) Generation of one input command files and alignment for modeller for every sequence generated in the step 1)
4) Generation of the pdb models and selection of the best structure(s) on the basis of the DOPE score

REQUIREMENTS: Python, bash, modeller v9.25 (other versions should work but you have to change executable name in this script)


NOTE: this is versio 0.0, which means no options: it only generates 160000 linear tetrapeptides

I would appreciate if you could cite this project in your reserch report|paper|thesis.

