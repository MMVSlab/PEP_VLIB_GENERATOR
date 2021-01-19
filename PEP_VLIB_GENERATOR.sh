#!/bin/bash
# Version 0.0
# 
#funzione di parallelizzazione dei processi.
function forky () {
   local num_par_procs
    if [[ -z $1 ]] ; then
       num_par_procs=8
    else
        num_par_procs=$1
    fi
     while [[ $(jobs | wc -l) -ge $num_par_procs ]] ; do
        sleep 1
     done
 }

#GENERATES A PDB TRIPEPTIDE TEMPLATE
echo "HEADER
ATOM    676  N   ALA A  1       52.566  13.216  18.100  1.00  7.47           N  
ATOM    677  CA  ALA A  1       51.848  13.797  16.973  1.00  8.15           C  
ATOM    678  C   ALA A  1       50.394  13.860  17.412  1.00  7.79           C  
ATOM    679  O   ALA A  1       50.019  13.249  18.412  1.00  8.22           O  
ATOM    680  CB  ALA A  1       51.940  12.923  15.687  1.00  8.29           C  
ATOM    684  N   ALA A  2       49.571  14.604  16.684  1.00  8.07           N  
ATOM    685  CA  ALA A  2       48.167  14.680  17.048  1.00  9.23           C  
ATOM    61   C   ALA A  2       47.475  13.403  16.587  1.00  9.70           C  
ATOM    62   O   ALA A  2       47.172  13.219  15.405  1.00 11.16           O  
ATOM    63   CB  ALA A  2       47.530  15.925  16.434  1.00 10.43           C  
ATOM    693  N   ALA A  3       47.254  12.510  17.547  1.00  9.01           N  
ATOM    694  CA  ALA A  3       46.629  11.223  17.291  1.00  8.2            C  
ATOM    695  C   ALA A  3       45.109  11.265  17.261  1.00  9.04           C  
ATOM    696  O   ALA A  3       44.476  12.058  17.955  1.00 10.20           O  
ATOM    697  CB  ALA A  3       47.062  10.214  18.361  1.00  8.14           C  
END" > 1tri.pdb

# LIST OF SEQUENCE PEPTIDES
echo "#!/usr/bin/env python
#from array import array
import sys
a1 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q','R', 'S', 'T', 'V', 'W', 'Y']
#print a1[2]

# GENERATING ALL POSSIBLE COMBINATIONS OF AMINO ACIDS IN A PENTA PEPTIDE THAT HAS C IN POSITION 1 AND 5

for aa in a1:
 for ab in a1:
  for ac in a1:
   for ad in a1:
    print aa, ab, ac, ad
" > tripep.py

chmod 755 tripep.py

./tripep.py | sed -e 's/ //g' > targets.dat

# HERE START INSTRUCTION FOR THE MODELLER PROGRAM TO BUILD THE MODELS
# THE CODE WITH INSTRUCTIONS WILL BE PRINTED ON A GENERIC FILE "modeller.dat"
# BY DEFAULT A PENTAPEPTIDE WITH GENERIC SEQUENCE C...C, CYCLIZED BY A DISULFIDE BETWEEN RES 1 AND 5
# TO MODIFY MODELING PARAMETERS CHECK THE SALI'S WEBSTE https://salilab.org/modeller/manual/

echo "from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
import sys
# INSTRUCTIONS TO ADD DISULFIDE, N- AND C-TERMINAL PATCHES (DISU, NTER, CTER)
class mymodel(allhmodel):
    def special_patches(self, aln):
       #self.patch(residue_type='ACE', residues=self.residues['1'])
       #Standard C terminal patch
       #self.patch(residue_type='CT2', residues=self.residues['5'])
       #self.patch(residue_type='DISU', residues=(self.residues['1'],
       #                                          self.residues['5']))
       self.patch(residue_type='CTER', residues=self.residues['4']) 
       self.patch(residue_type='NTER', residues=self.residues['1'])
log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in
env.patch_default = False
# directories for input atom files
env.io.atom_files_directory = './:../atom_files'
a = mymodel(env, alnfile='XXXX.ali',
              knowns='1triA', sequence='XXXX',
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
#GENERATES 20 MODELS FOR EACH PEPTIDE
a.starting_model = 1
a.ending_model =  20
# Very thorough VTFM optimization:
# a.library_schedule = autosched.slow
# minimizzazione della struttura col conjugate gradient
a.md_level=refine.fast
a.library_schedule = autosched.slow
a.max_var_iterations = 20000
a.repeat_optimization = 2
a.make()

# Get a list of all successfully built models from a.outputs
ok_models = [x for x in a.outputs if x['failure'] is None]

# Rank the models by DOPE score
key = 'DOPE score'
if sys.version_info[:2] == (2,3):
    # Python 2.3's sort doesn't have a 'key' argument
    ok_models.sort(lambda a,b: cmp(a[key], b[key]))
else:
    ok_models.sort(key=lambda a: a[key])

# GET TOP MODEL
m1 = ok_models[0]
m2 = ok_models[1]
m3 = ok_models[2]
m4 = ok_models[3]
m5 = ok_models[4]
m6 = ok_models[5]
m7 = ok_models[6]
m8 = ok_models[7]
m9 = ok_models[8]
m10 = ok_models[9]
m11 = ok_models[10]
m12 = ok_models[11]
m13 = ok_models[12]
m14 = ok_models[13]
m15 = ok_models[14]
m16 = ok_models[15]
m17 = ok_models[16]
m18 = ok_models[17]
m19 = ok_models[18]
m20 = ok_models[19]
print(\"ranocchia %s ( %.3f)\" % (m1['name'], m1[key]))
print(\"ranocchia %s ( %.3f)\" % (m2['name'], m2[key]))
print(\"ranocchia %s ( %.3f)\" % (m3['name'], m3[key]))
print(\"ranocchia %s ( %.3f)\" % (m4['name'], m4[key]))
print(\"ranocchia %s ( %.3f)\" % (m5['name'], m5[key]))
print(\"ranocchia %s ( %.3f)\" % (m1['name'], m6[key]))
print(\"ranocchia %s ( %.3f)\" % (m2['name'], m7[key]))
print(\"ranocchia %s ( %.3f)\" % (m3['name'], m8[key]))
print(\"ranocchia %s ( %.3f)\" % (m4['name'], m9[key]))
print(\"ranocchia %s ( %.3f)\" % (m5['name'], m10[key]))
print(\"ranocchia %s ( %.3f)\" % (m1['name'], m11[key]))
print(\"ranocchia %s ( %.3f)\" % (m2['name'], m12[key]))
print(\"ranocchia %s ( %.3f)\" % (m3['name'], m13[key]))
print(\"ranocchia %s ( %.3f)\" % (m4['name'], m14[key]))
print(\"ranocchia %s ( %.3f)\" % (m5['name'], m15[key]))
print(\"ranocchia %s ( %.3f)\" % (m1['name'], m16[key]))
print(\"ranocchia %s ( %.3f)\" % (m2['name'], m17[key]))
print(\"ranocchia %s ( %.3f)\" % (m3['name'], m18[key]))
print(\"ranocchia %s ( %.3f)\" % (m4['name'], m19[key]))
print(\"ranocchia %s ( %.3f)\" % (m5['name'], m20[key]))
" > modeller.dat

# GENERIC ALIGNMENT
echo ">P1;1triA
structureX:1tri:1    :A:3  :A:templtrialanin:every organism:     :     
AAA*

>P1;XXXX
sequence:XXXX:1    : :3    : :targettripep:no organism:     :     
XXXX*
" > align.dat


for file in `cat targets.dat`; do sed -e 's/XXXX/'$file'/' align.dat > $file.ali ; done
for file in `cat targets.dat`; do sed -e 's/XXXX/'$file'/' modeller.dat > $file.py ; done

# 4 runs of modeller (depend on the CPUs you've in your ws)
# modeller executable is mod9.25, modify this if you have another version)
for file in *.py
do
mod9.25 $file &
forky 4

done
wait

rm -rf targets.dat modeller.dat #align.dat
rm -rf ????.D????????
rm -rf ????.V????????
rm -rf *ini *sch *rsr *ali
rm -rf *D.py tripep.py tripep.log D.log

# GENERATES PDB STRUCTURES OF ALL THE PEPTIDES

# mkdir delenda
mkdir goods
for file in *.py
  do
   base=`basename $file .py`
   buone=(`cat $base".log" | grep ranocchia | awk '{print $2}'| head -n 1`)
   # tutte=(`ls $base.B999900??.pdb`)

   mv ${buone[@]} goods
   # mv ${tutte[@]} delenda
done
rm -rf ????.B999900??.pdb
mv goods/* .
rm -rf goods
#rm -rf delenda
rm -rf *log


for py in ????.py
do
   num=1
   base=`basename $py .py`
     for file in `ls $base.B999900??.pdb`
     do
     mv $file $base"."pdb
     ((num++))
   done
done
rm -rf ????.py 1tri.pdb align.dat 
