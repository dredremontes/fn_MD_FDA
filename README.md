# fn-MD-FDA-demo
Fibronectin (fn) steered molecular dynamics (MD) and force distribution analysis (FDA) tutorial
==========
This is an instruction manual for running steered molecular dynamics simulations of a single fibronectin subunit. All input files required to run the simulations are in `./fn-MD-FDA-demo/input/`. Completed simulation trajectories for postprocessing can be found in: `./fn-MD-FDA-demo/output/`. This tutorial walks you through how to get these outputs from scratch. However, you may choose to download the trajectory files and view the results in VMD. 

Downloading Gromacs
===========
This implementation was run in Ubuntu 20.04. The molecular dynamics software used was Gromacs 2020.4, which can be downloaded from https://ftp.gromacs.org/pub/gromacs-2020.4.tar.gz and is also available in the `./fn-MD-FDA-demo/.` The source code can be downloaded from https://manual.gromacs/org/2020.4/download.html. The tarball can be unpacked using:
```
$> tar -xvzf gromacs-2020.4.tar.gz
```
Please refer to the gromacs installation guide found in https://manual.gromacs.org/2020.4/install_guide/index.html. Check for Clang by running:
```
$> clang –version 
```
to check if Clang is already installed. If not, it can be installed by following the instructions here: https://www.ics.uci.edu/~pattis/common/handouts/macclion/clang.html. Next, download cmake from https://cmake.org/download/. 

First, we will move the gromacs-2020.4 folder from its current location to the folder where our programs are kept. For my PC, this folder is `C:/Program Files`. For MacOS, this will be the Applications folder. Next, we’re going to open up the terminal and navigate to the folder containing gromacs-2020.4. Now we can follow the "Quick and dirty installation" steps from the GROMACS website:
```
$> cd gromacs-2020.4
$> mkdir build
$> cd build
$> cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
$> make
$> make check
$> sudo make install
$> source /usr/local/gromacs/bin/GMXRC
```
Setting up the simulation
==========

Reminder: all input files required to run the simulations are in `./fn-MD-FDA-demo/input/`

If starting from scratch, download the `1fna.pdb` from https://www.rcsb.org/structure/1FNA. The amino acids 6-ARG and 7-ASP are incomplete, which creates problems later on. You can replace the 6-ARG and 7-ASP with the full ARG and ASP amino acids, respectively by using MODELLER's mutate.py function (https://salilab.org/modeller/wiki/Mutate_model). The file `./fn-MD-FDA-demo/input/fn.pdb` has this change done already. Convert the pdb file to a gro file:
```
$> gmx pdb2gmx -f fn.pdb -o fn.gro
```
Select AMBER99SB-ildn force field and TIP3P. Rotate the molecule, add a box, solvate, and add ions with:
```
$> gmx_fda editconf -f fn.gro -rotate 0 -25 10 -o fn.gro
$> gmx_fda editconf -f fn.gro -o fn.gro -box 12.5 5.0 5.0 -center 4 2.25 2.75
$> gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
$> gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
$> gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```
Select the SOL option.

Minimization
==========

NOTE: Minimization can take ~10-20 minutes on a modern laptop depending on how many cores you have. 
```
$> gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
$> gmx mdrun -v -deffnm em
```
Energy can be extracted to an xvg file using:
```
$> gmx energy -f em.edr -s em.tpr -o em.xvg
```
Equilibration
==========

NVT:
```
$> gmx grompp -f nvt.mdp -c em.gro -p topol.top -r em.gro -o nvt.tpr
$> gmx mdrun -deffnm nvt
```
NPT:
```
$> gmx grompp -f npt.mdp -c nvt.gro -p topol.top -r nvt.gro -t nvt.cpt -o npt.tpr
$> gmx mdrun -deffnm npt
```
The RMSD can be calculated with:
```
$> gmx rms -s npt.tpr -f npt.trr -o rmsd.xvg
```
The energy, temperature, and pressure can be extracted with:
```
$> gmx energy -f npt.edr -s npt.tpr -o output.xvg
```
Steered MD
==========

Define pull groups on fibronectin:
```
$> gmx make_ndx -f npt.gro
$> a1-26
$> name 20 pull
$> a1348-1367
$> name 21 hold
```
Create constraints
```
$> gmx genrestr -f npt.gro -n index.ndx -o R_A.itp
$> gmx genrestr -f npt.gro -n index.ndx -o R_B.itp
```
These position restraints need to be defined in their respective protein chain .itp files. For the restraint on 6-ARG, we can add the following to the bottom of the topol_Protein_chain_A.itp file:
```
#ifdef POSRES_R_A
#include "R_A.itp"
#endif
```
For the restraint on 96-ILE, we can add the following to the bottom of the topol_Protein_chain_B.itp file:
```
#ifdef POSRES_R_B
#include "R_B.itp"
#endif
```
We must open up the R_A.itp and R_B.itp files and manually renumber the left-most column, i. This column will contain the original atom numbers, which confuses GROMACS because position restraint (.itp) files always start with the number 1. So we have to go in and manually change ALL the atom numbers in R_A.itp and R_B.itp to 1, 2, 3 and so on. 

Add the following to the pull.mdp at the top if it is not already there:
```
define = -DPOSRES_R_A -DPOSRES_R_B
```
Running the steered MD.
```
$> gmx grompp -f pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr
$> gmx mdrun -deffnm pull -s pull.tpr -pf pullf.xvg -px pullx.xvg
```
Post-processing
===============

After running the simulations, we can extract the COM coordinates of the pull groups to calculate extension:
```
$> gmx trajectory -f sim_name.xtc -s sim_name.tpr -n index.ndx -ox pull-coords.xvg -seltype res_com -y yes
```
Select the two pull groups defined earlier to complete the selection and extraction. Gromacs will output a sim-namef.xvg file with the forces. Python can parse through the data files and after some manipulation, we can plot the force vs extension. Gromacs can also condense the trajectory files for faster viewing in VMD by extracting only the protein:
```
$> gmx trjconv -f sim-name.xtc -s sim-name.tpr -o new-traj.xtc
```
Then selecting option 1 to select the protein. We can then use VMD to load npt.gro and sim-name.xtc to create movies.

Force Distribution Analysis
==============
To run the force distribution analysis, we used the input file found in `./fn-MD-FDA-demo/input/input.pfi`, the steered MD trajectory of choice, and gromacs-fda from: https://github.com/HITS-MBM/gromacs-fda.
