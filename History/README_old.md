# fn_MD_FDA
Fibronectin (fn) steered molecular dynamics (MD) and force distribution analysis (FDA) tutorial
==========
This is an instruction manual for running steered a molecular dynamics simulation (constant pulling force of 700 pN) of a single fibronectin subunit. All input files required to run the simulations are in `./fn_MD_FDA/input/`. Completed simulation trajectories for postprocessing can be found in: `./fn_MD_FDA/output/`. This tutorial walks you through how to get these outputs from scratch and all the simulations can be done on a modern laptop with 1-3 hour runtime depending on specifications. However, you may choose to download the trajectory files and view the results in VMD. 

Downloading Gromacs
===========
On Linux, the most straightforward way to download gromacs is by entering `sudo apt install gromacs` in the command line. However, another option is to download the source code and build it manually, which is described below for Gromacs 2020.4.

This implementation was run in Ubuntu 20.04. The molecular dynamics software used was Gromacs 2020.4, which can be downloaded from https://ftp.gromacs.org/pub/gromacs/ and is also available in the `./fn_MD_FDA/.` The source code can be downloaded from https://manual.gromacs.org/2020.4/download.html. The tarball can be unpacked using:
```
tar -xvzf gromacs-2020.4.tar.gz
```
Please refer to the gromacs installation guide found in https://manual.gromacs.org/2020.4/install-guide/index.html. 

If you run into installation errors, you may need to install compilers and cmake. To get the latest version of the C and C++ compilers, run:
```
sudo apt-get install g++
```
Then:
```
sudo apt install build-essential
```
To install the latest version of cmake:
```
sudo apt get install cmake
```
First, we will move the gromacs-2020.4 folder from its current location to the folder where our programs are kept. For my PC, this folder is `C:/Program Files`. For MacOS, this will be the Applications folder. Next, we’re going to open up the terminal and navigate to the folder containing gromacs-2020.4. Now we can follow the "Quick and dirty installation" steps from the GROMACS website:
```
cd gromacs-2020.4
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make
make check
sudo make install
source /usr/local/gromacs/bin/GMXRC
```
Setting up the simulation
==========

Reminder: all input files required to run the simulations are in `./fn_MD_FDA/input/`

If starting from scratch, download `1fna.pdb` from https://www.rcsb.org/structure/1FNA. The amino acids 6-ARG and 7-ASP are incomplete, which creates problems later on. You can replace the 6-ARG and 7-ASP with the full ARG and ASP amino acids, respectively by using MODELLER's mutate.py function (https://salilab.org/modeller/wiki/Mutate_model). The file `./fn_MD_FDA/input/fn.pdb` has this change done already. Convert the pdb file to a gro file:
```
gmx pdb2gmx -f fn.pdb -o fn.gro
```
Select AMBER99SB-ildn force field and TIP3P. 
Rotate the molecule
```
gmx editconf -f fn.gro -rotate 0 -25 10 -o fn.gro
```
Add a box
```
gmx editconf -f fn.gro -o box.gro -box 12.5 5.0 5.0 -center 4 2.25 2.75
```
Solvate
```
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
```
Add ions:
```
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
```
```
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```
Select the SOL option.

Minimization
==========

NOTE: Minimization can take ~10-20 minutes on a modern laptop depending on how many cores you have. 
```
gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
```
```
gmx mdrun -v -deffnm em
```
Energy can be extracted to an xvg file using:
```
gmx energy -f em.edr -s em.tpr -o em.xvg
```
Equilibration
==========

NVT:
```
gmx grompp -f nvt.mdp -c em.gro -p topol.top -r em.gro -o nvt.tpr
```
```
gmx mdrun -v -deffnm nvt
```
NPT:
```
gmx grompp -f npt.mdp -c nvt.gro -p topol.top -r nvt.gro -t nvt.cpt -o npt.tpr
```
```
gmx mdrun -v -deffnm npt
```
The RMSD can be calculated with:
```
gmx rms -s npt.tpr -f npt.trr -o rmsd.xvg
```
And then selecting the protein twice. The energy, temperature, and pressure can be extracted with:
```
gmx energy -f npt.edr -s npt.tpr -o output.xvg
```
Steered MD
==========

Define pull groups on fibronectin:
```
gmx make_ndx -f npt.gro
```
Type:
```
a1-26
```
Press enter/return. Check to see the next available number slot. For example, if the group number stops at 18, use the next available number, which would be 19:
```
name 19 pull
```
Press enter/return. Type:
```
a1348-1367
```
Press enter/return. This selected the residue 6-ARG (Arginine) as the amino acid group to pull, we now have to select the other amino acid group to hold. Type:
```
name 20 hold
```
Press enter/return twice. This selected the residue 96-ILE (Isoleucine) as the hold group. Press q and then enter to save and quit. 
Next, we create constraints
```
gmx genrestr -f npt.gro -n index.ndx -o R_pull.itp
```
Type the number next to the pull group and press enter.
```
gmx genrestr -f npt.gro -n index.ndx -o R_hold.itp
```
Type the number next to the hold group and press enter.

These position restraints also need to be defined in their respective protein chain .itp files. For the restraint on 6-ARG (residue 6, Arginine), we can add the following to the bottom of the topol.top and posre.itp files:
```
#ifdef R_hold
#include "R_hold.itp"
#endif
```
For the restraint on the hold group, 96-ILE (residue 96, Isoleucine), we can add the following to the bottom of the topol.top and posre.itp file:
```
#ifdef R_pull
#include "R_pull.itp"
#endif
```
We must open up the R_hold.itp file and manually renumber the left-most column, i. This column will contain the original atom numbers, which confuses GROMACS because position restraint (.itp) files always start with the number 1. So we have to go in and manually change ALL the atom numbers in R_hold.itp to 1, 2, 3 and so on. 

Add the following to the pull.mdp at the top if it is not already there:
```
define = -DR_hold -DR_pull
```
Running the steered MD.
```
gmx grompp -f pull_fn.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull_fn.tpr
```
```
gmx mdrun -v -deffnm pull_fn -s pull_fn.tpr -pf pullforce.xvg -px pullx.xvg
```

Force Distribution Analysis
==============
To run the force distribution analysis, gromacs-fda will need to be installed from https://github.com/HITS-MBM/gromacs-fda. Use the input file found in `./fn_MD_FDA/input/input.pfi` and the steered MD trajectory of choice: 
```
gmx_fda mdrun -nt 1 -rerun pull_fn.xtc -pfi input.pfi -pfn index.ndx -s pull_fn.tpr -psr output.psr
```
Post-processing - Trajectory
===============

Download and open VMD 1.9.4a https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD. Click `file > New molecule > browse` and open `./fn_MD_FDA/output/npt.gro`. In main VMD window, right click the `npt.gro` molecule (which should be highlighted) and click `Load data into molecule` and select `./fn_MD_FDA/output/pull_fn.xtc`  Then click `Load` in the Molecule File Browser window. The trajectory should now be loaded in the display. 

Go to `Graphics > Representations` and under `Drawing Method` select `NewCartoon`. The cartoon representation should be displayed now and the water and ions should be hidden. 

Now go to `Extensions > VMD Preferences` and go to the Custom tab. Click `New` and name it "pf_loaduser". In the "code" box, write the following lines, but replacing the "path/to/pf_loaduser.tcl" with the location of the `pf_loaduser.tcl` script. This script is available in the `fn_MD_FDA` directory. In my local PC, I moved the `pf_loaduser.tcl` file to `"C:/Program Files/VMD/plugins/WIN64/tcl/pf_loaduser.tcl"`.
```
source "path/to/pf_loaduser.tcl"
package require pf_loaduser
```
Click "Update" and then "Push All Settings to VMD". Close out of the VMD Preferences Panel. 

Now go to `Extensions > TkConsole` and type `pbc box` and hit `enter` to display the simulation cell. Now, in the TkConsole, using cd and ls, navigate to the folder containing the cloned repository and find the output.psr file in `path/to/fn_MD_FDA/output/`. For my local PC, I use `cd "C:/Users/andre/Documents/GitHub/fn_MD_FDA/output"`. Remove the first frame by typing the command:
```
animate delete beg 0 end 0 skip 0 0
```
Source the pf_loaduser.tcl file in the TKConsole, replacing `path/to/pf_loaduser.tcl` with the directory that contains your `pf_loaduser.tcl` file:
```
source "path/to/pf_loaduser.tcl"
```
Load the pf_loaduser script:
```
package require pf_loaduser
```
The following command will overlay the FDA results onto the molecule:
```
pf_loaduser "output.psr" true 1 BWR
```
The TkConsole should print out "now loading per-residue data" and after a few seconds to a minute, the results should be shown on the molecule itself. You can scroll through the trajectory and see how the different regions within the fibronectin beta sheets activate to resist the unwinding. The final video can be seen in `fn_FDA_traj.mp4`.