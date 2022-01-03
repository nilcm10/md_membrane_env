#!/bin/bash

### INDEX GENERATION: Coupling groups + restraints groups
gmx make_ndx -f system_min.gro -o index.ndx<<EOF 
r PA | r PC | r OL | r CHL
name 9 Membrane
r PC | r CHL & a P31 | a O1
name 10 headgroups
q
EOF

### FIRST EQUILIBRATION STEP
gmx genrestr -f system_GMX.gro -n index.ndx -o posre.itp -fc 0 0 1000 <<EOF
headgroups
EOF
gmx grompp -f mdp/equi1.mdp -r system_min.gro -c system_min.gro -n index.ndx -p system_GMX.top -o system_equi1.tpr -maxwarn 1
gmx mdrun -deffnm system_equi1 -v

### SECOND EQUILIBRATION STEP
gmx genrestr -f system_GMX.gro -n index.ndx -o posre.itp -fc 0 0 500 <<EOF
headgroups
EOF
gmx grompp -f mdp/equi2.mdp -r system_equi1.gro -c system_equi1.gro -n index.ndx -p system_GMX.top -o system_equi2.tpr -maxwarn 1
gmx mdrun -deffnm system_equi2 -v

### THIRD EQUILIBRATION STEP
gmx genrestr -f system_GMX.gro -n index.ndx -o posre.itp -fc 0 0 0 <<EOF
headgroups
EOF
gmx grompp -f mdp/equi3.mdp -r system_equi2.gro -c system_equi2.gro -n index.ndx -p system_GMX.top -o system_equi3.tpr -maxwarn 1
gmx mdrun -deffnm system_equi3 -v