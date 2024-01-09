#!/bin/bash

### INDEX GENERATION: Coupling groups + restraints groups
gmx make_ndx -f system_min.gro -o index.ndx<<EOF
r PA | r PC | r OL
name 8 Membrane
ri 2 & a P31
q
EOF

### FIRST EQUILIBRATION STEP
gmx genrestr -f system_GMX.gro -n index.ndx -o posre.itp -fc 0 0 1000 <<EOF
9
EOF
gmx grompp -f mdp/equi1.mdp -r system_min.gro -c system_min.gro -n index.ndx -p system_GMX.top -o system_equi1.tpr -maxwarn 2
gmx mdrun -deffnm system_equi1 -v

### SECOND EQUILIBRATION STEP
gmx genrestr -f system_GMX.gro -n index.ndx -o posre.itp -fc 0 0 500 <<EOF
9
EOF
gmx grompp -f mdp/equi2.mdp -r system_equi1.gro -c system_equi1.gro -n index.ndx -p system_GMX.top -o system_equi2.tpr -maxwarn 2
gmx mdrun -deffnm system_equi2 -v

### THIRD EQUILIBRATION STEP
gmx genrestr -f system_GMX.gro -n index.ndx -o posre.itp -fc 0 0 0 <<EOF
9
EOF
gmx grompp -f mdp/equi3.mdp -r system_equi2.gro -c system_equi2.gro -n index.ndx -p system_GMX.top -o system_equi3.tpr -maxwarn 2
gmx mdrun -deffnm system_equi3 -v
