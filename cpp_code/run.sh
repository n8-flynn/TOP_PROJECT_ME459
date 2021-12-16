#!/usr/bin/env zsh
#SBATCH --job-name=huz_final_project
#SBATCH --partition=wacc
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:30:00
#SBATCH -o run.out
#SBATCH --mem=8G

nelx=32
nely=20
volfrac=0.4
penal=3
rmin=1.2
./run $nelx $nely $volfrac $penal $rmin
echo "Calculations complete.To export result to PNG, running post prossesing script"
module load mamba
bootstrap_conda
conda activate py459
python3 post_pros.py $nelx $nely
echo "Done. View density_field.png for results"
