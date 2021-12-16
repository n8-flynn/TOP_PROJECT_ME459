#!/usr/bin/env zsh
nelx=32
nely=20
volfrac=0.4
penal=3
rmin=1.2
wh=2
./trial $nelx $nely $volfrac $penal $rmin $wh
echo "Running post prossesing script"
python3 post_pros.py $nelx $nely
echo "Done. View results in density_field.png or density_field.csv"
