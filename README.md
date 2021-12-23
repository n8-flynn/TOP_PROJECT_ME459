# Topology optimization of a 2D Cantilivered Beam

## How do you run our code?

CLone the entire repo into your local system. The first folder cpp_code contains all of the necessary files to build, run, and export a final output as a .png file.


Once the repository is appropriately cloned,there are two ways the project can be built and run; 1) Using Cmake, 2) Using the shell scripts provided. First cd into cpp\_code.Inputs can be verified in the run.sh file. Once inputs are set, follow the steps given below.

### Instructions using CMake
```bash
mkdir build
cd build
cmake ..
make
cp ../{run.sh,post_pros.py} .
sbatch run.sh
```

### Instructions using Shell scripts
```bash
sbatch compile.sh
sbatch run.sh
```

   
