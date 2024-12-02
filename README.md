# Identifying yield surfaces

## Important notes
Currently, the script only runs on Linux based operating systems. 

## Description
This project provides an interface for identifying yield functions. On the HPC, many serial simulations can be executed simultaneously. 

## Installation
The Python scripts are tested with Python interpreter **3.8.6** and **3.8.10**. The installation is seperated for local machine and HPC. The individual Python files require different packages to be installed. The file *libmulti.py* is not intendeted to be executed by the user. The execution of *id-fy.py* requires the packages
 * numpy
 * matplotlib
 * optuna
 * mpl_toolkits
 * multiprocessing
 * pymoo
 * compilation of the PythonFE

The execution of stl-nice.py requires:
 * numpy

Furthermore, the *liblapack-dev* package needs to be installed on the machine.

### Installation on local machine
Firstly check if Python **3.8.10** or **3.8.6** are installed. It is recommended to use virtual environments. After activating the virtual environments, the installation of the packages can be handled by pip.

### Installation on HPC
On the HPC, it is efficient to only run the simulations with the corresponding python scripts. The plotting of the data can be handled on the local machine better. For the usage on the HPC, virtual environments are highly recommended. Create a virtual environment in the directory, where the simulations will be started. Afterwards, the modules **release/23.04**, **GCCcore/10.2.0** and **Python/3.8.6** need to be loaded in this specific order (Names can vary depending on the specific HPC). Subsequently, the required packages can be installed via pip. Alternatively, the packages can be loaded as seperate modules. This increases runtime, since each processors has to load all packages. Loading all required packages is not tested.

### Installation of PythonFE
The provided finite element source code, which solves one element, needs to be compiled via f2py from numpy. The numpy versions 1.17.4, and 1.23.5 are tested for compilation. Within the directory, there is a shell script, which compiles the Fortran source code and generates the required .so and .pyf files. These need to be copied into the directory, where the execution of the script takes place. Since the solution routine `dposv` from lapack is used, the *liblapack-dev* package needs to be installed on the compiling machine.

## Usage
The usage depends whether the script is executed on a local machine or on the HPC. For both methods, the parameters at the top of id-fy.py need to be set accordingly. The job to run can either be set in the python file or be given on execution via the command line. The minimal files, which are always necessary for identification contains id-fy.py, libmulti.py, classforpymoo.py and the compiled executables from the PythonFE.

### Local machine
For execution on the local machine, the parameters after *#Python* block do not need to be set anymore. If all packages are available, the command `python3 id-fy.py sample-for-run` will create files containing the points that need to be tested on the hydrostatic axis for each processor and the search paths in the deviatoric plane. An alternative sampling, where the search paths are divided onto the processors, is also possible be setting `hydroforprocs = False` in id-fy.py. The command `python3 id-fy.py parallel-<optimiser>` will start the identification procedure on all processors at the same time. `<optimiser>` needs to be replaced with the identifier for the corresponding optimiser. *ES* is will use the genetic algorithm from the pymoo package. *optuna* uses the optimiser from the optuna package. *nr* uses a Newton Raphson search with numerically approximated gradients. The Newton-Raphson search is self implemented. After the identification endedn, `python3 id-fy.py joinfiles` will gather the results saved per processor and generate one final result file. 

For plotting the results as scatter plots, the file *pure-plot.py* can be used. Firstly, the file name of the result file of the identification procedure needs to be set in *pure-plot.py*. When all Python packages are loaded, issuing the command `python3 pure-plot.py` will produce a plot from the isometric view and a plot along the hydrostatic axis as png files. The specific settings can be adjusted in *pure-plot.py*. 

For generating stl files *stl-output.py* can be used. Firstly, the file name, containing all the results of the identification in one file, needs to be set in *stl-output.py*. When all Python packages are loaded, the command `python3 stl-nice.py` will generate the stl file.

### HPC
On HPC environments, running the algorithm differs. It is only sipported to run everything by using virtual environments. After activating the virtual environment, the command `python3 id-fy.py sample-for-run` needs to be executed to create input files for each processor. Afterwards, a running script for all new HPC clusters of TU Dresden can be created by `python3 id-fy.py create-input-barnard`. This will produce an .sh file that can be submitted via `sbatch <name of the script>` to SLURM. Then, the HPC will execute the simulations. After obtaining the result files for each processor, the command `python3 id-fy.py joinfiles` produces the final result file. The final result file should be copied to the local machine for plotting purposes.

### Mininmal and maximal yield surface
If anisotropic material models are present, a search for the minimal and maximal yield surface can be conducted. Currently, the anisotropy needs to be handled via rotation angles. Extensions are possible. Therefore, the flag doanisotropic needs to be set to true. The variable rotationspecs contains the identifiers, which are replaced in the input file with the anlges of rotation. If the Python FE is used, the positions in the materialarray for the rotation angles need to be set in pferotations. Anglesteps contains the amount of steps, which are generated between the minimum and maximum per angle.

### Interface to pymoo package
The interface to the pymoo package is handled via a class contained in classforpymoo.py

## Extension to other yield functions
If other macroscopic yield functions should be incorporated, only modifications of `plastpyieldfunctions.f` are required. In the subroutine in this file, there can be other cases for idfy. There only fy needs to be calculated. The material parameters are set in id-fy.py

## Extension to other FE software
The source code can be extended to work with other FE software. Currently, usage with the FE software FEAP is supported. For the incorporation of other software, the function dofeap2 needs to be adjusted. This function takes the point in the eigenstress space as input (vall) and returns the kind of behavior and the value of the yield function (resi). The API to the software needs to be adjusted. The input file needs to be modified so that the point on the hydrostatic axis is tested. Furthermore, the FE software needs to provide the values of the yield function. This API needs to be programmed by the user.

## Theory
A scientific paper about the theory of this algorithm is available at the [International Journal of Plasticity](https://doi.org/10.1016/j.ijplas.2024.104183).

## Future Features
A list of things that can be included in the future.
 * support for other simulation software (FEM or MPM)
 * inclusion of further optimizers
 * more efficient minimal and maximal search implementation
 * rotation of boundary conditions for min-max runs

## License
When using this framework in publications, it is required to cite the publication of the algorithm.

## Known issues
Starating with numpy versions of 1.24.0, there are incompatibilities with package pymoo.

## Project status
Project is under partially under development.
