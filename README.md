This README file contains instructions about how to compile the DYKE-CODE_GEOTWINS on a Linux OS (1), and run the simulation in the sub-folder “RUNS/TEST_01” (2). Here you will also find a description of the content of the repository “GEOMOD-T2.1_DYKE-CODE” (3).

IMPORTANT NOTE: In order to compile the code and run the executable file on a Linux operating system, gfortran compiler and LAPACK libraries must be installed.

1) Compile the FORTRAN90 source files for the DYKE-CODE_GEOTWINS on a Linux OS (tested on Ubuntu20.04).

1.1) Move to the folder “PROGRAM_FILES/DYKE-CODE_GEOTWINS”

1.2) From the terminal run the bash script “compile_DYKE-CODE_GEOTWINS.sh”:
> sh compile_DYKE-CODE_GEOTWINS.sh
This will generate the executable file “DYKE-SIMULATION_GEOTWINS” in the sub-folder
“MAIN”

2) Run the executable file “DYKE-SIMULATION_GEOTWINS”.
   
2.1) Copy the executable file “DYKE-SIMULATION_GEOTWINS” to the folder “RUNS/TEST_01”;

2.2) Delete the previous output files from the sub-folder "output";

2.3) Execute the “DYKE-SIMULATION_GEOTWINS” file. From the terminal, in the folder “RUNS/TEST_01” run the command:
> ./DYKE-SIMULATION_GEOTWINS
This will produce the output files in the folder “output”. Please refer to the
“outputfile_description.txt” for a description of the outputs.

3) Description of the content of the directory “GEOTWINS-T2.1_DYKE-CODE”.

- inputfile_description.txt and outputfile_description.txt: These files contain the description of the input-files needed to run the DYKE-SIMULATION_GEOTWINS executable file, and the description of the output-files generated the DYKE-SIMULATION_GEOTWINS.

- RUNS: contains one example of a DYKE-SIMULATION_GEOTWINS run (TEST_01).
  
  -- TEST_01: contains input and output files (in the “input” and “output” subfolders respectively) relative to the DYKE-SIMULATION_GEOTWINS execution, the DYKE-SIMULATION_GEOTWINS executable file, the c-basch script “plot_crack_shape.csh” that produces the figure “crack_shape.ps” provided that GMT4 software is installed.

    --- input: please refer to inputfile_description.txt.
  
    --- output: please refer to outputfile_description.txt.

- PROGRAM_FILES: contains all files needed to compile and run the DYKE-CODE_GEOTWINS
  
  -- DYKE-CODE_GEOTWINS: contains the FORTRAN90 program files for the DYKE-CODE_GEOTWINS, and the script to compile the code.
  
    --- MAIN: contains the main file of the DYKE-CODE_GEOTWINS (MAIN_PROPAGATION_GEOTWINS.f90), and the executable file (DYKE-SIMULATION_GEOTWINS).
  
    --- MODULES_F90: contains the module files of the DYKE-CODE_GEOTWINS (DISL2D.f90, CRACK2D.f90, EXTERNAL_FIELD.f90, BELEMENT_GEOTWINS.f90, COMP_FIELD_GEOTWINS.f90, PROPAGATION_GEOTWINS.f90).
  
    --- compile_DYKE-CODE_GEOTWINS.sh: bash script to compile the DYKE-CODE_GEOTWINS and produce the executable file “DYKE-SIMULATION_GEOTWINS”.
