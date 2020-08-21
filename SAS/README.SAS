
       README.SAS

Python scripts in this directory can be used to process and plot LES and SCM output. These scripts have been tested on NCAR's casper platform.

To use these scripts follow the steps below:

1. Log on cheyenne
2. To start an interactive session on caseper execute the following command:
   % execdav --account=NWSA0002 --time=12:00:00 --constraint=x11 --mem=150G
3. Once the interactive session has started execute commands:
   % module load ncarenv
   % ncar_pylib
4. Modify or add the path and filename in <NAME>analyis.py file to point to 
   the directory and file where LES output is stored.
   Here NAME can be one of the following: FE, NCAR, WRF, etc.
   Modify also path to and the name of the output file. 
5. To create a pickle file from LES or SCM output execute for example:
   % ./<NAME>analysis.py 
   Depending on the size of LES output this may take awhile.
6. Once pickle (*.pkl) file is created modify path to and file name of
   the pickle file in SASplotEvolution and SASplotProfile.  
   For output on the screen set 
   show_fig = True
   in python scripts.
7. Run plotting scripts by executing the following commands.
   % ./SASplotEvolution
   % ./SASplotEvolution
