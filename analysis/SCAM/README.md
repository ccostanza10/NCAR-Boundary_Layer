## Scripts and code for the analysis of Single Column CAM (SCAM) output for the NCAR PBL project modeling cases
* Structured as python in a main jupyter notebook with separate python compute and plotting functions
* It is intended for the cmparison either of multiple SCAM cases, or comparison between SCAM and LES cases
* Currently there is 1D and 2D timeseries and vertical profile time snapshot plotting
* In the future will add time mean vertical profiles, profile budgets, animations and CLUBB PDF illustrations 

### Files (for basic workflow: Parameter selections, multiple SCAM runs, analysis scripts, output)
* scam_plots.ipynb : Main driving notebook
* scam_fplot.py : Plotting routines
* scam_var_defn.py : Variable definitions for plotting (metadata look ups) 
* scam_utils.py : Mulitple utilities (reading SCAM/LES data, derived variables etc.)

### Other routines
* scam_forcings_plev_flex.ipynb : Generate forcings for IOPs (SAS yes, PERDIGAO setting up, RICO soon)
* scam_forcings_iop.py : IOP specific settings
* scam_nc2sas.ipynb : Interpreter from NCAR-LES format to SCAM netcdf file format
* scam_nc2sas_vars.py : Variables mapping from NCAR-LES to SCAM names where possible.
* scam_plot_cases.ipynb : Stand aone simple SCAM case comparison tool, for the above and other standard IOP cases.
* scam_plot_pdfs.ipynb : (Under development) Scripts to reconstruct Gaussian PDFs using CLUBB output in SCAM

