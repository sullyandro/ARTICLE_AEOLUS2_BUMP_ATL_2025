
## Materials for the article 



# An Intermediate Complexity Approach to the Dynamics of Localized Extreme Heatwaves in the Mid-Latitude Atmosphere using Aeolus2.0


### Authors

Sullyandro O. Guimarães<sup>1,2</sup>, Masoud Rostami<sup>1,3</sup>, Stefan Petri<sup>1</sup>

1. Potsdam Institute for Climate Impact Research (PIK), Member of the Leibniz Association, Potsdam, Germany
2. University of Potsdam, Potsdam, Germany
3. Laboratoire de Météorologie Dynamique (LMD), Sorbonne University (SU), Ecole Normale Supérieure (ENS), Paris, France


### Goals

Heatwaves in mid-latitudes stem from localized heating that can persist and spread through the atmosphere, but the underlying processes remain unclear. Using the Aeolus 2.0 model, built on a Moist-convective Thermal Rotating Shallow Water (mcTRSW) framework, we examine how buoyancy anomalies evolve in dry and moist settings. Local heating rapidly drives surface convergence, upper-level divergence, and paired circulations distorted by Earth’s rotation, yielding asymmetries, spiral rainbands, and feedbacks. Moist convection amplifies responses through latent heat release, sustaining rainfall and instability. Inertia–gravity and Rossby waves redistribute energy and reshape anomalies. Aeolus 2.0 thus links convection, wave dynamics, and circulation, illuminating mechanisms sustaining extreme heatwaves.

The routines and main data for the article production are available in this repository for reproducibility and transparency.<br/>
 <br/>
 <br/>


### Organizational structure of codes and data
 <br/>
 
### > Directory  "Aeolus2.0_model"

Here you find the core of the Atmospheric model, with small adaptations from <br/>
the published version in Rostami et al. (2024; https://doi.org/10.5281/zenodo.13987667).
			

+ <b>Aeolus2.0 model guide (adapted for this study)</b>

  - Requirements

    Anaconda3 (e.g. Miniconda3 64-Bit with Python 3.10) <br/>
    Python    (until 3.10 is recommended; Dedalus v2 showed conflicts with newer Python versions) <br/>
    Dedalus   (dedalus==2.2006) <br/>
    Python libs: climlab, xarray, pooch, future, mpmath, netCDF4 <br/>

  - Installation on Linux
    
    <b>[1]</b> Considering the current directory is "Aeolus2.0_model", create a directory for requirements installation: <br/>
 
    $ mkdir requirements_installation/ <br/>
    $ cd requirements_installation/ <br/>
			   
    <b>[2]</b> Download the Miniconda3 64-Bit (x86) Installer (https://docs.conda.io/projects/miniconda/en/latest/) <br/>
    
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-py310_24.3.0-0-Linux-x86_64.sh <br/>
    $ bash Miniconda3-py310_24.3.0-0-Linux-x86_64.sh <br/>

    <b>[3]</b> Instructions for Dedalus installation: <br/>
			
    [Detailed installation is given in https://dedalus-project.readthedocs.io/en/latest/pages/installation.html] <br/>
			   
    Download the Dedalus v2 conda installation script from this link: <br/>
    https://raw.githubusercontent.com/DedalusProject/dedalus_conda/master/conda_install_dedalus2.sh <br/>
			   
    Or use: <br/>
    $ curl https://raw.githubusercontent.com/DedalusProject/dedalus_conda/master/conda_install_dedalus2.sh --output conda_install_dedalus2.sh <br/>
    
    In this shell script conda_install_dedalus2.sh: <br/>

    Modify the options at the top of the script to change the name of the resulting conda environment,  <br/>
    link against custom MPI/FFTW/HDF5 libraries, choose between OpenBLAS and MKL-based numpy/scipy, and more. <br/>
    Change where PYTHON_VERSION="3.10" to match the Python version chosen. <br/>

    Activate the base conda environment and run the script to build a new conda environment with Dedalus and its dependencies, as requested: <br/>

    $ conda activate base <br/>
    $ bash conda_install_dedalus2.sh <br/>
			   
    To use Dedalus, you need to activate the new environment. You can test the installation using the command-line interface: <br/>
		   
    $ conda activate dedalus2 <br/>
    $ python3 -m dedalus test <br/>
			   
    <b>[4]</b> Additional required packages: <br/>

    $ pip install xarray pooch future mpmath netCDF4 <br/>
			
    Climlab (https://climlab.readthedocs.io/en/latest/installation.html) <br/>
    $ conda install -c conda-forge climlab <br/>

    If that yields strange messages about version conflicts, compile from source: <br/>

    $ git clone https://github.com/climlab/climlab.git <br/>
    $ cd climlab <br/>
    $ pip install . --no-deps -vv <br/>
 
    At this point, the requirements are completed. <br/>

  
    <b>[5]</b> Go back to the directory "Aeolus2.0_model": <br/>
    
    $ cd ../Aeolus2.0_model <br/>
    $ ls * <br/>
    ```
    Output:
    
    aeolus2main.py            [Main model code]
    config.py                 [Configuration script to set up the model]
    equations_aeolus2main.py  [Routine with equations used by the model]
    jacobi128.py              [Routine from Dedalus package adapted for Aeolus2.0]
    sphere128.py              [Routine from Dedalus package adapted for Aeolus2.0]
    sphere_wrapper.py         [Routine from Dedalus package adapted for Aeolus2.0]
    timesteppers.py           [Routine from Dedalus package adapted for Aeolus2.0]
				   
    NetCDFOutput.py           [Function for netcdf creation]
    npytomat.py               [Function for npy to mat conversion]
    npytonc.py                [Function for npy to netcdf conversion]
				   
    Topo_non_unigrid.mat      [Topographic data that the model can read and is adapted to the grid]
    albedo_smooth.mat         [Climatological monthly averaged albedo adapted to a smooth grid]
    
    README.txt

    aeolus2_postprocessing.py [Post-process the output (npz) into a single netcdf4]   
    ```

    The files listed above can be found without modifications in the Aeolus2.0 package (https://doi.org/10.5281/zenodo.13987667). <br/>


  - To run Aeolus2.0 (in directory "Aeolus2.0_model"): <br/>

    $ conda activate dedalus2 <br/>
			
    Using 1 core: <br/>
    $ python3 aeolus2main.py <br/>
			
    Using 4 cores: <br/>
    $ mpiexec -n 4 python3 aeolus2main.py <br/>
			
    Important note:  <br/>
    Use multiples of 4 for the number of cores. The internal grid discretization requires that configuration for the solver. <br/>
		
		
  - Post-processing: <br/>
		
    $ python aeolus2_postprocessing.py <br/>
    <br/>


+ <b>Experiment Cases of this study:</b>

    The different cases of buoyancy anomaly were generated by varying the model configuration as follows.
		
    Moist-convective configuration will be frequently referred to as MC.
				
    In <b>config.py</b>:
    ```
    output_folder    = '../Data/Aeolus2.0_Output_[...]'
		
    moist_convection = 0 or 1 (w.r.t. dry case and moist-convective case)
		
    external_forcing = 0 or 'barotropic' or 'baroclinic' (w.r.t. no forcing; or in all layers; or in lower layer only)
		
    external_forcing_epsilon = 0.1 or 0.2 (w.r.t. maximum amplitude for the artificial buoyancy anomaly forcing)
    ```
		
    In total, we performed 10 cases simulations:
    ```
    Dry Control: control without external forcing, in dry situation (moist_convection = 0, external_forcing = 0)
    MC  Control: control without external forcing, in MC  situation (moist_convection = 1, external_forcing = 0)

    Dry Barotropic Weak:   weak   barotropic external forcing, in dry situation
                           (moist_convection = 0, external_forcing = 'barotropic', external_forcing_epsilon = 0.1)
    
    Dry Barotropic Strong: strong barotropic external forcing, in dry situation
                           (moist_convection = 0, external_forcing = 'barotropic', external_forcing_epsilon = 0.2)
		
    Dry Baroclinic Weak:   weak   baroclinic external forcing, in dry situation
                           (moist_convection = 0, external_forcing = 'baroclinic', external_forcing_epsilon = 0.1)
    
    Dry Baroclinic Strong: strong baroclinic external forcing, in dry situation
                           (moist_convection = 0, external_forcing = 'baroclinic', external_forcing_epsilon = 0.2)

    MC  Barotropic Weak:   weak   barotropic external forcing, in MC  situation
                           (moist_convection = 1, external_forcing = 'barotropic', external_forcing_epsilon = 0.1)
    
    MC  Barotropic Strong: strong barotropic external forcing, in MC  situation
                           (moist_convection = 1, external_forcing = 'barotropic', external_forcing_epsilon = 0.2)
		 
    MC  Baroclinic Weak:   weak   baroclinic external forcing, in MC  situation
                           (moist_convection = 1, external_forcing = 'baroclinic', external_forcing_epsilon = 0.1)
    
    MC  Baroclinic Strong: strong baroclinic external forcing, in MC  situation
                           (moist_convection = 1, external_forcing = 'baroclinic', external_forcing_epsilon = 0.2)
    ```

    To reproduce them, modify the <b>config.py</b> with the above-mentioned parameters for each case, and run the <b>aeolus2main.py</b> accordingly.
  
 <br/>
 
This repository has a Shell script called <b>Run_All_Steps.sh</b>, summarizing all steps to be executed to reproduce the data and figures. 

 <br/>
 
### > Directory "Data"

Data prepared or used for the study. NetCDF4 format is the most used. <br/>
	
#### Subdirectories:

+ <b>"Aeolus2.0_Input_ERA5"</b>:  ERA5 temperature and wind preprocessed to be used as start condition for Aeolus2.0.
	
+ <b>"Aeolus2.0_Output_Dry_Control"</b>: Aeolus2.0 output from case "Dry Control".
+ <b>"Aeolus2.0_Output_Dry_Baroclinic_Strong"</b>: Aeolus2.0 output from case "Dry Baroclinic Strong".
+ <b>"Aeolus2.0_Output_Dry_Baroclinic_Weak"</b>:  Aeolus2.0 output from case "Dry Baroclinic Weak".
+ <b>"Aeolus2.0_Output_Dry_Barotropic_Strong"</b>: Aeolus2.0 output from case "Dry Barotropic Strong".
+ <b>"Aeolus2.0_Output_Dry_Barotropic_Weak"</b>: Aeolus2.0 output from case "Dry Barotropic Weak".
 
+ <b>"Aeolus2.0_Output_MC_Control"</b>: Aeolus2.0 output from case "MC Control".
+ <b>"Aeolus2.0_Output_MC_Baroclinic_Strong"</b>: Aeolus2.0 output from case "MC Baroclinic Strong".
+ <b>"Aeolus2.0_Output_MC_Baroclinic_Weak"</b>: Aeolus2.0 output from case "MC Baroclinic Weak".
+ <b>"Aeolus2.0_Output_MC_Barotropic_Strong"</b>: Aeolus2.0 output from case "MC Barotropic Strong".
+ <b>"Aeolus2.0_Output_MC_Barotropic_Weak"</b>: Aeolus2.0 output from case "MC Barotropic Weak".
 <br/>
 
> [!NOTE]<br/>
> Given the amount of data, the simulations' output placed here were regridded from ~0.5 degrees to ~1 degree using <br/>
CDO (Climate Data Operator) "remapbil,r360x180". <br/>
> <br/>
> The data was checked after the interpolation, where, besides small differences, it has no significant deviation from the original. <br/>
> <br/>
> The original data can be produced by running the model as explained above, or obtained by contacting the authors. <br/>
	
 <br/>


### > Directory "Analysis_Scripts"

Production of the data and figures placed in this repository.

+ <b>Aeolus2_Data_Diff_Remap.py</b>: compute difference between Bump and Control cases, and regrid to ~1 degree (360x180).
  
+ <b>Aeolus2_Plot_Mapas.py</b>: prepare Divergence, Hamiltonian, Wind, and produce the maps. 

+ <b>Aeolus2_Plot_Series.py</b>: prepare Divergence, Hamiltonian, Wind, and produce the timeseries plots.
  
+	<b>Aeolus2_Plot_Hovmoller.py</b>: prepare Divergence, Hamiltonian, Wind, and produce the Hovmmoller plots.
	
The scripts listed here were made and executed in Python version 3.10.
	
Python libraries used are listed at the beginning of each Python script.

Execute Aeolus2_Difference_and_Remap.py first before the Plots.
	
 <br/>
 
### > Directory "Figures"

Generated figures for the article.

+ <b>"Plot_Maps"</b>

+ <b>"Plot_Hovmoller"</b>

+ <b>"Plots_Series"</b>
 
 <br/>
 
### Additional information

For additional material or information, contact the authors. <br/>
In this project, the authors used data from ERA5. <br/>
Citation of each data source is referred to in the article. <br/>
The authors are grateful to the institutes and projects that make the data and software available. <br/>
	
 <br/>

### Main references

Rostami, M. (2024). Open-Source Stand-Alone Version of Atmospheric Model Aeolus 2.0 Software. Zenodo. https://doi.org/10.5281/zenodo.13987667 <br/>

Rostami, M., Petri, S., Guimarães, S. O., Fallah, B. Open-source stand-alone version of atmospheric model Aeolus 2.0 Software. Geoscience Data Journal, 11, 1086–1093. DOI: 10.1002/gdj3.249 (2024). <br/>

Rostami, M., Severino, L., Petri, S., Hariri, S. Dynamics of localized extreme heatwaves in the mid-latitude atmosphere: A conceptual examination. Atmospheric Science Letters, 25(1), e1188. DOI: 10.1002/asl.1188 (2023). <br/>

Cao, Y., Kurganov, A., Liu, Y., Rostami, M., Zeitlin, V. On the dynamics of equatorial excited dipolar systems. Physics of Fluids, 37 (5), 056618. DOI: 10.1063/5.0270628 (2025). <br/>















