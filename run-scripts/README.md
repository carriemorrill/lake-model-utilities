# Running multiple iterations of Lake Model using shell scripting
This procedure is useful for calibrating the lake model across a range of possible parameter values, or for generating an ensemble of lake model simulations to sample the impact of parameter uncertainty. Instructions below explain how to run a script in either a Unix shell or a Windows Powershell.

1. In lake.inc, comment out the statements for parameters that you wish to vary. Save the new lake.inc file as lake.inc.save. For example, if you want to test different values for eta, the parameter statement should look like:

   > !      parameter (eta = 0.2) 

2. Generate scaling factors for parameters you wish to vary using Latin Hypercube sampling in R. To do this in R, use the following commands. Set nsim and nparam to be the number of lake model simulations you plan to run and the number of parameters that you plan to vary.

   > install.packages("lhs")
   
   > library(lhs)
   
   > nsim = 1000; nparam = 7
   
   > improvedLHS(nsim,nparam)
   
   Save output as lake-params.txt. Output will have nsim rows and nparam columns. The command above will generate 1000 sets of 7 parameter scalings. 
   
   Documentation for the function improvedLHS is at: https://www.rdocumentation.org/packages/lhs/versions/1.0.1/topics/improvedLHS
   
   Note that this step is not needed if you already know what parameter values you want to use for an ensemble (i.e., lake model calibration is already complete). Instead, create a file called lake-params.txt and save actual parameter values there. Each row of lake-params.txt should be the set of parameter values to be used in one lake model simulation.

3. Edit the fortran program writeParameter.f to reflect the parameters you wish to vary and their possible ranges of values. The Latin Hypercube sampling generates values in the range 0 to 1. This fortran program will be called by the shell script to apply these scalings to parameter ranges and then write out parameter values to be concatenated to the lake.inc file.

   Note that if lake-params.txt contains actual parameter values, and not Latin Hypercube values, then the lines in writeParameter.f that apply scalings should be commented out.

4. Edit shell script (Unix: runLakeModel.sh, Windows Powershell: runLakeModel.ps1). First, match the number of loop iterations to the number of sets of parameters (= number of rows in lake-params.txt). For example, when lake-params.txt contains 1000 rows, then for Unix and Powershell, respectively:

    > while [$j -le 1000]

    > while($j -le 1000)
   
    Second, edit the fortran compiler and the name of the lake model version you are running, for example: 
 
    > gfortran -o lake HB-windy-sed-varD2.f90
    
    Third, the script assumes that the lake model generates one output file called surface.dat. The script appends the simulation number to the name of this file, so that it will not be overwritten by the next lake model simulation. If the output file has a different name, the line below will need editing to reflect that. Also, if there are multiple output files, additional mv commands are necessary to rename each additional output file.
    
    > mv surface.dat surface-"$j".txt

5. Run script using the appropriate command:

    > ./runLakeModel.sh
    
    > ./runLakeModel.ps1
   
