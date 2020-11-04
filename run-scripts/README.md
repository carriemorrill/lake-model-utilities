# Running multiple iterations of Lake Model using shell scripting
This procedure is useful for calibrating the lake model across a range of possible parameter values, or for generating an ensemble of lake model simulations to sample the impact of parameter uncertainty.

1. In lake.inc, comment out the statements for parameters that you wish to vary. For example, if you want to test different values for eta, the parameter statement should look like:

   > !      parameter (eta = 0.2) 

2. Generate scaling factors for parameters you wish to vary using Latin Hypercube sampling in R. To do this in R, use the following commands. Set nsim and nparam to be the number of lake model simulations you plan to run and the number of parameters that you plan to vary.

   > install.packages("lhs")
   
   > library(lhs)
   
   > nsim = 1000; nparam = 7
   
   > improvedLHS(nsim,nparam)
   
   Save output as lake-params.txt
   
   Documentation for the function improvedLHS is at: https://www.rdocumentation.org/packages/lhs/versions/1.0.1/topics/improvedLHS

   Briefly, the arguments for this function are: number of sets of parameter scalings to generate, number of parameters. The command above will generate 1000 sets of 7 parameter scalings. 
   
   Note that this step is not needed if you already know what parameter values you want to use for an ensemble (i.e., lake model calibration is already complete).

3. Edit the fortran program writeParameter.f to match the parameters you wish to vary and their possible ranges of values. The Latin Hypercube sampling generates values in the range 0 to 1. This fortran program will be called by the shell script to apply these scalings to parameter ranges and then write out parameter values to be concatenated to the lake.inc file.

4. Edit shell script (Unix: runLakeModel.sh, Windows Powershell: runLakeModel.ps1) as necessary. The main edit needed will be to match the number of loop iterations to the number of sets of parameters (= first argument of improvedLHS).

5. Run script using the appropriate command:

    > ./runLakeModel.sh
    
    > ./runLakeModel.ps1
