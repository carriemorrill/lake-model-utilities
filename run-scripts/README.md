# Running multiple iterations of Lake Model using shell scripting
This procedure is useful for calibrating the lake model across a range of possible parameter values, or for generating an ensemble of lake model simulations to sample the impact of parameter uncertainty.
1. Generate scaling factors for parameters using Latin Hypercube sampling in R. To do this in R, use the following commands:

   install.packages("lhs")
   
   library(lhs)
   
   improvedLHS(1000,7)
   
Documentation for the function improvedLHS is at: https://www.rdocumentation.org/packages/lhs/versions/1.0.1/topics/improvedLHS

Briefly, the arguments for this function are: number of sets of parameter scalings to generate, number of parameters. The command above will generate 1000 sets of 7 parameter scalings. Note that this step is not needed in the case of generating an ensemble of lake model simulations with pre-determined parameter values (i.e., lake model calibration is already complete).

2. Specify ranges of parameter values in the fortran program writeParameter.f. The Latin Hypercube sampling generates values in the range 0 to 1. This fortran program will be called by the shell script to apply these scalings to parameter ranges and then write out parameter values to be concatenated to the lake.inc file.

3. Edit shell script runLakeModel.sh as necessary. The main edits needed will be to match the number of loop iterations to the number of sets of parameters (= first argument of improvedLHS = number of rows in the Latin Hypercube sampling) and to specify the directory to which lake model output files should be moved.

4. Run bash script using the following command:

    ./runLakeModel.sh
