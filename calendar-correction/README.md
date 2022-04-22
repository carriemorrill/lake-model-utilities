# Calendar correction for climate model output
These steps document how to run the Bartlein and Shafer (https://doi.org/10.5194/gmd-12-3889-2019) calendar correction code in python.  

You will need:  
- <b>Fortran compiler:</b> the example below uses mingw32 
- <b>Python numpy:</b> necessary for using f2py  
- <b>bartlein_sub.f90:</b> the Bartlein and Shafer (2019) code edited for use in python
- <b>(optional) spatialCC.py:</b> example python code showing how to read a PMIP netcdf, calendar correct, and write a corrected netcdf

Steps: 
1. In a system shell (e.g., Windows Powershell, NOT in a python shell) and in the directory that contains bartlein_sub.f90, run the following replacing "mingw32" with the name of your fortan compiler:
    > f2py -c -m bartlein_calendar bartlein_sub.f90 -- compiler=mingw32   
2. Check that a python extension module called "bartlein_calendar" was created. If so, you are ready to go to your python environment and run the calendar correction.
3. In python, run the calendar correction code by first importing the module you just made:
    > import bartlein_calendar  
4. Then call the cal_adjust_pmip function in your python code, for example as in spatialCC.py:
    > bartlein_calendar.cal_adjust_pmip("tas",-999.,-6000.,-6000.,1.,1.,ntime/12,"365_day",model_data[x,y,:].expand_dims(("lat","lon")),ntime)  

   The arguments of this function are:   
   
    - <b>variable name:</b> Use the variable name as found in the title or metadata of the netcdf file, in this example "tas". The calendar correction code will ensure that some relevant variables will not have unphysical negative values.  
    - <b>missing value:</b> In this example, the model output does not have any missing values and this argument is arbitrarily set to -999.  
    - <b>end age BP:</b> In this example, we are correcting a mid Holocene time slice simulation (6000 years BP). Note the use of negative for years before present.  
    - <b>begin age BP:</b> For a time slice simulation, this will have the same value as end age BP since time is not progressing. When correcting a transient time series, this will have a value different from end age BP.  
    - <b>begin year CE:</b> This is relevant only when calendar type (see below) is "gregorian" or "proleptic gregorian" meaning that the model has leap years. You can find calendar type and begin year CE in the netcdf metadata of the climate model file. In this example, the model has a constant 365-day year (no leap years) and this argument is arbitrarily set to 1.  
    - <b>agestep:</b> Interval between age calculations. The code expects a time series consisting of monthly averages. The agestep argument specifies how much time advances for each set of 12 monthly averages. In this example, the monthly means are for one year and agestep is equal to 1. If, however, you are correcting decadal-averaged monthly-mean transient model files, this value would be 10.  
    - <b>nsimyears:</b> Number of years in the model simulation to calendar convert. In the example above, the variable ntime is the length of the time dimension of the model input and ntime/12 is the number of years for the model input.  
    - <b>calendar type:</b> Different climate models use different calendar types, you can find this info in the metadata of the model netcdf file. For this example, it is "365_day", or a constant 365-day year.   
    - <b>input array reshaped to proper size:</b> In this example, the input array was read from the netcdf file and named model_data.  
    - <b>ntime:</b> Total number of months being converted. In this example, the variable ntime is the length of the time dimension of the model input.  
