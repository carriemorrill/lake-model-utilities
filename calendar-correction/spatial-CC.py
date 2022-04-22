import numpy as np
import bartlein_calendar
import xarray as xr

filename = 'C:\\Users\\carrie.morrill\\CMIP6\\MH-non-CC\\tas_Amon_NorESM2-LM_midHolocene_r1i1p1f1_gn_210101-220012.nc'
with xr.open_dataset(filename) as ds:
#    print(ds.keys())
    nlat = ds.lat.shape[0]
    nlon = ds.lon.shape[0]
    ntime = ds.time.shape[0]
    hold = ds['tas'].transpose('lat','lon','time')
    for x in range(0,nlat):
       for y in range(0,nlon):
            hold2 = bartlein_calendar.cal_adjust_pmip("tas",-999.,-6000.,-6000.,1.,1.,ntime/12,"365_day",hold[x,y,:].expand_dims(("lat","lon")),ntime)
            ds['tas'][:,x,y] = hold2[0,0,:]
    print(hold2)
    ds.to_netcdf("C:\\Users\\carrie.morrill\\CMIP6\\test.nc")