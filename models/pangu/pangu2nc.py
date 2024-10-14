import xarray as xr
import numpy as np
import os

if __name__ == "__main__":
    # TOP DIR and variable
    ens_size=50
    DART_DIR="/share/home/lililei4/haoxing/DART/" # where you install DART
    PANGU_DIR="/scratch/lililei4/haoxing/wrfchem_exps/eastAsia/Pangu/Pangu_24h_noDA_smallpert/" # where the pangu output/input is
    # ----------------------------------------------------------
    time = 2023093000
    for i in range(1,ens_size+1): # this is for ensemble members
    # for i in range(1): # this is for OSSE truth
        print(f"now working on member {i}")
        input_fname = f"{PANGU_DIR}/{time}/advance_ensemble/{i}/"
        # DEFINITION of netcdf variable/dimension
        upper = np.load(os.path.join(input_fname, 'output_upper.npy'))
        surf = np.load(os.path.join(input_fname, 'output_surface.npy'))
        lev=[1000,925,850,700,600,500,400,300,250,200,150,100,50]
        lat=np.linspace(90,-90,721)
        lon=np.linspace(0,359.75,1440)
        # bgrid ask the QTY_U/QTY_V to be on the velocity grid. 
        # For now, just transform U/V to velocity grid,
        # do the assimilation, then transform back. Ugly, but it runs.
        # NOTE: need to to some interpolation after the assimition!
        # lat_v=np.arange(90-0.25/2,-90,-0.25)
        # lon_v=np.arange(0+0.25/2,360,0.25)
        # ----------------------------------------------------------
        ds = xr.Dataset(
            data_vars=dict(
                Z=(["time","member","lev", "lat", "long"], 
                    upper[np.newaxis,np.newaxis,0,:,:,:], 
                    {
                        'long_name': "geopential",
                        'units': "gpm",
                    }),
                Q=(["time","member","lev", "lat", "long"], 
                    upper[np.newaxis,np.newaxis,1,:,:,:], 
                    {
                        'long_name': "specific humidity",
                        'units': "kg/kg",
                    }),
                T=(["time","member","lev", "lat", "long"], 
                    upper[np.newaxis,np.newaxis,2,:,:,:], 
                    {
                        'long_name': "temperature",
                        'units': "degrees Kelvin",
                    }),
                U=(["time","member","lev", "lat", "long"], 
                    upper[np.newaxis,np.newaxis,3,:,:,:], 
                    {
                        'long_name': "zonal wind component",
                        'units': 'm/s',
                    }),
                V=(["time","member","lev", "lat", "long"], 
                    upper[np.newaxis,np.newaxis,4,:,:,:], 
                    {
                        'long_name': "meridional wind component",
                        'units': 'm/s'
                        }),
                MSLP=(["time","member","lat", "long"], 
                    surf[np.newaxis,np.newaxis,0,:,:], 
                    {
                        'long_name': 'surface pressure',
                        'units': 'Pa',
                        'units_long_name': 'pascals'
                    }),
                U10=(["time","member","lat", "long"], 
                    surf[np.newaxis,np.newaxis,1,:,:], 
                    {                    
                        'long_name': "10m zonal wind",
                        'units': 'm/s',
                    }),
                V10=(["time","member","lat", "long"], 
                    surf[np.newaxis,np.newaxis,2,:,:], 
                    {                    
                        'long_name': "10 meridional wind",
                        'units': 'm/s',
                    }),            
                T2M=(["time","member","lat", "long"], 
                    surf[np.newaxis,np.newaxis,3,:,:], 
                    {                    
                        'long_name': "2m temperature",
                        'units': "degrees Kelvin",
                    }),
                # TIME=(["time"],[14.5],{'unit': 'days'}),
                # MemberMetadata=(['member', 'metadatalength'],np.array([1,1]))
            ),
            coords=dict(
                # member=("member",[1]),
                lev=("lev", lev, {
                                    'long_name': "level",
                                    'cartesian_axis': "Z",
                                    'units': "hPa",
                                    }),
                lat=("lat", lat, {
                                        'long_name': 'latitude',
                                        'cartesian_axis': 'Y',
                                        'units': 'degrees_north',
                                        'valid_range': [-90.0, 90.0]
                                    }),
                long=("long", lon, {
                                        'long_name': 'longtitude',
                                        'cartesian_axis': 'X',
                                        'units': 'degrees_east',
                                        'valid_range': [0.0, 360.0]
                                    }),
                # tracers=("tracers", [1], {
                #                         'long_name': 'identifier ',
                #                         })
                time=("time",["2011-01-01_00:00:00"], {
                                        'long_name': 'valid time of the model state',
                                        'axis': "T",
                                        'cartesian_axis': 'T',
                                        'calendar': 'none',
                                        # 'units': 'days'
                                    }),
            ),
            attrs=dict(
                description="Pangu nc file for DART, refer to Bgrid model",
                model="pangu",
                history="Haoxing edited in NCAR",
            )
        )
        # ds=ds.assign(MemberMetadata=temp_ds['MemberMetadata'])
        # ds=ds.assign_coords(time=temp_ds['time'])

        encoding = {var: {'_FillValue': None} for var in ds.variables}
        # print(ds['T'].shape) # just show the structure of netcdf file

        # ds.to_netcdf(f"./{time}/NR_{time}.nc",'w',unlimited_dims='time',encoding=encoding)
        ds.to_netcdf(f"./{time}/prior_{time}_mem{i}.nc",'w',unlimited_dims='time',encoding=encoding)
