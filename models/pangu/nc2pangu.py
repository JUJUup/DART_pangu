import xarray as xr
import numpy as np
import os

if __name__ == "__main__":
    # TOP DIR and variable
    ens_size=50
    DART_DIR="/share/home/lililei4/haoxing/DART/" # where you install DART
    PANGU_DIR="/scratch/lililei4/haoxing/wrfchem_exps/eastAsia/Pangu/Pangu_24h_noDA_smallpert/" # where the pangu output/input is
    template_fname = f"{DART_DIR}/models/bgrid_solo/work/perfect_input.nc"
    temp_ds=xr.open_dataset(template_fname)
    # ----------------------------------------------------------
    time = 2023092000
    # for i in range(1,ens_size+1): # this is for ensemble members

    for i in range(1): # this is for OSSE truth
        print(f"now working on member {i}")
        fname = f"NR_{time}.nc"
        upper = np.zeros((5,13,721,1440))
        surf = np.zeros((4,721,1440))
        ds = xr.open_dataset(fname)
        for i,item in enumerate("Z,Q,T,U,V".split(",")):
            upper[i,:,:,:] = ds[item].sequence().astype(np.float32) # should be single precision!!!
        for i,item in enumerate("MSLP,U10,V10,T2M".split(",")):
            surf[i,:,:] = ds[item].sequence().astype(np.float32)
        np.save(f"output_upper.npy",upper)
        np.save(f"output_surface.npy",surf)
