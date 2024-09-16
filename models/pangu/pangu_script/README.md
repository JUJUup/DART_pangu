# How to make it work

1. Generate netcdf file using ./pangu2bgrid.py. You may need to change some directory or filename based on your setting.

2. Link the output netcdf to ../work/perfect_input.nc

3. Create the synthetic observation by running ./perfect_model_obs. Possibly need to copy the input.nml in pangu_script dir to work dir.

4. Run filter. Use the assimilate.csh script for mpi running if you like.

## have problem on quickbuild?

There is a mkmf.template in this directory. It is tested on NJU server and should work on same envirement setting.