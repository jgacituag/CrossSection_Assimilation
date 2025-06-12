import numpy as np
from netCDF4 import Dataset
from wrf import getvar, vertcross, CoordPair, to_np


# Define ensemble members and variables
ensemble_members = [f"{i:03d}" for i in range(1, 31)]
variables = ["QGRAUP", "QRAIN", "QSNOW", "temp", "pressure","ua","va","wa"]
nvar = len(variables)

# Define cross-section coordinates
lat = -39.5
cross_start = CoordPair(lat=-39.2, lon=-65.5)
cross_end = CoordPair(lat=-39.2, lon=-62)

# Placeholder for the first file to get dimensions
sample_file = f"/home/jorge.gacitua/datosmunin2/EXPERIMENTS_UNWEATHER/DATA/PREVENIR_LOWRESOLUTION_HYDRA_2023121612/HIST/FCST/20231216120000/{ensemble_members[0]}/wrfout_d01_2023-12-16_19:00:00"
wrf_file = Dataset(sample_file)
ht = getvar(wrf_file, "z", timeidx=-1)

# Get dimensions
cross_sample = vertcross(ht, ht, wrfin=wrf_file, start_point=cross_start, end_point=cross_end, latlon=True, meta=True)
nz, nx = cross_sample.shape  # nx: horizontal, nz: vertical levels
ny = 1  # Vertical profile
nbv = len(ensemble_members)  # Number of ensemble members

# Initialize storage array
cross_sections = np.zeros((nx, ny, nz, nbv, nvar))

# Loop over ensemble members
for member_idx, member in enumerate(ensemble_members):
    print(f"Processing member {member_idx+1}/{nbv}")
    filename = f"/home/jorge.gacitua/datosmunin2/EXPERIMENTS_UNWEATHER/DATA/PREVENIR_LOWRESOLUTION_HYDRA_2023121612/HIST/FCST/20231216120000/{member}/wrfout_d01_2023-12-16_19:00:00"
    wrf_file = Dataset(filename)
    ht = getvar(wrf_file, "z", timeidx=-1)
    
    for var_idx, var_name in enumerate(variables):
        var_data = getvar(wrf_file, var_name, timeidx=-1)
        cross_data = vertcross(var_data, ht, wrfin=wrf_file, start_point=cross_start, end_point=cross_end, latlon=True, meta=True)
        cross_sections[:, 0, :, member_idx, var_idx] = to_np(cross_data).transpose()
    
    wrf_file.close()

# Save data to compressed npz file
np.savez_compressed("Data/ensemble_cross_sections_39-2S.npz", cross_sections=cross_sections)