# -*- coding: utf-8 -*-
"""
Driver for a single assimilation step using LETKF with tempering,
realistic synthetic RHI observations, and KL divergence analysis.
"""

import sys
sys.path.append('common_python/common_functions/')
sys.path.append('common_python/common_modules/')
sys.path.append('common_python/common_letkf/')

import numpy as np
from scipy.stats import entropy
from cletkf_wloc import common_da as cda
from skimage.draw import line

def compute_total_attenuation(qr, qs, qg, rho_air=1.2):
    wc_r = qr * rho_air * 1e3
    wc_s = qs * rho_air * 1e3
    wc_g = qg * rho_air * 1e3
    k_r = 0.05 * wc_r**0.8
    k_s = 0.01 * wc_s**0.7
    k_g = 0.02 * wc_g**0.7
    return k_r + k_s + k_g, wc_r + wc_s + wc_g

def add_noise(Z, sigma=1.0, seed=None):
    if seed is not None:
        np.random.seed(seed)
    noise = np.random.normal(0, sigma, Z.shape)
    Z_noisy = Z + noise
    Z_noisy[np.isnan(Z)] = np.nan
    return Z_noisy

def simulate_rhi_obs(data, truth_member=0, y_fixed=0, dz_km=0.25, radar_location=(0, 0), 
                     angles_deg=np.linspace(5, 60, 30), beamwidth_deg=1.0,
                     apply_attenuation=True, apply_blocking=True, apply_noise=True,
                     sigma_dbz=1.0, block_threshold_wc=15, block_threshold_att=25, seed=None):
    nx, ny, nz, _, _ = data.shape
    truth = data[:, y_fixed, :, truth_member, :]
    Z = np.full((nx, nz), np.nan)
    att_map = np.zeros((nx, nz))

    radar_x, radar_z = radar_location
    theta_beam_rad = np.radians(beamwidth_deg)

    for angle_deg in angles_deg:
        angle_rad = np.radians(angle_deg)
        beam_length = int(np.hypot(nx, nz))  # maximum diagonal length
        dx = np.cos(angle_rad)
        dz = np.sin(angle_rad)

        cumulative = 0.0
        blocked = False

        for step in range(1, beam_length):
            x = int(round(radar_x + dx * step))
            z = int(round(radar_z + dz * step))

            if x < 0 or x >= nx or z < 0 or z >= nz:
                break

            r = np.hypot(x - radar_x, z - radar_z)
            beam_radius = r * theta_beam_rad
            grid_radius = int(round(beam_radius))

            for i in range(-grid_radius, grid_radius + 1):
                for j in range(-grid_radius, grid_radius + 1):
                    xi, zi = x + i, z + j
                    if xi < 0 or xi >= nx or zi < 0 or zi >= nz:
                        continue
                    if (i**2 + j**2) > grid_radius**2:
                        continue

                    qr = truth[xi, zi, 1]
                    qs = truth[xi, zi, 2]
                    qg = truth[xi, zi, 0]
                    tt = truth[xi, zi, 3]
                    pp = truth[xi, zi, 4]

                    if np.any(np.isnan([qr, qs, qg, tt, pp])) or blocked:
                        Z[xi, zi] = np.nan
                        att_map[xi, zi] = cumulative
                        continue

                    Z0 = calc_reflectivity(qr, qs, qg, tt, pp)

                    if apply_attenuation:
                        k_tot, wc_tot = compute_total_attenuation(qr, qs, qg)
                        cumulative += k_tot * dz_km
                    else:
                        k_tot, wc_tot = 0.0, 0.0

                    if apply_blocking and (wc_tot > block_threshold_wc or cumulative > block_threshold_att):
                        Z[xi, zi] = np.nan
                        blocked = True
                    else:
                        Z[xi, zi] = Z0 - cumulative if apply_attenuation else Z0

                    att_map[xi, zi] = cumulative

    if apply_noise:
        Z = add_noise(Z, sigma=sigma_dbz, seed=seed)

    return Z, att_map

def calc_kl_divergence(observed, modeled):
    observed = np.asarray(observed)
    modeled = np.asarray(modeled)

    if observed.size == 0 or modeled.size == 0:
        return np.nan

    bins = np.linspace(min(observed.min(), modeled.min()), max(observed.max(), modeled.max()), 50)
    obs_hist, _ = np.histogram(observed, bins=bins, density=True)
    mod_hist, _ = np.histogram(modeled, bins=bins, density=True)

    obs_prob = obs_hist / np.sum(obs_hist)
    mod_prob = mod_hist / np.sum(mod_hist)

    epsilon = 1e-10
    obs_prob += epsilon
    mod_prob += epsilon

    kl_div = entropy(obs_prob, mod_prob)
    return kl_div

def calc_reflectivity(qr, qs, qg, t, p):
    #Tong and Xue (2006, 2008) and Smith et al. (1975)
    pi = np.pi
    mindbz = -20.0

    rd = 287.05  # gas constant for dry air [J/kg/K]
    ro = p / (rd * t)  # air density

    nor = 8.0e6   # [m^-4]
    nos = 2.0e6
    nog = 4.0e6
    ror = 1000.0  # [kg/m3]
    ros = 100.0
    rog = 913.0
    roi = 917.0
    roo = 1.0
    ki2 = 0.176
    kr2 = 0.930
    pip = pi ** 1.75

    cf = 1.0e18 * 720 / (pip * (nor ** 0.75) * (ror ** 1.75))
    cf2 = 1.0e18 * 720 * ki2 * (ros ** 0.25) / (pip * kr2 * (nos ** 0.75) * (roi ** 2))
    cf3 = 1.0e18 * 720 / (pip * (nos ** 0.75) * (ros ** 1.75))
    cf4 = (1.0e18 * 720 / (pip * (nog ** 0.75) * (rog ** 1.75))) ** 0.95

    zr = cf * ((ro * qr) ** 1.75) if qr > 0.0 else 0.0
    if qs > 0.0:
        if t <= 273.16:
            zs = cf2 * ((ro * qs) ** 1.75)
        else:
            zs = cf3 * ((ro * qs) ** 1.75)
    else:
        zs = 0.0

    zg = cf4 * ((ro * qg) ** 1.6625) if qg > 0.0 else 0.0

    ref = zr + zs + zg
    if ref > 0:
        ref = 10.0 * np.log10(ref)
    else:
        ref = mindbz

    return ref


if __name__ == "__main__":
    #=========================================================
    # Initial Config
    #=========================================================
    root_path = '/home/jorge.gacitua/experimentos/CrossSection_Assimilation/'
    truth_ens_member = 9
    obs_error_std = 5.0
    loc_scales = np.array([5, 5, 5])
    range_ntemp = np.arange(1, 4)
    range_alpha = [1, 2]
    qg_index, qr_index, qs_index, tt_index, pp_index = 0, 1, 2, 3, 4

    # Load data
    input_ens = np.load(f"{root_path}Data/ensemble_cross_sections_39-5S.npz")["cross_sections"]
    tmp_mask = np.zeros(input_ens.shape[3]).astype(bool)
    tmp_mask[truth_ens_member] = True
    true_state = input_ens[:, :, :, tmp_mask, :][:, :, :, 0, :]
    xf = input_ens[:, :, :, ~tmp_mask, :]
    nx, ny, nz, nbv, nvar = xf.shape

    # Generate observations (synthetic truth)

    Z_rhi,att = simulate_rhi_obs(input_ens, truth_member=truth_ens_member, radar_location=(0, 0),
                                 angles_deg=np.linspace(5, 60, 40), beamwidth_deg=1.0,
                                 apply_attenuation=True, apply_blocking=True, apply_noise=True,
                                 sigma_dbz=1.0, seed=42)
    # Select obs points (non-NaN)
    obs_locs = np.argwhere(~np.isnan(Z_rhi))
    obs_loc_x = obs_locs[:, 0].tolist()
    obs_loc_z = obs_locs[:, 1].tolist()
    obs_loc_y = [0] * len(obs_loc_x)
    nobs = len(obs_loc_x)
    yo = Z_rhi[obs_loc_x, obs_loc_z]

    # Forecast obs ensemble (hxf)
    hxf = np.zeros((nobs, nbv))
    for i in range(nobs):
        ox, oy, oz = obs_loc_x[i], obs_loc_y[i], obs_loc_z[i]
        for jj in range(nbv):
            qr = xf[ox, oy, oz, jj, qr_index]
            qg = xf[ox, oy, oz, jj, qg_index]
            qs = xf[ox, oy, oz, jj, qs_index]
            tt = xf[ox, oy, oz, jj, tt_index]
            pp = xf[ox, oy, oz, jj, pp_index]
            hxf[i, jj] = cda.calc_ref(qr, qs, qg, tt, pp)

    obs_error = obs_error_std * np.ones(nobs)

    # Loop over tempering steps
    for NTemp in range_ntemp:
        for Alpha in range_alpha:
            dt = 1.0 / float(NTemp + 1)
            steps = np.exp(Alpha / np.arange(dt, 1.0 - dt/100.0, dt))
            steps = steps / np.sum(steps)
            steps = (1.0 / steps) / np.sum(1.0 / steps)

            xatemp = np.zeros((nx, ny, nz, nbv, nvar, NTemp + 1))
            xatemp[:, :, :, :, :, 0] = xf
            deps = np.zeros((NTemp, nobs))

            for it in range(NTemp):
                dep = yo - np.mean(hxf, axis=1)
                deps[it, :] = dep

                xatemp = np.asfortranarray(xatemp.astype('float32'))
                hxf_f = np.asfortranarray(hxf.astype('float32'))
                dep = dep.astype('float32')
                obs_error_temp = (obs_error * (1.0 / steps[it])).astype('float32')

                xatemp[:, :, :, :, :, it + 1] =cda.simple_letkf_wloc(nx=nx,ny=ny,nz=nz,
                                                                    nbv=nbv,nvar=nvar,nobs=nobs,
                                                                    hxf=hxf,xf=xatemp[:,:,:,:,:,it],
                                                                    dep=dep,ox=obs_loc_x,
                                                                    oy=obs_loc_y ,oz=obs_loc_z,
                                                                    locs=loc_scales,
                                                                    oerr=obs_error_temp
                                                                    ).astype('float32')

            kl_divergence = calc_kl_divergence(yo, np.mean(hxf, axis=1))

            output_file = f"{root_path}/Data/output_Ntemp{NTemp}_alpha{Alpha}_kl.npz"
            np.savez_compressed(output_file,
                                xatemp=xatemp, xf=xf, yo=yo,
                                hxf=hxf, obs_error=obs_error,
                                obs_loc_x=obs_loc_x,
                                obs_loc_y=obs_loc_y,
                                obs_loc_z=obs_loc_z,
                                true_state=true_state,
                                deps=deps,
                                kl_divergence=kl_divergence)

            print(f"Finished NTemp={NTemp}, Alpha={Alpha}, KL={kl_divergence:.3f}")