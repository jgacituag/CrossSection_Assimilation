import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Increase font sizes
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=14)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=18)   # fontsize of the figure title


# Define output directory
output_dir = "figures"
os.makedirs(output_dir, exist_ok=True)  # Create folder if it doesn’t exist

path = "./"

# Define cases
cases = ["mean", "below", "above"]

# Variable names
variables = {
    0: 'Q Graupel',
    1: 'Q Rain',
    2: 'Q Snow',
    3: 'Temperature'#,
    #5: 'u',
    #6: 'v',
    #7: 'w'
}

r_loc = 5  # Localization radius

for case in cases:
    print(f"Processing case: {case}")

    # Define file paths
    file_paths = {
        '1 Iterations': f'{path}output_Ntemp1_alpha1_8var_Rloc_5_{case}.npz',
        '2 Iterations': f'{path}output_Ntemp2_alpha1_8var_Rloc_5_{case}.npz',
        '3 Iterations': f'{path}output_Ntemp2_alpha1_8var_Rloc_5_{case}.npz'
    }

    for var_index, var_name in variables.items():
        ncols = len(file_paths)
        fig, axs = plt.subplots(3, ncols, figsize=(6 * ncols, 12), dpi = 300)
        if var_index < 3:
            factor = 1e3
        else:
            factor = 1

        contour_levels = {}  # Store contour levels per row
        vmin_fields = []
        vmax_fields = []
        diff_maxs = []
        err_maxs = []

        # First pass to determine shared levels
        for col_index, (col_name, file_path) in enumerate(file_paths.items()):
            data = np.load(file_path)

            obs_loc_x = data['obs_loc_x'][0]
            obs_loc_z = data['obs_loc_z'][0]

            x_start, x_end = obs_loc_x - 3 * r_loc, obs_loc_x + 3 * r_loc + 1
            z_start, z_end = obs_loc_z - 3 * r_loc, obs_loc_z + 3 * r_loc + 1

            xa = data['xatemp'][x_start:x_end, 0, z_start:z_end, :, :, -1]*factor
            xf = data['xf'][x_start:x_end, 0, z_start:z_end, :, :]*factor
            true_state = data['true_state'][x_start:x_end, 0, z_start:z_end, :]*factor

            mean_xa = np.mean(xa, axis=2)
            mean_xf = np.mean(xf, axis=2)

            forecast_field = mean_xf[:, :, var_index]
            analysis_field = mean_xa[:, :, var_index]
            true_field = true_state[:, :, var_index]

            diff_forecast_analysis = forecast_field - analysis_field
            error_analysis = analysis_field - true_field

            vmin_fields.append(np.nanmin(forecast_field))
            vmax_fields.append(np.nanmax(forecast_field))
            diff_maxs.append(np.nanmax(np.abs(diff_forecast_analysis)))
            err_maxs.append(np.nanmax(np.abs(error_analysis)))

        vmin_field = np.nanmin(vmin_fields)
        vmax_field = np.nanmax(vmax_fields)
        contour_levels['field'] = np.linspace(vmin_field, vmax_field, 10)

        diff_max = np.nanmax(diff_maxs)
        contour_levels['diff'] = np.linspace(-diff_max, diff_max, 20)

        err_max = np.nanmax(err_maxs)
        contour_levels['error'] = np.linspace(-err_max, err_max, 20)

        # Second pass to plot with the determined levels
        for col_index, (col_name, file_path) in enumerate(file_paths.items()):
            data = np.load(file_path)

            obs_loc_x = data['obs_loc_x'][0]
            obs_loc_z = data['obs_loc_z'][0]

            x_start, x_end = obs_loc_x - 3 * r_loc, obs_loc_x + 3 * r_loc + 1
            z_start, z_end = obs_loc_z - 3 * r_loc, obs_loc_z + 3 * r_loc + 1

            xa = data['xatemp'][x_start:x_end, 0, z_start:z_end, :, :, -1]*factor
            xf = data['xf'][x_start:x_end, 0, z_start:z_end, :, :]*factor
            true_state = data['true_state'][x_start:x_end, 0, z_start:z_end, :]*factor

            mean_xa = np.mean(xa, axis=2)
            mean_xf = np.mean(xf, axis=2)

            forecast_field = mean_xf[:, :, var_index]
            analysis_field = mean_xa[:, :, var_index]
            true_field = true_state[:, :, var_index]

            diff_forecast_analysis = forecast_field - analysis_field
            error_analysis = analysis_field - true_field

            # Forecast
            im0 = axs[0, col_index].contourf(forecast_field.T, cmap='Spectral_r', levels=contour_levels['field'])
            fig.colorbar(im0, ax=axs[0, col_index])
            if col_index == 0:
                axs[0, col_index].set_title(f"LETKF \n Forecast")
            elif col_index == 1:
                axs[0, col_index].set_title(f"LETKF T2 \n Forecast")
            else:
                axs[0, col_index].set_title(f"LETKF T3 \n Forecast")

            # Forecast - Analysis Difference
            im1 = axs[1, col_index].contourf(diff_forecast_analysis.T, cmap='RdBu_r', levels=contour_levels['diff'])
            fig.colorbar(im1, ax=axs[1, col_index])
            axs[1, col_index].set_title(f"Forecast - Analysis")

            # Compute RMSE inside the square
            square_x_start, square_x_end = r_loc, 6 * r_loc
            square_z_start, square_z_end = r_loc, 6 * r_loc
            square_error = error_analysis[square_x_start:square_x_end, square_z_start:square_z_end]
            rmse = np.sqrt(np.nanmean(square_error**2))

            # Analysis Error
            im2 = axs[2, col_index].contourf(error_analysis.T, cmap='RdBu_r', levels=contour_levels['error'])
            fig.colorbar(im2, ax=axs[2, col_index])
            axs[2, col_index].set_title(f"Error (Analysis - True) | RMSE: {rmse:.3f}")

            # Grid every 2 units and add square
            for ax in axs[:, col_index]:
                ax.grid(True, which='major', linestyle='--', linewidth=0.6, alpha=0.6)
                ax.set_xticks(np.arange(0, forecast_field.shape[0], 5))
                ax.set_yticks(np.arange(0, forecast_field.shape[1], 5))

                # Scatter point with white border
                ax.scatter(3 * r_loc, 3 * r_loc, color="DarkRed", linewidth=1.5, marker='x', s=100)

                # Square region
                #rect = plt.Rectangle((r_loc, r_loc), 6 * r_loc, 6 * r_loc,
                #                     linewidth=2, edgecolor='black', facecolor='none', linestyle='--', alpha=0.6)
                #ax.add_patch(rect)

        # General title
        fig.suptitle(f"Variable: {var_name} | Case: {case}", fontsize=18)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        # Save figure in `figures/`
        output_filename = f"{output_dir}/{case}_{var_name.replace(' ', '_')}.png"
        plt.savefig(output_filename, dpi=300, bbox_inches="tight")
        print(f"Saved: {output_filename}")

        plt.close(fig)  # Close figure to save memory

########### tables

path = "."

# Define cases
cases = ["mean", "below", "above"]

# Variable names
variables = {
    0: 'Q Graupel',
    1: 'Q Rain',
    2: 'Q Snow',
    3: 'Temperature',
    4: 'Pressure',
    5: 'u',
    6: 'v',
    7: 'w'
}
all_rmse = np.zeros((len(cases), len(variables), 3))
rmse_values = {case: np.zeros((len(variables), 3)) for case in cases}
r_loc = 5  # Localization radius

for case in cases:
    print(f"Processing case: {case}")

    # Define file paths
    file_paths = {
        '1 Iterations': f'{path}output_Ntemp1_alpha1_8var_Rloc_5_{case}.npz',
        '2 Iterations': f'{path}output_Ntemp2_alpha1_8var_Rloc_5_{case}.npz',
        '3 Iterations': f'{path}output_Ntemp2_alpha1_8var_Rloc_5_{case}.npz'
    }

    for var_index, var_name in variables.items():
        if var_index < 3:
            factor = 1e3
        else:
            factor = 1
        # Second pass to plot with the determined levels
        for col_index, (col_name, file_path) in enumerate(file_paths.items()):
            data = np.load(file_path)

            obs_loc_x = data['obs_loc_x'][0]
            obs_loc_z = data['obs_loc_z'][0]

            x_start, x_end = obs_loc_x - 3 * r_loc, obs_loc_x + 3 * r_loc + 1
            z_start, z_end = obs_loc_z - 3 * r_loc, obs_loc_z + 3 * r_loc + 1

            xa = data['xatemp'][x_start:x_end, 0, z_start:z_end, :, :, -1]*factor
            true_state = data['true_state'][x_start:x_end, 0, z_start:z_end, :]*factor

            mean_xa = np.mean(xa, axis=2)
            analysis_field = mean_xa[:, :, var_index]
            true_field = true_state[:, :, var_index]

            error_analysis = analysis_field - true_field

            # Compute RMSE inside the square
            square_x_start, square_x_end = r_loc, 6 * r_loc
            square_z_start, square_z_end = r_loc, 6 * r_loc
            square_error = error_analysis[square_x_start:square_x_end, square_z_start:square_z_end]
            rmse = np.sqrt(np.nanmean(square_error**2))
            all_rmse[cases.index(case), var_index, col_index] = rmse
            rmse_values[case][var_index, col_index] = rmse

variables = ['Q Graupel', 'Q Rain', 'Q Snow', 'Temperature','Pressure', 'u', 'v', 'w']
cases = ["mean", "below", "above"]
iterations = ["1 Iter", "2 Iter", "3 Iter"]

# Simulated RMSE values (Replace this with actual RMSE calculations)
#rmse_values = {case: np.zeros(len(variables), len(iterations)) for case in cases}

# Compute reduction percentages
rmse_reductions = {case: np.zeros((len(variables), 2)) for case in cases}

for case in cases:
    for i in range(len(variables)):
        rmse_reductions[case][i, 0] = ((rmse_values[case][i, 1] - rmse_values[case][i, 0]) / rmse_values[case][i, 0]) * 100
        rmse_reductions[case][i, 1] = ((rmse_values[case][i, 2] - rmse_values[case][i, 0]) / rmse_values[case][i, 0]) * 100

# Creating the structured table
columns = [
    "LETKF", "LETKF T2 %", "LETKF T3 %",
    "LETKF", "LETKF T2 %", "LETKF T3 %",
    "LETKF", "LETKF T2 %", "LETKF T3 %"
]

table_values = np.hstack([
    rmse_values["mean"][:, [0]], rmse_reductions["mean"],
    rmse_values["below"][:, [0]], rmse_reductions["below"],
    rmse_values["above"][:, [0]], rmse_reductions["above"]
])

# Create DataFrame
df = pd.DataFrame(table_values, columns=columns, index=variables)

# Print table in a readable format
print(df)

# Save table as CSV
df.to_csv("rmse_table.csv", index=True)
print("Table saved as rmse_table.csv")