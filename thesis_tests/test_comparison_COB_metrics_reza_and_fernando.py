"""
Test for comparing the COB results from Reza
and my code.

Fernando June 16.
"""
import numpy as np # type: ignore

Reza_results_dir = '/home/fercho/double-aperture-photometry-tests-directory/test_results_multiprocessing/'
Fernando_results_dir = '/home/fercho/double-aperture-photometry/simulation_results/1000_targets_per_magnitude_bin_fixed_dback_85000ppm_and_td_4hr/'
#Fernando_results_dir = '/home/fercho/double-aperture-photometry/test_results/metrics_comparison/'
reza_metrics_nominal_mask = np.load(Reza_results_dir+'targets_P5.npy')
fernando_metrics_nominal_mask = np.load(Fernando_results_dir+'targets_P5.npy')
reza_metrics_extended_mask = np.load(Reza_results_dir+'targets_P5_extended.npy')
fernando_metrics_extended_mask = np.load(Fernando_results_dir+'targets_P5_extended.npy')
reza_metrics_eta_nominal_mask = np.load(Reza_results_dir+'eta_bt_24_cameras.npy')
#fernando_metrics_eta_nominal_mask = np.load(Fernando_results_dir+'eta_bt_24_cameras.npy')
#reza_metrics_eta_extended_mask = np.load(Reza_results_dir+'eta_ext_bt_24_cameras.npy')
#fernando_metrics_eta_extended_mask = np.load(Fernando_results_dir+'eta_ext_bt_24_cameras.npy')
#reza_metrics_eta_cob_nominal_mask = np.load(Reza_results_dir+'eta_cob_bt_24_cameras.npy')
#reza_metrics_eta_cob_extended_mask = np.load(Reza_results_dir+'eta_ext_cob_bt_24_cameras.npy')
fernando_metrics_secondary_mask = np.load(Fernando_results_dir+'targets_P5_secondary.npy') 
reza_metrics_secondary_mask = np.load(Reza_results_dir+'targets_P5_2ndmask.npy')

# Now we print some useful stuff
print("Shape of Réza's metrics:", reza_metrics_nominal_mask.shape)
print("Shape of Fernando's metrics:", fernando_metrics_nominal_mask.shape)

# Now I can use dictionaries and write a .txt file with the comparison
nominal_metrics = {
    'Target ID': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 0], 1),
        'Fernando': np.round(fernando_metrics_nominal_mask[171, 0], 1)
    },
    'Centroid_shift_errors': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 42:52], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[171, 56:66], 6)
    },
    'Gamma_values': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 32:42], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[171, 106:116], 6)
    },
    'SPRk_values': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 22:32], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[171, 17:27], 6)
    },
    'eta_COB_values':{
 #       'Reza': np.round(reza_metrics_eta_cob_nominal_mask[0, :], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[171, 46:56], 6)
    },
    'eta_flux_values':{
        'Reza': np.round(reza_metrics_eta_nominal_mask[0, :], 6),
 #       'Fernando': np.round(fernando_metrics_eta_nominal_mask[171, :], 6)
    }
}

extended_metrics = {
    'Target ID': {
        'Reza': np.round(reza_metrics_extended_mask[0, 0], 1),
        'Fernando': np.round(fernando_metrics_extended_mask[0, 0], 1)
    },
    'Centroid_shift_errors': {
        'Reza': np.round(reza_metrics_extended_mask[0, 41:51], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[171, 55:65], 6)
    },
    'Gamma_values': {
        'Reza': np.round(reza_metrics_extended_mask[0, 31:41], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[171, 105:115], 6)
    },
    'SPRk_values': {
        'Reza': np.round(reza_metrics_extended_mask[0, 21:31], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[171, 14:24], 6)
    },
    'eta_COB_values':{
  #      'Reza': np.round(reza_metrics_eta_cob_extended_mask[0, :], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[171, 45:55], 6)
    },
    'eta_flux_values':{
  #      'Reza': np.round(reza_metrics_eta_extended_mask[0,:], 6),
  #      'Fernando': np.round(fernando_metrics_eta_extended_mask[171, :], 6)
    }
}

secondary_metrics = {
    'Target ID': {
        'Reza': np.round(reza_metrics_secondary_mask[0, 0], 1),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 0], 1)
    },
    'Centroid_shift_errors': {
        'Reza': np.round(reza_metrics_secondary_mask[0, 10], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 10], 6)
    },
    'Gamma_values': {
        'Reza': np.round(reza_metrics_secondary_mask[0, 13], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 17], 6)
    },
    'SPRk_values': {
        'Reza': np.round(reza_metrics_secondary_mask[0, 12], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 5], 6)
    },
    'abs_COB_values':{
        'Reza': np.round(reza_metrics_secondary_mask[0, 9], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 8], 6)
    },
    'sigma_COB_values':{
        'Reza': np.round(reza_metrics_secondary_mask[0, 10], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 10], 6)
    },    
    'eta_COB_values':{
        'Reza': np.round(reza_metrics_secondary_mask[0, 11], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 9], 6)
    },
    'NSR_1h_24':{
        'Reza': np.round(reza_metrics_secondary_mask[0, 4], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 4], 6)
    },
    'eta_flux_values':{
        'Reza': np.round(reza_metrics_secondary_mask[0, 7], 6),
        'Fernando': np.round(fernando_metrics_secondary_mask[171, 6])
    }
}

# Function to format and write nominal and extended metrics
def write_metrics(file_path, nominal_metrics, extended_metrics):
    with open(file_path, 'w') as file:
        file.write("="*65 + "\n")
        file.write(f"{'Nominal Mask Metrics':<40}{'Reza':<15}{'Fernando':<25}\n")
        file.write("="*65 + "\n")
        file.write(f"{'Target ID':<40}{nominal_metrics['Target ID']['Reza']:<15}{nominal_metrics['Target ID']['Fernando']:<25}\n")
        file.write("-"*65 + "\n")
        
        # Write nominal metrics
        for metric, values in nominal_metrics.items():
            if metric != 'Target ID':
                for i, (r, f) in enumerate(zip(values['Reza'], values['Fernando'])):
                    file.write(f"{metric}[{i}]:{'':<15}{r:<15}{f:<15}\n")
                file.write("-"*65 + "\n")

        file.write("="*65 + "\n")
        file.write(f"{'Extended Mask Metrics':<40}{'Reza':<15}{'Fernando':<25}\n")
        file.write("="*65 + "\n")
        file.write("-"*65 + "\n")

        # Write extended metrics
        for metric, values in extended_metrics.items():
            if metric != 'Target ID':
                for i, (r, f) in enumerate(zip(values['Reza'], values['Fernando'])):
                    file.write(f"{metric}[{i}]:{'':<15}{r:<15}{f:<15}\n")
                file.write("-"*65 + "\n")
        file.write("="*65 + "\n")
        file.write(f"{'Secondary Mask Metrics':<40}{'Reza':<15}{'Fernando':<25}\n")
        file.write("="*65 + "\n")
        file.write("-"*65 + "\n")
        
        for metric, values in secondary_metrics.items():
            reza_value = values['Reza']
            fernando_value = values['Fernando']
            if isinstance(reza_value, np.ndarray):
                reza_value = ', '.join(map(str, reza_value))
            if isinstance(fernando_value, np.ndarray):
                fernando_value = ', '.join(map(str, fernando_value))
            file.write(f"{metric}\t{reza_value}\t{fernando_value}\n")

# Call the functions to write metrics to a text file
write_metrics('COB_metrics_comparison.txt', nominal_metrics, extended_metrics)

print("Table saved as COB_metrics_comparison.txt")