"""
Test for comparing the COB results from Reza
and my code.

Fernando June 16.
"""
import numpy as np # type: ignore

Reza_results_dir = '/home/fercho/double-aperture-photometry-tests-directory/test_results_multiprocessing/'
Fernando_results_dir = '/home/fercho/double-aperture-photometry/test_results/'

reza_metrics_nominal_mask = np.load(Reza_results_dir+'targets_P5.npy')
fernando_metrics_nominal_mask = np.load(Fernando_results_dir+'targets_P5.npy')
reza_metrics_extended_mask = np.load(Reza_results_dir+'targets_P5_extended.npy')
fernando_metrics_extended_mask = np.load(Fernando_results_dir+'targets_P5_extended.npy')
reza_metrics_eta_nominal_mask = np.load(Reza_results_dir+'eta_bt_24_cameras.npy')
fernando_metrics_eta_nominal_mask = np.load(Fernando_results_dir+'eta_bt_24_cameras.npy')
reza_metrics_eta_extended_mask = np.load(Reza_results_dir+'eta_ext_bt_24_cameras.npy')
fernando_metrics_eta_extended_mask = np.load(Fernando_results_dir+'eta_ext_bt_24_cameras.npy')
reza_metrics_eta_cob_nominal_mask = np.load(Reza_results_dir+'eta_cob_bt_24_cameras.npy')
reza_metrics_eta_cob_extended_mask = np.load(Reza_results_dir+'eta_ext_cob_bt_24_cameras.npy')


# Now we print some useful stuff
print("Shape of Réza's metrics:", reza_metrics_nominal_mask.shape)
print("Shape of Fernando's metrics:", fernando_metrics_nominal_mask.shape)

# Now I can use dictionaries and write a .txt file with the comparison
nominal_metrics = {
    'Target ID': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 0], 1),
        'Fernando': np.round(fernando_metrics_nominal_mask[0, 0], 1)
    },
    'Centroid_shift_errors': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 42:52], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[0, 56:66], 6)
    },
    'Gamma_values': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 32:42], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[0, 106:116], 6)
    },
    'SPRk_values': {
        'Reza': np.round(reza_metrics_nominal_mask[0, 22:32], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[0, 17:27], 6)
    },
    'eta_COB_values':{
        'Reza': np.round(reza_metrics_eta_cob_nominal_mask[0, :], 6),
        'Fernando': np.round(fernando_metrics_nominal_mask[0, 46:56], 6)
    },
    'eta_flux_values':{
        'Reza': np.round(reza_metrics_eta_nominal_mask[0, :], 6),
        'Fernando': np.round(fernando_metrics_eta_nominal_mask[0, :], 6)
    }
}

extended_metrics = {
    'Target ID': {
        'Reza': np.round(reza_metrics_extended_mask[0, 0], 1),
        'Fernando': np.round(fernando_metrics_extended_mask[0, 0], 1)
    },
    'Centroid_shift_errors': {
        'Reza': np.round(reza_metrics_extended_mask[0, 41:51], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[0, 55:65], 6)
    },
    'Gamma_values': {
        'Reza': np.round(reza_metrics_extended_mask[0, 31:41], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[0, 105:115], 6)
    },
    'SPRk_values': {
        'Reza': np.round(reza_metrics_extended_mask[0, 21:31], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[0, 14:24], 6)
    },
    'eta_COB_values':{
        'Reza': np.round(reza_metrics_eta_cob_extended_mask[0, :], 6),
        'Fernando': np.round(fernando_metrics_extended_mask[0, 45:55], 6)
    },
    'eta_flux_values':{
        'Reza': np.round(reza_metrics_eta_extended_mask[0,:], 6),
        'Fernando': np.round(fernando_metrics_eta_extended_mask[0, :], 6)
    }
}

# Function to format and write the metrics to a nice text file
def write_metrics_to_txt(file_path, nominal_metrics, extended_metrics):
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

# Call the function to write metrics to a text file
write_metrics_to_txt('COB_metrics_comparison.txt', nominal_metrics, extended_metrics)

print("Table saved as metrics_comparison.txt")