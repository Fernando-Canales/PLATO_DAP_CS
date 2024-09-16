"""
Test for creating a table with all the relevant metrics
for a given target

Fernando July 8.
"""
import numpy as np # type: ignore

target_id = 117
Reza_results_dir = '/home/fercho/double-aperture-photometry-tests-directory/test_results_multiprocessing/'
Fernando_results_dir = '/home/fercho/double-aperture-photometry/test_results/'

reza_metrics_nominal_mask = np.load(Reza_results_dir+'targets_P5.npy')
fernando_metrics_nominal_mask = np.load(Fernando_results_dir+'targets_P5.npy')
reza_metrics_extended_mask = np.load(Reza_results_dir+'targets_P5_extended.npy')
fernando_metrics_extended_mask = np.load(Fernando_results_dir+'targets_P5_extended.npy')
fernando_metrics_secondary_mask = np.load(Fernando_results_dir+'targets_P5_secondary.npy')
reza_metrics_eta_nominal_mask = np.load(Reza_results_dir+'eta_bt_24_cameras.npy')
fernando_metrics_eta_nominal_mask = np.load(Fernando_results_dir+'eta_bt_24_cameras.npy')
reza_metrics_eta_extended_mask = np.load(Reza_results_dir+'eta_ext_bt_24_cameras.npy')
fernando_metrics_eta_extended_mask = np.load(Fernando_results_dir+'eta_ext_bt_24_cameras.npy')
reza_metrics_eta_cob_nominal_mask = np.load(Reza_results_dir+'eta_cob_bt_24_cameras.npy')
reza_metrics_eta_cob_extended_mask = np.load(Reza_results_dir+'eta_ext_cob_bt_24_cameras.npy')

# Now I write a dictionary with all the corresponding metrics
important_metrics = {
   'Target ID': fernando_metrics_nominal_mask[target_id, 0],
   'Target Magnitude': fernando_metrics_nominal_mask[target_id, 1],
   'Number of contaminants surrounding the target': fernando_metrics_nominal_mask[target_id, 2],
   'Magnitude of the most significant contaminant': fernando_metrics_nominal_mask[target_id, 3],
   'Distance from the target t.t.m. sign. contaminant': fernando_metrics_nominal_mask[target_id, 4],
   'SPRk_10first': fernando_metrics_nominal_mask[target_id, 17:27]*1000000,
   'eta_sec': fernando_metrics_secondary_mask[target_id, 6],
   'eta_sec_cob': fernando_metrics_secondary_mask[target_id, 9],
   'eta_cob_10first_contaminants': fernando_metrics_nominal_mask[target_id, 46:56],
   'abs_cob_shift_10first': fernando_metrics_nominal_mask[target_id, 66:76],
   'eta_10first_contaminants': fernando_metrics_eta_nominal_mask[target_id, :],
   'eta_ext_10first_contaminants': fernando_metrics_eta_extended_mask[target_id, :],
   'eta_ext_cob_10first_contaminants': fernando_metrics_extended_mask[target_id, 45:55],
   'IDs_from_10first_contaminants': fernando_metrics_nominal_mask[target_id, 189:199],
   'delta_x_contaminant_star': fernando_metrics_nominal_mask[target_id, 199:209],
   'delta_y_contaminant_star': fernando_metrics_nominal_mask[target_id, 209:219],
   'magnitude_contaminant_star': fernando_metrics_nominal_mask[target_id, 169:179]
}

# Function to format and write the important metrics to a LaTeX file for A&A template
def write_important_metrics_to_tex(file_path, important_metrics, decimal_places=3):
    format_str = f"{{:.{decimal_places}f}}" 
    with open(file_path, 'w') as file:
        file.write("\\documentclass{article}\n")
        file.write("\\usepackage{booktabs}\n")   # For cleaner table lines
        file.write("\\usepackage{tabularx}\n")   # For responsive table layout
        file.write("\\usepackage{siunitx}\n")    # For aligning numbers in the table
        file.write("\\usepackage{caption}\n")    # For better captioning
        file.write("\\captionsetup{font=small,labelfont=bf}\n")
        file.write("\\begin{document}\n\n")
        
        file.write("\\begin{table*}[htbp]\n")
        file.write("\\centering\n")
        file.write("\\caption{Comparison of Important Metrics}\n")
        file.write("\\label{tab:important_metrics}\n")
        
        # Defining the table columns with alignment for numbers
        file.write("\\begin{tabularx}{\\textwidth}{|X|S[table-format=2.2]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|}\n")
        
        # Top rule
        file.write("\\toprule\n")
        file.write("\\multicolumn{10}{c}{\\textbf{Target GAIA ID DR3: 5561175130047431680}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{10}{c}{\\textbf{R.A. and Decl. [deg]: 107.1353823470124 and -42.00714666593229}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{10}{c}{\\textbf{Number of contaminants: 141}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{10}{c}{\\textbf{Magnitude: 10.02}} \\\\\n")
        file.write("\\midrule\n")
        
        # Table headers with eta_sec
        file.write("Contams. & \\textbf{P mag.} & {$\\Delta x_{\\rm CCD}$} & {$\\Delta y_{\\rm CCD}$} & {$\\rm SPR_{k}^{nom}$ [ppm]} & {$\\eta_k^{\\rm nom}$} & {$\\eta_k^{\\rm ext}$} & {$\\eta_{\\rm kmax}^{\\rm sec}$} & {$\\eta_k^{\\rm nom, \\Delta C}$} & {$\\eta_k^{\\rm ext, \\Delta C}$} & {$\\eta_{\\rm kmax}^{\\rm sec, \\Delta C}$} \\\\\n")
        file.write("\\midrule\n")
        
        # Contaminant rows
        for i in range(10):
            mc = format_str.format(important_metrics['magnitude_contaminant_star'][i])
            xc = format_str.format(important_metrics['delta_x_contaminant_star'][i])
            yc = format_str.format(important_metrics['delta_y_contaminant_star'][i])
            sprk = format_str.format(important_metrics['SPRk_10first'][i])
            eta_k_nom = format_str.format(important_metrics['eta_10first_contaminants'][i])
            eta_k_ext = format_str.format(important_metrics['eta_ext_10first_contaminants'][i])
            eta_k_nom_cob = format_str.format(important_metrics['eta_cob_10first_contaminants'][i])
            eta_k_ext_cob = format_str.format(important_metrics['eta_ext_cob_10first_contaminants'][i])

            # Add 'eta_sec' only for the first contaminant
            if i == 0:
                eta_kmax_sec = format_str.format(important_metrics['eta_sec'])
                eta_kmax_sec_cob = format_str.format(important_metrics['eta_sec_cob'])
                file.write(f"\\textbf{{Cont. \\#{i + 1}}} & {mc} & {xc} & {yc} & {sprk} & {eta_k_nom} & {eta_k_ext} & {eta_kmax_sec} & {eta_k_nom_cob} & {eta_k_ext_cob} & {eta_kmax_sec_cob} \\\\\n")
            else:
                # Leave the eta_sec column empty for other rows
                file.write(f"\\textbf{{Cont. \\#{i + 1}}} & {mc} & {xc} & {yc} & {sprk} & {eta_k_nom} & {eta_k_ext} & {eta_kmax_sec} & {eta_k_nom_cob} & {eta_k_ext_cob} & {eta_kmax_sec_cob} \\\\\n")
            file.write("\\midrule\n")
        
        # Bottom rule and end of table
        file.write("\\bottomrule\n")
        file.write("\\end{tabularx}\n")
        file.write("\\end{table*}\n")
        
        file.write("\\end{document}\n")
        
# Call the function to write important metrics to a LaTeX file with rounding to 3 decimal places
write_important_metrics_to_tex('important_metrics_comparison_artistic.tex', important_metrics, decimal_places=3)
print("Table saved as important_metrics_comparison_artistic.tex")