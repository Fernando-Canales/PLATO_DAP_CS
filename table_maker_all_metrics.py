"""
Test for creating a table with all the relevant metrics
for a given target

Fernando July 8.
"""
import numpy as np # type: ignore


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

target_id = 1911666
target_ids = fernando_metrics_nominal_mask[:, 0]

index_target_id = np.where(target_ids == target_id)[0]

contaminant_stars_IDs = np.array([5561175061324980736, 5561178119344664192,  5561175061324982528, 5561175130047426944, 5561175164408520832, 5561175057030410624, 5561175061324982400, 5561175130044433152 , 5561175061327958400, 5561175095687700864])

# Function to format and write the important metrics to a LaTeX file for A&A template
def write_important_metrics_to_tex(file_path, contaminant_stars_IDs, fernando_metrics_nominal_mask, fernando_metrics_ext_mask, 
                                   fernando_metrics_secondary_mask, target_id, decimal_places=3):
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
        file.write("\\begin{tabularx}{\\textwidth}{|S[table-format=1.0]|S[table-format=2.2]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|S[table-format=1.3]|}\n")
        
        # Top rule
        file.write("\\toprule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{Target GAIA ID DR3: 5561175130047431680}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{R.A. and Decl. [deg]: 107.1353823470124 and -42.00714666593229}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{Number of contaminants: 141}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{Magnitude: 10.02}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{Units for $\\rm \Delta x_{\\rm CCD} and \Delta y_{\\rm CCD} $: px.}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{Units for $\\rm SPR_{\\rm k}$: ppm}} \\\\\n")
        file.write("\\midrule\n")
        
        
        # Table headers with eta_sec
        file.write("\ Contams & \ P mag. & {$\\Delta x_{\\rm CCD}$} & {$\\Delta y_{\\rm CCD}$} & {$\\rm SPR_{k}^{nom}$} & {$\\eta_k^{\\rm nom}$} & {$\\eta_k^{\\rm ext}$} & {$\\eta_{\\rm kmax}^{\\rm sec}$} & {$\\eta_k^{\\rm nom, \\Delta C}$} & {$\\eta_k^{\\rm ext, \\Delta C}$} & {$\\eta_{\\rm kmax}^{\\rm sec, \\Delta C}$} \\\\\n")
        file.write("\\midrule\n")
        
        # Contaminant rows with GAIA IDs
        for i in range(10):
            cont_id = format_str.format(contaminant_stars_IDs[i])
            mc = format_str.format(fernando_metrics_nominal_mask[target_id, 3].item())    # magnitude contaminant
            xc = format_str.format(fernando_metrics_nominal_mask[target_id, 4].item())    # delta x contaminant
            yc = format_str.format(fernando_metrics_nominal_mask[target_id, 5].item())    # delta y contaminant
            sprk = format_str.format(fernando_metrics_nominal_mask[target_id, 6].item())  # SPR_k nominal
            eta_k_nom = format_str.format(fernando_metrics_nominal_mask[target_id, 7].item())  # eta_k nominal
            eta_k_ext = format_str.format(fernando_metrics_ext_mask[target_id, 7].item())  # eta_k extended
            eta_k_nom_cob = format_str.format(fernando_metrics_nominal_mask[target_id, 8].item())  # eta_k nominal COB
            eta_k_ext_cob = format_str.format(fernando_metrics_ext_mask[target_id, 8].item())  # eta_k extended COB

            # Add 'eta_sec' only for the first contaminant
            if i == 0:
                eta_kmax_sec = format_str.format(fernando_metrics_secondary_mask[target_id, 6].item())  
                eta_kmax_sec_cob = format_str.format(fernando_metrics_secondary_mask[target_id, 9].item())  
                file.write(f"{cont_id} & {mc} & {xc} & {yc} & {sprk} & {eta_k_nom} & {eta_k_ext} & {eta_kmax_sec} & {eta_k_nom_cob} & {eta_k_ext_cob} & {eta_kmax_sec_cob} \\\\\n")
            else:
                # Leave the eta_sec column empty for other rows
                file.write(f"{cont_id} & {mc} & {xc} & {yc} & {sprk} & {eta_k_nom} & {eta_k_ext} & & {eta_k_nom_cob} & {eta_k_ext_cob} \\\\\n")
            file.write("\\midrule\n")
        
        # Bottom rule and end of table
        file.write("\\bottomrule\n")
        file.write("\\end{tabularx}\n")
        file.write("\\end{table*}\n")
        
        file.write("\\end{document}\n")

# Call the function to write important metrics to a LaTeX file with rounding to 3 decimal places
write_important_metrics_to_tex('important_metrics_comparison_artistic.tex', contaminant_stars_IDs, fernando_metrics_nominal_mask, 
                               fernando_metrics_extended_mask, fernando_metrics_secondary_mask, target_id=index_target_id, decimal_places=3)
print("Table saved as important_metrics_comparison_artistic.tex")