"""
Test for creating a table with all the relevant metrics
for a given target

Fernando July 8.
"""
import numpy as np # type: ignore

results_dir = '/home/fercho/double-aperture-photometry/test_results/'

metrics_nominal_mask = np.load(results_dir+'targets_P5.npy')
metrics_extended_mask = np.load(results_dir+'targets_P5_extended.npy')
metrics_secondary_mask = np.load(results_dir+'targets_P5_secondary.npy')
metrics_eta_nominal_mask = np.load(results_dir+'eta_bt_24_cameras.npy')
metrics_eta_extended_mask = np.load(results_dir+'eta_ext_bt_24_cameras.npy')

#target_id = 1911666
#target_id = 1135002
target_id = 2427
target_ids = metrics_nominal_mask[:, 0]

index_target_id = np.where(target_ids == target_id)[0][0]

contaminant_stars_IDs = metrics_nominal_mask[index_target_id, 189:199]
magnitude_10first_contaminants = metrics_nominal_mask[index_target_id, 169:179]
delta_x_10first_contaminants = metrics_nominal_mask[index_target_id, 199:209]
delta_y_10first_contaminants = metrics_nominal_mask[index_target_id, 209:219]


sprk_nominal_mask_10first = metrics_nominal_mask[index_target_id, 17:27]*1000000
eta_k_nominal_mask_10first = metrics_eta_nominal_mask[index_target_id, :]
eta_k_cob_nominal_mask_10first = metrics_nominal_mask[index_target_id, 46:56]
nsr_1h_24_cameras_nominal_mask = metrics_nominal_mask[index_target_id, 7]
print('nsr_1h_nom:', nsr_1h_24_cameras_nominal_mask)

eta_k_extended_mask_10first = metrics_eta_extended_mask[index_target_id, :]
eta_k_cob_extended_mask_10first = metrics_extended_mask[index_target_id, 45:55]
nsr_1h_24_cameras_extended_mask = metrics_extended_mask[index_target_id, 4]
print('nsr_1h_ext:', nsr_1h_24_cameras_extended_mask)

nsr_1h_24_cameras_secondary_mask = metrics_secondary_mask[index_target_id, 4]
print('nsr_1h_sec:', nsr_1h_24_cameras_secondary_mask)

# From gaia_crossmatcher.py
#contaminant_stars_GAIA_IDs = np.array([5552642370059981312, 5552642404418780928, 5552643121676968448, 5552643121676965504, 5552642365762705536, 5552643125970942208,  5552642400122460032, 5552643121676971392, 5552643125973284096, 5552642404418783488])
#contaminant_stars_GAIA_IDs = np.array([491200766946790400,  5491200762648590848, 5490450178461270912, 5491200797008328832, 5491200801306529792, 5491200766943729280, 5491200831369388416, 5491200797008329088, 5490450178461271168, 5491200797008329344])
contaminant_stars_GAIA_IDs = np.array([5572932551479590656, 5572932517119851776, 5572932650261393792, 5572932512823744256, 5572932585839329152, 5572932448399294976, 5572932654558804480, 5572932581541918336, 5572932826357494016, 5572932650261397120])

# Function to format and write the important metrics to a LaTeX file for A&A template
def write_important_metrics_to_tex(file_path, fernando_metrics_secondary_mask, target_id, decimal_places=3):
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
        file.write("\\multicolumn{11}{c}{\\textbf{R.A. and Decl. [deg]: 101.84629226792525 and -46.173084507722656}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{Number of contaminants: 31}} \\\\\n")
        file.write("\\midrule\n")
        file.write("\\multicolumn{11}{c}{\\textbf{Magnitude: 12.08}} \\\\\n")
        file.write("\\midrule\n")
        
        
        # Table headers with eta_sec
        file.write("\ Contams & \ P mag. & {$\\Delta x_{\\rm CCD}$} & {$\\Delta y_{\\rm CCD}$} & {$\\rm SPR_{k}^{nom}$} & {$\\eta_k^{\\rm nom}$} & {$\\eta_k^{\\rm ext}$} & {$\\eta_{\\rm kmax}^{\\rm sec}$} & {$\\eta_k^{\\rm nom, \\Delta C}$} & {$\\eta_k^{\\rm ext, \\Delta C}$} & {$\\eta_{\\rm kmax}^{\\rm sec, \\Delta C}$} \\\\\n")
        file.write("\\midrule\n")
        
        # Contaminant rows with GAIA IDs
        for i in range(10):
            cont_id = format_str.format(contaminant_stars_GAIA_IDs[i])
            mc = format_str.format(magnitude_10first_contaminants[i])    # magnitude contaminant
            xc = format_str.format(delta_x_10first_contaminants[i])    # delta x contaminant
            yc = format_str.format(delta_y_10first_contaminants[i])    # delta y contaminant
            
            sprk =format_str.format(sprk_nominal_mask_10first[i]) # SPR_k nominal
            eta_k_nom = format_str.format(eta_k_nominal_mask_10first[i])  # eta_k nominal
            eta_k_ext = format_str.format(eta_k_extended_mask_10first[i])  # eta_k extended
            eta_k_nom_cob = format_str.format(eta_k_cob_nominal_mask_10first[i]) # eta_k nominal COB
            eta_k_ext_cob = format_str.format(eta_k_cob_extended_mask_10first[i]) # eta_k extended COB

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
write_important_metrics_to_tex('important_metrics_comparison_artistic.tex', metrics_secondary_mask, target_id=index_target_id, decimal_places=3)
print("Table saved as important_metrics_comparison_artistic.tex")