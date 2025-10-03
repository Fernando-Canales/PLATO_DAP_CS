
import subprocess

for i in range(10):
    n_targets = (i+1)*1000
    output_folder = f"output_{n_targets:06}_targets"
    subprocess.run(["python", "dap_metrics_dummy.py", output_folder, "--number_of_targets", str(n_targets)])
        