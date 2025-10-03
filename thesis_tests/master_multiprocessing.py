from multiprocessing import Pool
import subprocess

N_PROCESSES = 2

pool = Pool(N_PROCESSES)

def run_dap_metrics(index: int) -> None:
    output_folder = f"output_{index:03}"
    subprocess.run(["python", "dap_metrics.py", output_folder])

pool.map(run_dap_metrics, list(range(10)))