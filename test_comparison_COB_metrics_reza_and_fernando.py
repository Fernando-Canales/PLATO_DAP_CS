import numpy as np

Reza_results_dir = 'test_results_multiprocessing/'
Fernando_results_dir = '/home/fercho/double-aperture-photometry/test_results/'

reza_metrics = np.load(Reza_results_dir+'targets_P5.npy')
fernando_metrics = np.load(Fernando_results_dir+'targets_P5.npy')

print(reza_metrics.shape)
print(fernando_metrics.shape)

