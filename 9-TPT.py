import numpy as np
import pandas as pd
import pyemma.msm as msm
import pyemma
from multiprocessing import Pool
import itertools
import gc
import os
os.chdir('/home/tnguyen/FAH_MSM_4')


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


if os.path.isdir('TPT') == False:
    os.mkdir('TPT')


A = np.arange(2480, 2490)
B = np.arange(2490, 2500)
process = 'A-I_to_I-A'


def count_tpt(A, B, process):
    MSM = pyemma.load('MSM/pyemma/msm_enspara_tprobs.pyemma')
    tpt = msm.tpt(MSM, A, B)
    (paths, pathfluxes) = tpt.pathways(fraction=0.99)
    per_path = []
    of_total = []
    cumflux = 0
    for i in range(len(paths)):
        cumflux += pathfluxes[i]
        per_path.append('%3.3f' % (100.0*pathfluxes[i]/tpt.total_flux) + ' %')
        of_total.append('%3.3f' % (100.0*cumflux/tpt.total_flux) + ' %')
    print(len(paths))
    d = {'Path flux': pathfluxes, '%path': per_path,
         '%of total': of_total, 'Path': paths}
    df = pd.DataFrame(data=d)
    # df.to_csv(f'TPT_0/pathways_{process}.csv')
    np.save(f'TPT_0/pathways_{process}.npy', paths)
    np.save(f'TPT_0/pathfluxes_{process}.npy', pathfluxes)
    np.save(f'TPT_0/weight_of_paths_{process}.npy', per_path)
    Fsub = tpt.major_flux(fraction=0.99)
    np.save(f'TPT_0/fluxsub_{process}.npy', Fsub)


count_tpt(A, B, process)

print(f'Done')
