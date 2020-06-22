"""
average the Fst values from permutation results
"""
import glob
import numpy as np
gold = '/Users/ryan/Documents/Ryan_workplace/DelBay19/11_permutation/obs.fst'

g_arr = []
with open(gold, 'r') as f:
    for l in f.readlines():
        l = l.strip('\n')
        g_arr.append((l.split()[0]+'_'+l.split()[1], float(l.split()[4])))

files = glob.glob('/Volumes/cornell/DelBay19_Hopper/permutation/py_06182020/*.fst')



assert len(files)==999

cnt_array = [0]*291145

for filename in files:

    with open(filename, 'r') as f:
        for i, l in enumerate(f.readlines()):
            l = l.strip('\n')
            snp = l.split()[0] + '_' + l.split()[1]
            fst = float(l.split()[4])

            if fst > g_arr[i][1]:
                cnt_array[i] += 1
            
 
with open('out.txt', 'w') as w, open(gold, 'r') as g:

    for l, fst_cnt in zip(g.readlines(), cnt_array):
        #if fst_cnt < 1:
        w.write('{0}\t{1}\t{2:.4f}\n'.format(l.split()[0],l.split()[1],(fst_cnt*1.0)/(len(files)*1.0)))
        # Permutation P-values Should Never Be Zero (see https://genomicsclass.github.io/book/pages/permutation_tests.html)
        #w.write('{0}\t{1}\t{2:.4f}\n'.format(l.split()[0],l.split()[1],(fst_cnt+1.0)/(len(files)+1.0)))
