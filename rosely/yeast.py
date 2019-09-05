'''
Created on Apr 13, 2018

@author: jiangnan
'''
import pickle
from rosely import *
from numpy import *
from multiprocessing.pool import Pool


wt = CountReads('WT_raw.tsv',   sep='\t', has_groups=False, header=None, index_col=0)
for i in [21, 22, 25, 28, 34, 36]: del wt.data['all'][i]
wt.normalize_data(normalize_by='quantile3')
wt_keys = sorted(wt.data['all'])

idx0=0; nrep = 3; idx1 = idx0 + 21
l = len(wt_keys)

def analysis(args):
    n, i = args
    random.seed(i*n+n)
    keys = array(wt_keys)[random.permutation(l)]
    print('n:{} i:{}:'.format(n, i), keys)
    res = count_analysis([wt.subgroup('all', keys[idx0:idx0+n]),
                            wt.subgroup('all', keys[idx1:idx1+n]),],
                            transform='log1', minmean=10, debug=False, with_plot=False, 
                            pre_neutralize=False, fine_tune=False)
    fdr = false_discovery_rate(list(res.ps.values()))
    nf = sum(fdr<0.05)
    fdr = false_discovery_rate(list(res.nps.values))
    nnf = sum(fdr<0.05)
    return (n, i, nf, nnf)

pool = Pool(8)
idxes = [(n, i) for n in range(2, 21) for i in range(100)]
rs = pool.map(analysis, idxes)
nfs = []; nnfs = []
for r in rs: 
    n, i, nf, nnf = r
    if n-2 >= len(nfs): nfs.append([]); nnfs.append([])
    nfs[-1].append(nf); nnfs[-1].append(nnf)

filename = 'WT_vs_WT_bootstrap_no_neutralize_results_nfs_nnfs.pkl'
with open(filename, 'wb') as f:
    pickle.dump((nfs, nnfs), f)


# ko = CountReads('Snf2_raw.tsv', sep='\t', has_groups=False, header=None, index_col=0)
# for i in [6, 13, 25, 35]: del ko.data['all'][i]
# ko.normalize_data(normalize_by='quantile3')
# ko_keys = sorted(ko.data['all'])
