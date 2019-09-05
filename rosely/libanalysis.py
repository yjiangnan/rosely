'''
Created on Jul 3, 2017

@author: jiangnan
'''
import os, numpy as np
import scipy.stats as ss
from scipy.special import betaln, betainc
from scipy.integrate import quad # General purpose integration
import matplotlib as mpl
mpl.rc('font',family='Helvetica', weight='bold', size=14)
import pickle, operator, logging
np.seterr(all='ignore')
import pandas as pd
from .ascertainedttest import ascertained_ttest, ttest, mean_median_norm
from .neutralstats import *
from collections import OrderedDict
from .pathwayanalysis import getGeneId

__all__ = ['CountReads', 'freqcdf', 'count_analysis', 'presence_absence_analysis', 'combine_by_gene',
           'prepare_ttest_data', 'Result', 'FreqCDF']

default_cross_contamination_rate = 4.4186e-06; 
# From Luciferase count in failed samples: 
# mean([102, 71, 281, 161, 100, 1761, 122])/(sum(hsc(408,:))+sum(mpp(408,:))+sum(b220(408,:))+sum(cmp(408,:)))
class CountReads(object):
    """Class represent the count data and their basic transformations
    
    seqids: a list of sequence ids (e.g. gene names) in the data
    Important dict members with reference hierarchy:
        data['group1']['sample1']['gene1'] = count for gene1 in sample1
        normdata['group1']['sample1']['gene1'] = relative abundance (ratio) for gene1 in sample1 
                                                 compared to the control in sample1
        
        presences['group1']['gene1'] = number of samples gene1 is present in for group1
        sample_presences['group1']['sample1'] = number of genes that are present in sample1 (of group1)
        sample_sizes['group1'] = number of samples for group1
            **Note:**: 
                1. Presence is defined as true if the count for a gene is above 1 plus the total counts for  
                this gene in all samples multiplied by a cross contamination rate * 3 (e.g. barcode mutation).
                2. sample sizes could be non-integers; if a sample contains much fewer genes than the
                median of all samples, then it is considered a partial sample: min(1, p / (median - std))
        
    """
    def __init__(self, filename, sep='\t', has_groups=True, cross_contamination_rate=default_cross_contamination_rate, 
                 group_sample_sep='.', sample_name_transform=None, is_microarray=False, **kwargs):
        """filename: the path to the data file.
        The data file can either be a pickled tuple of (data, seqids) with file name ending with '.dat',
                          or a tab- or comma-(specified by sep) separated file of following format (header 
                              has one element less than other lines, i.e., no name for seqids in header):
group1.sample1,group1.sample2,group2.sample3,group2.sample4
gene1.1,10,15,0,20
gene1.2,15,25,5,30
...
                          or an Excel file with file name ending with '.xls' or '.xlsx' (sample name and 
                              its count data should align).
        """
        self.cross_contamination_rate = cross_contamination_rate
        self.total_count = {}; self.threshold ={}
        self.norm_samples = {}; self.references = {}; self.extra = {}
        self.filename = filename
        if '.dat' == filename[-4:]: self.read_pickled()
        else: self.read_csv(sep, has_groups=has_groups, group_sample_sep=group_sample_sep, 
                            sample_name_transform=sample_name_transform, is_microarray=is_microarray, **kwargs)
        
    def read_pickled(self):
        with open(os.path.expanduser(self.filename),'rb') as f: data, self.seqids = pickle.load(f)
        self.data = {}
        for group in data: 
            for bc in data[group]: # barcode or sample name
                self.data[group, bc] = data[group][bc]
        self.data = pd.DataFrame(self.data)
        self.data[self.data.isnull()] = 0
        
    def read_csv(self, sep='\t', igroup=1, has_groups=True, group_sample_sep='.', sample_name_transform=None, 
                 is_microarray=False, **kwargs):
        if '.xls' in self.filename[-5:]: 
            csv = pd.read_excel(self.filename, **kwargs)
        else: csv = pd.read_csv(self.filename, sep=sep, **kwargs)
        if is_microarray:
            cleaned = pd.DataFrame()
            keys = list(csv.keys())
            cleaned[keys[0]] = csv[keys[0]]
            for i in range(1, len(keys), 2):
                if 'pval' not in keys[i+1].lower(): 
                    raise Exception('Microarray data should have one column of Detection_Pval for each data column.')
                cleaned[keys[i]] = csv[keys[i]]
                cleaned.loc[csv[keys[i+1]] > 0.05, keys[i]] = np.NaN
            csv = cleaned
        if csv.first_valid_index() == 0: csv = csv.set_index(list(csv.keys())[0])
        self.seqids = list(csv.index); self.data = {}
        for sample in csv.keys():
            splname = sample
            if type(splname) is str: splname = sample.strip()
            if sample_name_transform is not None: splname = sample_name_transform(sample) 
            if has_groups: spl = splname.split(group_sample_sep)
            if has_groups and len(spl) <= 1:
                print('Warning:', sample, 'is not included (saved in self.extra) because it cannot be separated into group and sample.')
                self.extra[splname] = csv[sample]
            else:
                if has_groups: 
                    group = group_sample_sep.join(spl[:igroup])
                    spl   = group_sample_sep.join(spl[igroup:])
                else: group = 'all'; spl = splname
                try: self.data[group, spl] = pd.to_numeric(csv[sample])
                except: logging.exception((sample + ' cannot be parsed in numbers'))
        self.data = pd.DataFrame(self.data)
        if not len(self.data):
            raise Exception('No samples have been found.')
        
    def normalize_data(self, normcutoff, normalize_by='median', with_threshold=False):
        '''Normalize different samples to the same level and calculate important stats such as presences
        normalize_by: either 'median' (of all non-zero count seqids)
                      or 'quantileN' (1 <= N <= 9) to normalize by the mean of the 1/N proportion of data around the median
                      or a seqid (e.g. Luciferase with massive quantity) in the data
                      or a list of seqids whose median will be used for normalization.
         normcutoff: The minimum average value to be considered non-zero and included for calculating median 
                     and quantile during normalization. 
        '''
        self.total_count = self.data.sum(axis=1)
        self.threshold = (1 + self.cross_contamination_rate * self.total_count * 3) * with_threshold
        self.data[self.data.subtract(self.threshold, axis=0) < 0] = 0
        self.shared_presences = {}
        self.sample_presences = (self.data>0).sum(axis=0)
        median_presences = self.sample_presences.groupby(level=0).median()
        std_presences = self.sample_presences.groupby(level=0).std()
        p = self.sample_presences.divide(median_presences - std_presences, level=0)
        self.sample_weights = p * (p<1) + (p>=1)
        self.sample_sizes = self.sample_weights.groupby(level=0).sum()
        self.presences = {}
        for group in self.data.keys().levels[0]:
            sw = self.sample_weights[group]
            self.presences[group] = (self.data[group]>0).multiply(sw).sum(axis=1)
        self.nRefs = nRefs = {}
        nby = normalize_by
        if nby == 'median' or 'quantile' in nby: nby = list(self.seqids)
        allsamples = self.data.mean(axis=1)
        allsamples = allsamples[allsamples > normcutoff]
        if 'quantile' in normalize_by: 
            qt = int(normalize_by[-1]) * 2; l = len(allsamples)
            idx0 = l//2 - l//qt; idx1 = l//2 + l//qt
            if 'upper' in normalize_by: idx0 += l//qt; idx1 += l//qt
            mdall = sorted(allsamples)[idx0 : idx1]
            allsamples = allsamples[(allsamples > mdall[0]) & (allsamples < mdall[-1])]
        for key in self.data: 
            if type(nby) is list:
                sample = self.data[key]
                if 'quantile' in normalize_by: 
                    nRefs[key] = np.mean(sample.loc[allsamples.index])
                else: nRefs[key] = np.median(sample.loc[allsamples.index])
                if nRefs[key]==0: nRefs[key] = np.NaN
            else:
                if normalize_by not in self.data[key]: nRefs[key] = np.NaN
                else: nRefs[key] = max(1, self.data[key][normalize_by])
        self.mean_reference = np.nanmean(list(nRefs.values()))
        self.normdata = self.data / pd.Series(nRefs) * self.mean_reference
        d = self.normdata
        dup = d[d.index.duplicated()]
        if len(dup.index):
            print('Only the first of the duplicated indices will be kept:\n', dup)
            self.normdata = d.groupby(d.index).first()
        if normalize_by in self.seqids: 
            self.seqids.remove(normalize_by)
            self.normdata = self.normdata.loc[self.seqids]
        
    def combine_technical_replicates(self, rename):
        """rename: a function to rename the sample names of technical replicates into independent replicates."""
        for group in self.normdata:
            nd = pd.DataFrame()
            samples = set([rename(k) for k in self.normdata[group].keys()])
            for s in samples:
                sdf = pd.DataFrame()
                for k in self.normdata[group].keys():
                    if rename(k) == s: 
                        sdf[k] = self.normdata[group][k]
                nd[s] = sdf.mean(axis=1)
            self.normdata[group] = nd
    
    def get_data(self, bc, seqid):
        if seqid in self.samples[bc]: return self.samples[bc][seqid]
        else: return 0
        
    def get_threshold(self, seqid):
        if seqid in self.threshold: return self.threshold[seqid]
        else: return 0
        
    def nRefs_as_weights(self, group, samples='all'):
        if samples is 'all': samples = self.normdata[group].keys()
        weights = {}
        ws = [self.nRefs[group, sample] for sample in samples if self.nRefs[group, sample] > 0]
        mw = np.mean(sorted(ws)[len(ws)//2:])
        for sample in samples:
            weights[sample] = min(1, self.nRefs[group, sample] / mw)
        return weights
    
    def subgroup(self, group, samples):
        normdata = {}
        for sample in samples: normdata[sample] = self.normdata[group][sample]
        return normdata
    
    def set_missing_to_0s(self, seqids):
        for seqid in seqids:
            for group in self.data:
                for sample in self.data[group]:
                    if seqid not in self.data[group][sample]: 
                        self.data[group][sample][seqid] = 0
                        if seqid not in self.seqids: self.seqids.append(seqid)
    
    def correlation_heatmap(self, data=None, groups='all', samples='all', vmin=0.8, transform='log', 
                    by_residual=False, with_plot=True, **kwargs):
        if data is None: data = self.normdata
        correlation_heatmap(data, groups=groups, samples=samples, vmin=vmin, 
                    transform=transform, by_residual=by_residual, with_plot=with_plot, **kwargs)
                        
def correlation_heatmap(data, groups='all', samples='all', vmin=0.8, transform='log', 
                    by_residual=False, with_plot=True, paired=False, **kwargs):
    try: 
        if groups == 'all': 
            try: groups = sorted(data.keys().levels[0])
            except: groups = data.keys()
        if samples == 'all': group_samples = [(g, s) for g in groups for s in data[g]]
        elif len(groups) == 1: group_samples = list(zip(groups * len(samples), samples))
        else: group_samples = list(zip(*[[groups[i]] * len(samples[i]) for i in range(len(groups))],
                                  [s for spl in samples for s in spl]))
    except: group_samples = sorted(data.keys())
    try: data = pd.DataFrame({gs:data[gs] for gs in group_samples})
    except: data = pd.DataFrame({(g,s):data[g][s] for g in groups for s in data[g]})
    if transform == 'log': 
        data = np.log(data)
        data[data==-np.Inf] = np.nan
        data = data.dropna()
    if paired:
        data = data[groups[1]] - data[groups[0]]
        group_samples = sorted(data.keys())
    if by_residual:
        data = pd.DataFrame(data.values - data.mean(axis=1).values[:, None], 
                            index=data.index, columns=data.columns)
        vmin = -1
        kwargs['cmap'] ="coolwarm"
        kwargs['center'] = 0
    if type(group_samples[0]) is tuple:
        M = pd.DataFrame(columns=pd.MultiIndex.from_tuples(group_samples), 
                         index=pd.MultiIndex.from_tuples(group_samples))
    else:
        M = pd.DataFrame(columns = group_samples, 
                         index   = group_samples)
    for i, si in enumerate(group_samples):
        M[si][si] = 1
        for sj in group_samples[i+1:]:
            r = np.corrcoef(data[si], data[sj])[0][1]
            M[sj][si] = M[si][sj] = r
    if with_plot:
        import seaborn
        seaborn.heatmap(M.astype('float64'), vmin=vmin, **kwargs)
        mpl.pyplot.xlabel(''); mpl.pyplot.ylabel('')
    return M
    
        
class Result:
    def __init__(self, ps, zs, LFDRthr0, data_name='data', nbins=50, df_fit=5, minr0=0, with_plot=True,
                 base_pop=1, fine_tune=True, neutralize=True, global_changes=False, **kwargs):
        self.ps = pd.Series(ps); self.zs = pd.Series(zs)
        self.base_pop = base_pop
        self.results = OrderedDict()
        if neutralize:
            self.neutralize_p_values(LFDRthr0, data_name=data_name, nbins=nbins, df_fit=df_fit, with_plot=with_plot,
                                     minr0=minr0, fine_tune=fine_tune, global_changes=global_changes, **kwargs)
        else: 
            self.results['LFDR'] = locfdr(zs=zs, p0s=[0])
            self.results['z-score'] = self.zs
            self.results['p-value'] = self.ps
            self.results = pd.DataFrame(self.results).sort_values(by='p-value')
            
        
    def neutralize_p_values(self, LFDRthr0=0.5, data_name='data', nbins=50, df_fit=5, minr0=0, 
                            fine_tune=True, with_plot=True, global_changes=False,
                            poisreg_method='SLSQP', **kwargs):
        self.nps, self.LFDR, self.pop = neup(self.ps, LFDRthr0=LFDRthr0, data_name=data_name,
                    nbins=nbins, df_fit=df_fit, minr0=minr0, base_pop=self.base_pop,
                    fine_tune=fine_tune, with_plot=with_plot, global_changes=global_changes,
                    poisreg_method=poisreg_method, **kwargs)
        self.nzs = np.sign(self.zs) * abs(ss.norm.ppf(self.nps/2))
        self.results['LFDR'] = self.LFDR
        self.results['Controlled z-score'] = pd.Series(self.nzs)
        self.results['Controlled p-value'] = self.nps
        self.results = pd.DataFrame(self.results).sort_values(by='Controlled p-value')
        
    def update_gene_ids(self, species='mouse', taxid=None):
        ids = getGeneId(self.results.index, species=species, taxid=taxid)
        rr = self.results.copy()
        rr['geneid'] = ids['geneid']
        rr['symbol'] = ids['symbol']
        self.results['symbol'] = ids['symbol']
        return rr.set_index('geneid').dropna()
        
    def plot_variance(self, sample_index='all', xname='mean', uselogV=False, figsize=None):
        if 'vbs' not in self.__dict__: 
            raise Exception('Variable vbs does not exist. Run count_analysis with debug=True to calculate it.')
        from matplotlib.pyplot import plot, xlabel, ylabel, subplot, legend, figure
        def plotx(xname, vbs):
            x = vbs[xname]
            logVs = vbs['logVs']
            ElogVs = vbs['ElogVs']
            S2all = vbs['S2all']; VlogV0s = vbs['VlogV0s']; Vadj = vbs['Vadj']
            plot(x, logVs , '.', label='all genes')
            plot(x, ElogVs, '.', label='mean') 
            plot(x, ElogVs + S2all**0.5, '*r', label='SD')
            plot(x, ElogVs - S2all**0.5, '*r')
            plot(x[:len(Vadj)], ElogVs + Vadj**0.5, '+c', label='Sampling SD')
            plot(x[:len(Vadj)], ElogVs - Vadj**0.5, '+c')
            plot(x, ElogVs + VlogV0s ** 0.5, '.m', label='prior gene SD')
            plot(x, ElogVs - VlogV0s ** 0.5, '.m') 
            xlabel(xname, fontsize=20, fontweight='bold')
            if uselogV: ylabel('log variance', fontsize=20, fontweight='bold')
            else:
                a = vbs['a']
                ylabel('variance$ ^{{{:.3f}}} $'.format(a), fontsize=20, fontweight='bold')
            legend()
        if sample_index is 'all': sample_index = self.vbs['idxes']
        elif type(sample_index) is int: sample_index = [sample_index]
        if figsize is None: figsize = (2 + 6*len(sample_index), 5)
        figure(figsize=figsize)
        for i, idx in enumerate(sample_index):
            vbs = self.vbs[idx]
            for v in vbs: 
                if v != 'a': vbs[v] = pd.Series(vbs[v])
            subplot(101+10*len(sample_index) + i)
            plotx(xname, vbs)
#         if np.nanmin(vbs['n']) < np.nanmax(vbs['n']): 
#             subplot(122)
#             plotx('n', vbs)

def freqcdf(N1, X, N2, Y, r=1):
    '''Calculate the p-value (2-tailed) for the difference between frequencies X and Y in total trails of N1 and N2, 
    respectively, given a dilution ratio r (expected concentration of Y divided by X), using Bayesian inference.
    r is not simply the expected value of Y divided by that of X, but should be 
                log(1 - E(Y)/N2) / log(1 - E(X)/N1)
        i.e., when E(Y) is close to N1, the concentration of Y should be close to infinity.
        Consider infecting N cells with y viruses with a chance of p, i.e., yp viruses would enter cells. 
        Obviously, yp could increase to infinity but the total number Y of cells being infected can at 
        most be N. Then, we say the concentration of Y is yp for simplicity: yp = -log(1 - E(Y)/N)
    '''
    lnb1 = betaln(X+1, N1-X+1)
    lnb2 = betaln(Y+1, N2-Y+1)
    Et1 = 1 - np.exp(betaln(X+1, N1-X+1+r  ) - betaln(X+1, N1-X+1)) # exp(log beta) increases numeric stability 
    Et2 = 1 - np.exp(betaln(Y+1, N2-Y+1+1/r) - betaln(Y+1, N2-Y+1)) # exp(log beta) increases numeric stability
    s = np.sign((Y/N2 - Et1) - (X/N1 - Et2))
    def integrand(x, N1, X, N2, Y, r, s, lnb):
        # ss.binom.cdf(X, N, p) == betainc(N-X, X+1, 1-p), but betainc accepts non-integer parameters.
        # However, betainc(0, ., .) and betainc(., 0, .) are nan because gamma(0) == 0
        if s <= 0:  # Decreased, accumulate from 0 to Y 
            b = betainc(max(1e-19, N2 - Y),  Y+1,  (1-x)**r)
        else:       # Increased, accumulate from Y to N2, or 1 - cdf(from 0 to Y-1)
            b = 1 - betainc(N2 - (Y-1),  max(1e-19, Y),  (1-x)**r)
        return b * np.exp( X * np.log(x) + (N1 - X) * np.log(1 - x) - lnb)
    points = [(X/N1 + Y/N2)/2, Et1, Et2]
    I1      = quad(integrand, 0, 1, args=(N1, X, N2, Y, r  ,  s, lnb1), points=points)[0]
    try: I2 = quad(integrand, 0, 1, args=(N2, Y, N1, X, 1/r, -s, lnb2), points=points)[0]
    except: print('N2, Y, N1, X, r, Et1, Et2 = ', N2, Y, N1, X, r, Et1, Et2, sep=', ')
    p = I1 + I2
    if p > 1: p = 2 - p
    pval, z = p, -ss.norm.ppf(p/2) * s
    return pval, z

class FreqCDF:
    name = 'rosely.freqcdf'
    def __init__(self): print('Switching to', self.name, 'for p-value calculation.')
    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        return freqcdf(study_n, study_count, pop_n, pop_count, 1)[0]

def presence_absence_analysis(presences1, N1, presences2, N2, controls = None, min_count=0, shared_presences=None, 
                              LFDRthr0=0.5, minr0=0, data_name='data', nbins=50, df_fit=5, neutralize=True, **kwargs):
    """Compare the difference between presences of two groups for each gene
    Args:
    presences1: dict: {'gene A':number of replicates gene A is present in for group1, ...}
    N1:         total number of samples for group1
    presences2: dict: {'gene A':number of replicates gene A is present in for group2, ...}
    N2:         total number of samples for group2
    controls:   list of genes to control for overall differences between group1 and group2 
    min_count:  minimal total number of presences and absences required for each gene. 
                Failed genes would be treated as missing data (set to NaN for the results)
    Other parameters are passed to neup and locfdr
    
    Returns:
    An object with the following members:
        ps:    The raw p-values
        zs:    The raw z-scores
        nps:   The neutralized p-values
        nzs:   The neutralized z-scores
        LFDR:  Local False Discovery Rate calculated from netralized p-values nps
        pop:   power of p-values
    """
    zs = {}; ps = {}
    seqids = set(list(presences1.keys()) + list(presences2.keys()))
    if shared_presences is None: shared_presences = {}
    for seqid in seqids:
        if seqid not in presences1: presences1[seqid] = 0
        if seqid not in presences2: presences2[seqid] = 0
        if seqid not in shared_presences: shared_presences[seqid] = 0
    ctrls = controls
    if controls is None: ctrls = seqids
    nc1 = np.array(sorted([presences1[seqid] - shared_presences[seqid] for seqid in ctrls]))
    nc2 = np.array(sorted([presences2[seqid] - shared_presences[seqid] for seqid in ctrls]))
    valid = np.logical_and(np.logical_and(nc1 > 0, nc2 > 0), np.logical_and(nc1 < N1*0.999, nc2 < N2*0.999))
    vnc1 = nc1[valid]; vnc2 = nc2[valid]
    l = len(vnc1)
    if l < 10:
        print('Two few constructs have >0 but not full sample counts. valid/total: {}/{}'.format(l, len(nc1))) 
        raise Exception('')
    mean1 = np.mean(vnc1[l//4:l//4*3]); mean2 = np.mean(vnc2[l//4:l//4*3])
    mdn = np.mean(list(shared_presences.values()))
    r = np.log(1 - mean2/(N2-mdn)) / np.log(1 - mean1/(N1-mdn))
    print('N1: {:.2f}  N2: {:.2f};  25~75% mean1: {:.2f}  25~75% mean2: {:.2f} in {} controls;   shared: {:.2f}  MOI ratio: {:.2f}'.format(
        N1, N2, mean1, mean2, l, mdn, r))
    for seqid in seqids:
        dn = shared_presences[seqid]
        n1 = presences1[seqid]; n2 = presences2[seqid]; z = 0
        if n1 + n2 - dn*2 >= min_count and (N1-n1) + (N2-n2) >= min_count: 
            pval, z = freqcdf(N1-dn, n1-dn, N2-dn, n2-dn, r)
        else: z = np.NaN; pval = np.NaN
        zs[seqid] = z
        ps[seqid] = pval
    res = Result(ps, zs, LFDRthr0, data_name=data_name, nbins=nbins, df_fit=df_fit, minr0=minr0, 
                 neutralize=neutralize, **kwargs)
    change = {}
    for seqid in ps: 
        change[seqid] = '{:.3f} --> {:.3f}'.format(
                        presences1[seqid], presences2[seqid])
    res.results['presence change'] = pd.Series(change)
    return res

def prepare_ttest_data(normdata_list, transform, minmean, minn, sample_weights, paired, drop_missmatch):
    a = None
    if transform in ['mmm', 'MMM', 'minimize_mean_median']:
        alldata = {}
        for i, data in enumerate(normdata_list): 
            for sample in sorted(data): 
                alldata[(i, sample)] = data[sample]
        alldata = pd.DataFrame(alldata)
        ms = alldata.mean(axis=1)
        ms = ms[ms > minmean]
        if not np.isnan(minmean): ms += minmean * 0.5
        _, a = mean_median_norm(ms)
    for i, data in enumerate(normdata_list): 
        data = pd.DataFrame(data.copy())
        if not np.isnan(minmean): data = data[data.mean(axis=1) > minmean]
        if transform == 'log': 
            data[data<=0] = np.nan
            data = np.log2(data)
        elif type(transform) is str and transform[:3] == 'log':
            n0 = float(transform[3:])
            if '.' not in transform[3:]: n0 = int(n0)
            if str(n0) != transform[3:]: 
                raise Exception('The string following log is not a well formated float number')
            data[data<=-n0] = np.nan
            data = np.log2(n0# * (np.random.uniform(size=data.shape) * 0.25 + 0.75)
                            + data)
        elif type(transform) in [float, int]:
            data = data ** transform
        elif transform in ['mmm', 'MMM', 'minimize_mean_median']:
            if not np.isnan(minmean): data += minmean * 0.5
            data = data ** a
        else: raise Exception('transform has to be either log or a float number (power)')
        normdata_list[i] = data
    ttestdata = {}; weights = {}
    for i, normdata in enumerate(normdata_list):
        for sample in sorted(normdata): 
            if sample_weights: 
                weights[(i, sample)] = sample_weights[i][sample]
            ttestdata[(i, sample)] = normdata[sample]
    ttestdata = pd.DataFrame(ttestdata)
    good = ttestdata.notnull().sum(axis=1) >= minn
    for i in ttestdata.keys().levels[0]:
        good = (good & (ttestdata[i].notnull().sum(axis=1) >= minn)) # Remove rows if numbers of NaN samples in a group is < minn
    ttestdata = ttestdata[good]
    if paired:
        valid_keys = [k for k in ttestdata.keys() 
                      if  k[1] in normdata_list[0].keys()
                      and k[1] in normdata_list[1].keys()]
        if len(valid_keys) < len(ttestdata.keys()):
            print('Some sample names are not paired and removed:')
            print('valid names:', valid_keys)
            print('all   names:', list(ttestdata.keys())) 
            print('unpaired names:', set(ttestdata.keys()).difference(valid_keys))
            if drop_missmatch: 
                ttestdata = ttestdata[valid_keys]
            else: raise Exception('Paired test requires the input column names to match.')
            'For paired test, set data to NaN if any sample for the seqid is NaN'
        s = 0
        for i in range(len(normdata_list)): s = s + normdata_list[i]
        good = (s.notnull()).sum(axis=1) > 1
        ttestdata = ttestdata.loc[good]
    nv = (~ttestdata.isnull()).sum(axis=0)
    col = nv[(nv>5) & (nv > ttestdata.shape[0] * 0.02)].index
    ttestdata = ttestdata[col]
    if sample_weights is None: weights = None
    else: weights = pd.Series(weights)[col]
    return pd.DataFrame(ttestdata), weights, a

def count_analysis(normdata_list, transform=1., minmean=np.nan, controls=None, sample_weights=None, weighted=True, minn=1.25,
                           paired=False, drop_missmatch=True, equalV=False,  debug=False, method='ascertained_ttest', 
                           pre_neutralize=True, with_plot=True, 
                           span=None, LFDRthr0=0.5, minr0=0, fine_tune=True, neutralize=True, global_changes=False,
                           data_name='data', nbins=50, df_fit=5, do_SVA=False, nSV=2, penaltyRatioToSV=0.2, 
                           ridge=True, normalize_min_pval=0.1, parallel=True, **kwargs):
    """Analyze the counts by ascertained_ttest
    
    Args:
    transform: how the data should be transformed for normalization,
               either 'log' for log-normal data (0s will be excluded),
               or     'logN for log-normal data including 0s (excluding data <= -1, x --> log(N+x)),
               or     a power (0~1) for intermediate between log-normal and linear data, i.e., data**power should
                          be normally distributed. You can determine the power by function mean_median_norm(data_list)
    minmean: the minimum mean values for data to be included for analysis, default 5 if transform is log1, or NaN else.
    Other args are passed to ascertained_ttest (or ttest),  neup, and locfdr.
    
    Returns:
    An object with the following members:
        ps:    The raw p-values
        zs:    The raw z-scores
        nps:   The neutralized p-values
        nzs:   The neutralized z-scores
        LFDR:  Local False Discovery Rate calculated from netralized p-values nps
        pop:   power of p-values
    """
#     print('Sample sizes:', [len(d.keys()) for d in normdata_list])
    if np.isnan(minmean) and transform == 'log1': minmean = 10
    ttestdata, weights, a = prepare_ttest_data(normdata_list, transform, minmean, 
                                               minn, sample_weights, paired=paired, drop_missmatch=drop_missmatch)
    if with_plot==True and a: print('Data MMM normalization power:', a)
    test = ttest
    if not weighted: weights = False
    testkwargs = {
        'controls':controls, 'paired': paired, 'weights': weights, 'equalV': equalV,
        'do_SVA': do_SVA, 'nSV': nSV, 'penaltyRatioToSV': penaltyRatioToSV,
        'normalize_min_pval': normalize_min_pval,
        'ridge': ridge, 'parallel': parallel
        }
    if method=='ascertained_ttest':
        test = ascertained_ttest
        testkwargs = {**testkwargs, 'span': span, 'pre_neutralize':pre_neutralize}
    ps, zs, vbs = test(ttestdata, **testkwargs)
    dxs = vbs['dx']
    base_pop = 1
    if 'base_pop' in vbs: base_pop = vbs['base_pop'][0]
    res = Result(ps, zs, LFDRthr0, data_name=data_name, nbins=nbins, df_fit=df_fit, minr0=minr0, with_plot=with_plot,
                 fine_tune=fine_tune, global_changes=global_changes, neutralize=neutralize, base_pop=base_pop, **kwargs)
    if debug: res.vbs = vbs
    nm = 'dx'
    if 'ttest_pops' in vbs: res.ttest_pops = vbs['ttest_pops']
    if type(transform) is str and 'log' in transform: nm = 'log2 fold change'
    res.results[nm] = pd.Series(dxs)
    if with_plot: print('shape (dimensions) of results:', res.results.shape)
    return res

def combine_by_gene(zs, sep='.', LFDR0=0.3, p0s=[0, ], nbins=50, df_fit=5, power_of_pval=1, 
                    data_name='data', method='Stouffer and Fisher'):
    """Combine different z-scores of the same gene into a single z-score and LFDR
    
    Args:
    zs:     a dict of z-scores of different seqids
    sep:    The separation character in each seqid to identify the gene name
            e.g., geneX can have two seqids: geneX.1 and geneX.2
    Other parameters are passed to locfdr
    
    Returns:
        * Combined LFDR by Stouffer's z-score method (z-scores are added and normalized)
        * Combined LFDR by Fisher's p-value method (signs of z-scores do not matter)
    """
    genes = []
    for seqid in zs.index: genes.append(seqid.split(sep)[0])
    genes = list(set(genes))
    gene_zs = {}; gene_chi2ps = {}; gzs = {}; gene_targets = {}
    for seqid in zs.index:
        g = seqid.split(sep)[0]; p = zs[seqid]
        if not np.isnan(p):
            if g not in gzs: gene_targets[g] = []; gzs[g] = []
            gzs[g].append(p)
            gene_targets[g].append(seqid)
    for g in gzs:
        if len(gzs[g]) > 1: 
            gene_zs[g] = np.sum(gzs[g]) / len(gzs[g]) ** 0.5
            logps = np.log([2 * ss.norm.cdf(-abs(z)) for z in gzs[g]])
            gene_chi2ps[g] = 1 - ss.chi2.cdf(-2 * logps.sum(), 2 * len(gzs[g]))
    df = {'combined z':gene_zs}
    if 'Stouffer' in method:
        dn = data_name
        if 'Fisher' in data_name: dn += ' Stouffer'
        SLFDR = locfdr(zs=gene_zs, p0s=p0s, nbins=nbins, df_fit=df_fit, zthreshold=1.5/power_of_pval, data_name=data_name)
        lfdrs = np.array(list(SLFDR.values()))
        print("Combine genes by Stouffer's Z-score method (opposite changes cancel out):")
        if sum(lfdrs<LFDR0): print('gene    LFDR    shRNA : z-score')
        else: print('No genes are found below LFDR of ', LFDR0)
        items = [item for item in SLFDR.items() if abs(item[1])>=0] # remove NaN
        LFDR = sorted(items, key=operator.itemgetter(1))
        for g, lfdr in LFDR:
            if lfdr < LFDR0: 
                print(g, ' '*max(0, 8-len(g)), round(SLFDR[g], 3), end='  ', sep='')
                for t in gene_targets[g]: print(t, ':', round(zs[t],2), ' '*max(0, 11-len(t)), end='')
                print('')
        df['Stouffer LFDR'] = SLFDR
        print('')
    elif 'Fisher' in method:
        dn = data_name
        if 'Stouffer' in data_name: dn += ' Fisher'
        FLFDR = locfdr(ps=gene_chi2ps, p0s=p0s, nbins=nbins, df_fit=df_fit, zthreshold=1.5/power_of_pval, data_name=data_name)
        lfdrs = np.array(list(FLFDR.values()))
        print("Combine genes by p-values using Fisher's method (signs of change do not matter):")
        if sum(lfdrs<LFDR0): print('gene    LFDR    shRNA : z-score')
        else: print('No genes are found below LFDR of ', LFDR0)
        items = [item for item in FLFDR.items() if abs(item[1])>=0] # remove NaN
        LFDR = sorted(items, key=operator.itemgetter(1))
        for g, lfdr in LFDR:
            if lfdr < LFDR0: 
                print(g, ' '*max(0, 8-len(g)), round(FLFDR[g], 3), end='  ', sep='')
                for t in gene_targets[g]: print(t, ':', round(zs[t],2), ' '*max(0, 11-len(t)), end='')
                print('')
        df['Fisher LFDR'] = FLFDR
        print('')
    return pd.DataFrame(df)
