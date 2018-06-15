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
import pickle, operator
np.seterr(all='ignore')
import pandas as pd
from .ascertainedttest import ascertained_ttest, ttest
from .neutralstats import *
from collections import OrderedDict

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
    def __init__(self, filename, sep='\t', has_groups=True, cross_contamination_rate=default_cross_contamination_rate, **kwargs):
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
        self.total_count = {}; self.threshold ={}; self.samples = {}; 
        self.data = {}; self.norm_samples = {}; self.references = {}; self.extra = {}
        self.filename = filename
        if '.dat' == filename[-4:]: self.read_pickled()
        else: self.read_csv(sep, has_groups=has_groups, **kwargs)
        for group in self.data: 
            for bc in self.data[group]: # barcode or sample name
                self.samples[bc] = self.data[group][bc]
                for seqid in self.seqids:
                    if seqid not in self.data[group][bc]: self.data[group][bc][seqid] = 0
        
    def read_pickled(self):
        with open(os.path.expanduser(self.filename),'rb') as f: self.data, self.seqids = pickle.load(f)
        
    def read_csv(self, sep='\t', igroup=1, has_groups=True, **kwargs):
        if '.xls' in self.filename[-5:]: 
            csv = pd.read_excel(self.filename)
        else: csv = pd.read_csv(self.filename, sep=sep, **kwargs)
        self.seqids = list(csv.index); self.data = {}
        for sample in csv.keys():
            if has_groups: spl = sample.split('.')
            if has_groups and len(spl) <= 1:
                print('Warning:', sample, 'is not included (saved in self.extra) because it cannot be separated into group and sample.')
                self.extra[sample] = csv[sample]
            else:
                if has_groups: group = '.'.join(spl[:igroup])
                else: group = 'all'; spl = [sample]
                if group not in self.data: self.data[group] = {}
                self.samples[sample] = self.data[group][spl[-1]] = csv[sample]
        
    def normalize_data(self, normalize_by='median', normcutoff = 1, get_shared_presences=False, with_threshold=True):
        '''Normalize different samples to the same level and calculate important stats such as presences
        normalize_by: either 'median' (of all non-zero count seqids)
                      or 'quantileN' (1 <= N <= 9) to normalize by the mean of the 1/N proportion of data around the median
                      or a seqid (e.g. Luciferase with massive quantity) in the data
                      or a list of seqids whose median will be used for normalization.
         normcutoff: The minimum average value to be considered non-zero and included for calculating medium 
                     and quantile during normalization. 
        '''
        self.total_count = pd.DataFrame(self.samples).sum(axis=1)
        self.threshold = (1 + self.cross_contamination_rate * self.total_count * 3) * with_threshold
        if normalize_by in self.seqids: self.seqids.remove(normalize_by)
        self.presences = {}; self.shared_presences = {}
        self.sample_presences = {}; self.sample_weights = {}; self.sample_sizes = {}
        for group in self.data:
            df = pd.DataFrame(self.data[group])
            self.sample_presences[group] = {}; self.sample_weights[group] = {}
            for bc in self.data[group]: 
                sample = df[bc]
                sample[sample <= self.threshold] = 0
                self.sample_presences[group][bc] = (sample>0).sum()
            median_presences = np.nanmedian(list(self.sample_presences[group].values()))
            std_presences = np.std(list(self.sample_presences[group].values()))
            for bc in self.data[group]: 
                p = self.sample_presences[group][bc] / (median_presences - std_presences)
                self.sample_weights[group][bc] = p * (p<1) + (p>=1)

            self.sample_sizes[group] = sum(list(self.sample_weights[group].values()))
            sw = pd.Series(self.sample_weights[group])
            self.presences[group] = (df>0).multiply(sw).sum(axis=1)
        if get_shared_presences:
            for group in self.data:
                subgroup = ''
                if len(group.split('.'))>1: subgroup = group.split('.')[-1]
                self.shared_presences[subgroup] = {}
                for g in self.data:
                    sg = ''
                    if len(g.split('.'))>1: sg = g.split('.')
                    if g == group or subgroup != sg: continue
                    for seqid in self.seqids:
                        self.shared_presences[sg][seqid] = 0; th = self.threshold[seqid]
                        for spl in self.data[group]:
                            sw = min(self.sample_weights[group][spl], self.sample_weights[g][spl])
                            self.shared_presences[sg][seqid] += (self.data[group][spl][seqid] > th and
                                                              self.data[g    ][spl][seqid] > th) * sw
        self.nRefs = nRefs = {}
        nby = normalize_by
        if nby == 'median' or 'quantile' in nby: nby = list(self.seqids)
        allsamples = 0
        for spl in self.samples: allsamples = allsamples + self.samples[spl] / len(self.samples)
        allsamples = allsamples[allsamples > normcutoff]
        if 'quantile' in normalize_by: 
            qt = int(normalize_by[-1]) * 2; l = len(allsamples)
            mdall = sorted(allsamples)[l//2 - l//qt : l//2 + l//qt]
            allsamples = allsamples[(allsamples > mdall[0]) & (allsamples < mdall[-1])]
        for group in self.data: 
            for bc in self.data[group]:
                if type(nby) is list:
                    sample = self.data[group][bc]
                    if 'quantile' in normalize_by: 
                        nRefs[group, bc] = np.mean(sample.loc[allsamples.index])
                    else: nRefs[group, bc] = np.median(sample.loc[allsamples.index])
                    if nRefs[group, bc]==0: nRefs[group, bc] = np.NaN
                else:
                    if normalize_by not in self.data[group][bc]: nRefs[group, bc] = np.NaN
                    else: nRefs[group, bc] = max(1, self.data[group][bc][normalize_by])
            self.mean_reference = np.nanmean(list(nRefs.values()))
        self.normdata = OrderedDict()
        for group in self.data:
            self.normdata[group] = OrderedDict()
            for bc in self.data[group]:
                nRef = nRefs[group, bc]
                self.normdata[group][bc] = self.data[group][bc] / nRef * self.mean_reference
            self.normdata[group] = pd.DataFrame(self.normdata[group])
        
    def get_data(self, bc, seqid):
        if seqid in self.samples[bc]: return self.samples[bc][seqid]
        else: return 0
        
    def get_threshold(self, seqid):
        if seqid in self.threshold: return self.threshold[seqid]
        else: return 0
        
    def nRefs_as_weights(self, group, samples='all'):
        if samples is 'all': samples = self.normdata[group].keys()
        weights = {}
        ws = [self.nRefs[group, sample] for sample in samples]
        mw = np.nanmedian(ws)
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
                        
    def quality_heatmap(self, groups='all', samples='all', vmin=0.8, **kwargs):
        if groups == 'all': groups = sorted(self.data)
        if samples == 'all': group_samples = [(g, s) for g in groups for s in self.data[g]]
        elif len(groups) == 1: group_samples = list(zip(groups * len(samples), samples))
        else: group_samples = list(zip([[groups[i]] * len(samples[i]) for i in range(len(groups))],
                                  [s for spl in samples for s in spl]))
        M = pd.DataFrame(columns=pd.MultiIndex.from_tuples(group_samples), 
                         index=pd.MultiIndex.from_tuples(group_samples))
        for i, si in enumerate(group_samples):
            M[si][si] = 1
            for sj in group_samples[i+1:]:
                vs = pd.DataFrame({0:self.normdata[si[0]][si[1]], 1:self.normdata[sj[0]][sj[1]]}).dropna()
                r = np.corrcoef(vs[0], vs[1])[0][1]
                M[sj][si] = M[si][sj] = r
        import seaborn
        seaborn.heatmap(M.astype('float64'), vmin=vmin, **kwargs)
        mpl.pyplot.xlabel(''); mpl.pyplot.ylabel('')
        return M
    
        
class Result:
    def __init__(self, ps, zs, LFDRthr0, data_name='data', nbins=30, df_fit=5, minr0=0, fine_tune=True, **kwargs):
        self.ps = ps; self.zs = zs
        self.results = OrderedDict()
        self.neutralize_p_values(LFDRthr0, data_name=data_name, nbins=nbins, df_fit=df_fit, 
                                 minr0=minr0, fine_tune=fine_tune, **kwargs)
        
    def neutralize_p_values(self, LFDRthr0, data_name='data', nbins=30, df_fit=5, minr0=0, fine_tune=True, **kwargs):
        self.nps, self.LFDR, self.pop = neup(self.ps, LFDRthr0=LFDRthr0, data_name=data_name,
                                            nbins=nbins, df_fit=df_fit, minr0=minr0, fine_tune=fine_tune, **kwargs)
        self.nzs = {}
        try: 
            for i in self.nps.index: self.nzs[i] = np.sign(self.zs[i]) * abs(ss.norm.ppf(self.nps[i]/2))
        except: print(i, self.zs, self.nps); return
        self.results['LFDR'] = self.LFDR
        self.results['Controlled z-score'] = pd.Series(self.nzs)
        self.results['Controlled p-value'] = self.nps
        self.results = pd.DataFrame(self.results).sort_values(by='Controlled p-value')
        
    def plot_variance(self, sample_index='all', xname='mean'):
        if 'vbs' not in self.__dict__: 
            raise Exception('Variable vbs does not exist. Run count_analysis with debug=True to calculate it.')
        from matplotlib.pyplot import plot, xlabel, ylabel, subplot, legend, figure
        def plotx(xname, vbs):
            x = vbs[xname]
            logVs = np.log(vbs['Vs'])
            ElogVs = vbs['ElogVs']
            S2all = vbs['S2all']; Vsmpls = vbs['Vsmpls']; VlogV0s = vbs['VlogV0s']; Vadj = vbs['Vadj']
            plot(x, logVs, '.', label='all genes')
            plot(x, ElogVs, '.', label='mean log variance') 
            plot(x, ElogVs + S2all**0.5, '*r', label='SD(log variance)')
            plot(x, ElogVs - S2all**0.5, '*r')
            plot(x[:len(Vsmpls)], ElogVs + Vadj**0.5, '+c', label='Sampling SD')
            plot(x[:len(Vsmpls)], ElogVs - Vadj**0.5, '+c')
            plot(x, ElogVs + VlogV0s ** 0.5, '.m', label='prior gene SD')
            plot(x, ElogVs - VlogV0s ** 0.5, '.m') 
            xlabel(xname)
            ylabel('log variance')
            legend()
        if sample_index is 'all': sample_index = range(len(self.vbs)-2)
        elif type(sample_index) is int: sample_index = [sample_index]
        figure(figsize=(2 + 6*len(sample_index), 5))
        for i in sample_index:
            vbs = self.vbs[i]
            for v in vbs: vbs[v] = pd.Series(vbs[v])
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
    Et = 1 - np.exp(betaln(X+1, N1-X+1+r) - betaln(X+1, N1-X+1)) # exp(log beta) increases numeric stability 
    s = np.sign(Y/N2 - Et)
    def integrand(x, N1, X, N2, Y, r, s, lnb):
        # ss.binom.cdf(X, N, p) == betainc(N-X, X+1, 1-p), but betainc accepts non-integer parameters.
        # However, betainc(0, ., .) and betainc(., 0, .) are nan because gamma(0) == 0
        if s <= 0:  # Decreased, accumulate from 0 to Y 
            b = betainc(max(1e-19, N2 - Y),  Y+1,  (1-x)**r)
        else:       # Increased, accumulate from Y to N2, or 1 - cdf(from 0 to Y-1)
            b = 1 - betainc(N2 - (Y-1),  max(1e-19, Y),  (1-x)**r)
        return b * np.exp( X * np.log(x) + (N1 - X) * np.log(1 - x) - lnb)
    points = [(X/N1 + Y/N2)/2, Et]
    I1 = quad(integrand, 0, 1, args=(N1, X, N2, Y, r, s, lnb1), points=points)[0]
    try: I2 = quad(integrand, 0, 1, args=(N2, Y, N1, X, 1/r, -s, lnb2), points=points)[0]
    except: print('N2, Y, N1, X, r, s = ', N2, Y, N1, X, r, s, sep=', ')
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
                              LFDRthr0=0.5, minr0=0, data_name='data', nbins=30, df_fit=5, **kwargs):
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
    nc1 = sorted([presences1[seqid] - shared_presences[seqid] for seqid in ctrls])
    nc2 = sorted([presences2[seqid] - shared_presences[seqid] for seqid in ctrls])
    l1 = len(nc1); l2 = len(nc2)
    mean1 = np.mean(nc1[l1//4:l1//4*3]); mean2 = np.mean(nc2[l2//4:l2//4*3])
    mdn = np.mean(list(shared_presences.values()))
    r = np.log(1 - mean2/(N2-mdn)) / np.log(1 - mean1/(N1-mdn))
    print('N1: {:.2f}  N2: {:.2f}  25~75% mean1: {:.2f}  25~75% mean2: {:.2f}  mean shared: {:.2f}  MOI ratio: {:.2f}'.format(
        N1, N2, mean1, mean2, mdn, r))
    for seqid in seqids:
        dn = shared_presences[seqid]
        n1 = presences1[seqid]; n2 = presences2[seqid]; z = 0
        if n1 + n2 - dn*2 >= min_count and (N1-n1) + (N2-n2) >= min_count: 
            pval, z = freqcdf(N1-dn, n1-dn, N2-dn, n2-dn, r)
        else: z = np.NaN; pval = np.NaN
        zs[seqid] = z
        ps[seqid] = pval
    res = Result(ps, zs, LFDRthr0, data_name=data_name, nbins=nbins, df_fit=df_fit, minr0=minr0, **kwargs)
    change = {}
    for seqid in ps: 
        change[seqid] = '{:.3f} --> {:.3f}'.format(
                        presences1[seqid], presences2[seqid])
    res.results['presence change'] = pd.Series(change)
    return res

def prepare_ttest_data(normdata_list, transform, minmean, sample_weights):
    seqids = []
    weights = None
    if sample_weights is not None: weights = {}
    for data in normdata_list: 
        for sample in data: seqids += list(data[sample].keys())
    seqids = set(seqids)
    ttestdata = {}
    for seqid in seqids:
        ttestdata[seqid] = []; total = 0; n = 0
        if sample_weights is not None: weights[seqid] = []
        for i, normdata in enumerate(normdata_list):
            ttestdata[seqid].append([])
            if sample_weights is not None: weights[seqid].append([])
            for sample in sorted(normdata):
                d1 = normdata[sample][seqid]
                if transform == 'log':
                    if d1 >  0: ttestdata[seqid][-1].append(np.log2(d1))
                elif type(transform) is str and transform[:3] == 'log':
                    n0 = float(transform[3:])
                    if '.' not in transform[3:]: n0 = int(n0)
                    if str(n0) != transform[3:]: 
                        raise Exception('The string following log is not a well formated float number')
                    if d1 > -n0: ttestdata[seqid][-1].append(np.log2(n0 + d1))
                elif type(transform) is float:
                    ttestdata[seqid][-1].append(d1 ** transform)
                else: raise Exception('transform has to be either log or a float number (power)')
                w = 1
                if sample_weights is not None and len(weights[seqid][-1]) < len(ttestdata[seqid][-1]): 
                    w = sample_weights[i][sample]
                    weights[seqid][-1].append(w)
                total += d1 * w; n += w
        if total/n < minmean: 
            del ttestdata[seqid]
            if sample_weights is not None: del weights[seqid]
    return ttestdata, weights

def count_analysis(normdata_list, transform=1., minmean=np.nan, controls=None, sample_weights=None,
                           paired=False, equalV=False,  debug=False, method='ascertained_ttest', pre_neutralize=True, 
                           span=0.5, LFDRthr0=0.5, minr0=0, fine_tune=True,
                           data_name='data', nbins=30, df_fit=5, **kwargs):
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
    if np.isnan(minmean) and transform == 'log1': minmean = 10
    ttestdata, weights = prepare_ttest_data(normdata_list, transform, minmean, sample_weights)
    if method=='ascertained_ttest':
        ps, zs, vbs = ascertained_ttest(ttestdata, [0, 1], controls, paired=paired, debug=True, weights=weights,
                                        equalV=equalV, pre_neutralize=pre_neutralize, span=span)
        dxs = vbs['dx']
    else: ps, zs, dxs, vbs = ttest(ttestdata, [0, 1], controls, paired=paired, weights=weights, equalV=equalV)
    res = Result(ps, zs, LFDRthr0, data_name=data_name, nbins=nbins, df_fit=df_fit, minr0=minr0, fine_tune=fine_tune, **kwargs)
    if debug: res.vbs = vbs
    nm = 'dx'
    if 'log' in transform: nm = 'log2 fold change'
    res.results[nm] = pd.Series(dxs)
    return res

def combine_by_gene(zs, sep='.', LFDR0=0.3, p0s=[0, ], nbins=30, df_fit=5, power_of_pval=1, data_name='data'):
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
    for seqid in zs: genes.append(seqid.split(sep)[0])
    genes = list(set(genes))
    gene_zs = {}; gene_chi2ps = {}; gzs = {}; gene_targets = {}
    for seqid in zs:
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
    SLFDR = locfdr(zs=gene_zs, p0s=p0s, nbins=nbins, df_fit=df_fit, zthreshold=1.5/power_of_pval, data_name=data_name + ' Stouffer')
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
    print('')

    FLFDR = locfdr(ps=gene_chi2ps, p0s=p0s, nbins=nbins, df_fit=df_fit, zthreshold=1.5/power_of_pval, data_name=data_name + ' Fisher')
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
    print('')
    return pd.DataFrame({'combined z':gene_zs, 'Stouffer LFDR':SLFDR, 'Fisher LFDR':FLFDR})
