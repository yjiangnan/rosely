'''
Created on Jul 3, 2017

@author: jiangnan
'''
import os, numpy as np
import scipy.stats as ss
from scipy.special import beta, betainc
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
           'prepare_ttest_data', 'Result']

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
                if spl[-1] not in self.data[group]: 
                    self.data[group][spl[-1]] = {}
                for seqid in self.seqids:
                    self.data[group][spl[-1]][seqid] = csv[sample][seqid]
                self.samples[sample] = self.data[group][spl[-1]]
        
    def normalize_data(self, normalize_by='median', get_shared_presences=False, with_threshold=True):
        '''Normalize different samples to the same level and calculate important stats such as presences
        normalize_by: either 'median' (of all non-zero count seqids)
                      or 'quantileN' (1 <= N <= 9) to normalize by the mean of the 1/N proportion of data around the median
                      or a seqid (e.g. Luciferase with massive quantity) in the data
                      or a list of seqids whose median will be used for normalization. 
        '''
        for seqid in self.seqids: 
            self.total_count[seqid] = sum([self.data[group][bc][seqid] for group in self.data 
                                           for bc in self.data[group]])
        for seqid in self.seqids: 
            self.threshold[seqid] = (1 + self.cross_contamination_rate * self.total_count[seqid] * 3) * with_threshold
        if normalize_by in self.seqids: self.seqids.remove(normalize_by)
        self.presences = {}; self.shared_presences = {}
        self.sample_presences = {}; self.sample_weights = {}; self.sample_sizes = {}
        for group in self.data:
            self.presences[group] = {}
            self.sample_presences[group] = {}; self.sample_weights[group] = {}; self.sample_sizes[group] = {}
            for bc in self.data[group]: 
                self.sample_presences[group][bc] = sum([self.data[group][bc][seqid]>self.threshold[seqid] 
                                                     for seqid in self.data[group][bc]])
            median_presences = np.nanmedian(list(self.sample_presences[group].values()))
            std_presences = np.std(list(self.sample_presences[group].values()))
            for bc in self.data[group]: 
                p = self.sample_presences[group][bc] / (median_presences - std_presences)
                self.sample_weights[group][bc] = p * (p<1) + (p>=1)

            self.sample_sizes[group] = sum(list(self.sample_weights[group].values()))
            sw = self.sample_weights[group]
            for seqid in self.seqids: 
                self.presences[group][seqid] = sum([(self.data[group][bc][seqid] > self.threshold[seqid]) * sw[bc]
                                                   for bc in self.data[group]])
        self.normdata = {}; self.meanctrl = {}
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
            self.shared_presences
        self.nRefs = nRefs = {}
        nby = normalize_by
        if nby == 'median' or 'quantile' in nby: nby = list(self.seqids)
        for group in self.data: 
            for bc in self.data[group]:
                if type(nby) is list:
                    values = np.array([self.data[group][bc][seqid] for seqid in nby])
                    values = values[values>0]
                    if len(values): 
                        if 'quantile' in normalize_by: 
                            qt = int(normalize_by[-1]) * 2; l = len(values)
                            nRefs[group, bc] = np.mean(values[l//2 - l//qt : l//2 + l//qt])
                        else: nRefs[group, bc] = np.median(values)
                    else: nRefs[group, bc] = np.NaN
                else:
                    if normalize_by not in self.data[group][bc]: nRefs[group, bc] = np.NaN
                    else: nRefs[group, bc] = max(1, self.data[group][bc][normalize_by])
            self.mean_reference = np.nanmean(list(nRefs.values()))
        for group in self.data:
            self.normdata[group] = {}
            for bc in self.data[group]:
                self.normdata[group][bc] = {}
                nRef = nRefs[group, bc]
                for seqid in self.seqids: 
                    cnt = 0
                    if seqid in self.data[group][bc] and self.data[group][bc][seqid]>self.threshold[seqid]: 
                        cnt = self.data[group][bc][seqid]
                    self.normdata[group][bc][seqid] = cnt / nRef * self.mean_reference
            self.normdata[group] = pd.DataFrame(self.normdata[group])
        
    def get_data(self, bc, seqid):
        if seqid in self.samples[bc]: return self.samples[bc][seqid]
        else: return 0
        
    def get_threshold(self, seqid):
        if seqid in self.threshold: return self.threshold[seqid]
        else: return 0
    
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

def freqcdf(N1, X, N2, Y, r):
    '''Calculate the p-value (2-tailed) for the difference between frequencies X and Y in total trails of N1 and N2, 
    respectively, given a dilution ratio r (expected concentration of Y divided by X), using Bayesian inference.
    r is not simply the expected value of Y divided by that of X, but should be 
                log(1 - E(Y)/N2) / log(1 - E(X)/N1)
        i.e., when E(Y) is close to N1, the concentration of Y should be close to infinity.
        Consider infecting N cells with y viruses with a chance of p, i.e., yp viruses would enter cells. 
        Obviously, yp could increase to infinity but the total number Y of cells being infected can at 
        most be N. Then, we say the concentration of Y is yp for simplicity: yp = -log(1 - E(Y)/N)
    '''
    p1 = 1. / beta(X+1, N1-X+1)
    p2 = 1. / beta(Y+1, N2-Y+1)
    Et = 1 - beta(X+1, N1-X+1+r) / beta(X+1, N1-X+1)
    s = np.sign(Y/N2 - Et)
    def integrand(x, N1, X, N2, Y, r, s):
        # ss.binom.cdf(X, N, p) == betainc(N-X, X+1, 1-p), but betainc accepts non-integer parameters.
        # However, betainc(0, ., .) and betainc(., 0, .) are nan because gamma(0) == 0
        if s <= 0:  # Decreased, accumulate from 0 to Y 
            b = betainc(max(1e-19, N2 - Y),  Y+1,  (1-x)**r)
        else:       # Increased, accumulate from Y to N2, or 1 - cdf(from 0 to Y-1)
            b = 1 - betainc(N2 - (Y-1),  max(1e-19, Y),  (1-x)**r)
        return b * x ** X * (1 - x) ** (N1 - X)
    I1 = quad(integrand, 0, 1, args=(N1, X, N2, Y, r, s))[0]
    try: I2 = quad(integrand, 0, 1, args=(N2, Y, N1, X, 1/r, -s))[0]
    except: print('N2, Y, N1, X, r, s = ', N2, Y, N1, X, r, s, sep=', ')
    p = p1 * I1 + p2 * I2
    if p > 1: p = 2 - p
    pval, z = p, -ss.norm.ppf(p/2) * s
    return pval, z

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

def prepare_ttest_data(normdata_list, transform, minmean):
    seqids = []
    for data in normdata_list: 
        for sample in data: seqids += list(data[sample].keys())
    seqids = set(seqids)
    ttestdata = {}
    for seqid in seqids:
        ttestdata[seqid] = []; total = 0; n = 0
        for normdata in normdata_list:
            ttestdata[seqid].append([])
            for sample in sorted(normdata):
                d1 = normdata[sample][seqid]; total += d1; n += 1
                if transform == 'log':
                    if d1 >  0: ttestdata[seqid][-1].append(np.log2(d1))
                elif transform == 'log1':
                    if d1 > -1: ttestdata[seqid][-1].append(np.log2(1 + d1))
                elif type(transform) is float:
                    ttestdata[seqid][-1].append(d1 ** transform)
                else: raise Exception('transform has to be either log or a float number (power)')
        if total/n < minmean: ttestdata[seqid] = [[] for _ in range(len(normdata_list))]
    return ttestdata

def count_analysis(normdata_list, transform='log', minmean=np.nan, controls=None, 
                           paired=False, equalV=False,  debug=False, method='ascertained_ttest', pre_neutralize=True,
                           LFDRthr0=0.5, minr0=0, fine_tune=True,
                           data_name='data', nbins=30, df_fit=5, **kwargs):
    """Analyze the counts by ascertained_ttest
    
    Args:
    transform: how the data should be transformed for normalization,
               either 'log' for log-normal data (0s will be excluded),
               or     'log1 for log-normal data including 0s (excluding data <= -1),
               or     a power (0~1) for intermediate between log-normal and linear data, i.e., data**power should
                          be normally distributed. You can determine the power by function mean_median_norm(data_list)
    Other args are passed to ascertained_ttest and neup and locfdr.
    
    Returns:
    An object with the following members:
        ps:    The raw p-values
        zs:    The raw z-scores
        nps:   The neutralized p-values
        nzs:   The neutralized z-scores
        LFDR:  Local False Discovery Rate calculated from netralized p-values nps
        pop:   power of p-values
    """
    ttestdata = prepare_ttest_data(normdata_list, transform, minmean)
    if method=='ascertained_ttest':
        ps, zs, vbs = ascertained_ttest(ttestdata, [0, 1], controls, paired=paired, debug=True, 
                                        equalV=equalV, pre_neutralize=pre_neutralize)
    else: ps, zs, vbs = ttest(ttestdata, [0, 1], controls, paired=paired, equalV=equalV)
    res = Result(ps, zs, LFDRthr0, data_name=data_name, nbins=nbins, df_fit=df_fit, minr0=minr0, fine_tune=fine_tune, **kwargs)
    if debug: res.vbs = vbs
    res.results['dx'] = pd.Series(vbs['dx'])
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
