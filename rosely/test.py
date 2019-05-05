
import imp, rosely, time, operator, pickle
from pylab import *
from numpy import *
import scipy.stats as ss
import pandas as pd

from rosely import *
import GEOparse
import pandas as pd
import rosely
imp.reload(rosely)
from rosely import *
from collections import OrderedDict
results = {'Mouse':{}, 'Human':{}}

with open('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/group_bc_target.dat','rb') as fl: data, shRNAs = pickle.load(fl)
shRNAs.remove('Luciferase')


controls = []
for shRNA in shRNAs:
    if 'Plate.6' in shRNA: controls.append(shRNA)
balance = {}
HEKgroups = ['HEK002500', 'HEK008000', 'HEK030000', 'HEK060000', 'HEK1e5']
ncells = [2500, 8000, 30000, 60000, 1e5]
ncycles = [35-np.log2(ncell/2500) for ncell in ncells]
group = 'HEK200'; nLucs = []
with open('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/ms_41hm_luc_200BJ.group_bc_target.dat','rb') as f: bdata, bshRNAs = pickle.load(f)
for bc in bdata[group]: 
    if 'Luciferase' not in bdata[group][bc]: bdata[group][bc]['Luciferase'] = 1
    nLucs.append(bdata[group][bc]['Luciferase'])
meanLuc = mean(nLucs)
ESS = array(nLucs)/meanLuc; ESS[ESS>1] = 1; sumESS = sum(ESS)
def get_balance_for(shRNA):
    b = 0; 
    for bc in bdata[group]:
        nLuc = bdata[group][bc]['Luciferase']; 
        if shRNA in bdata[group][bc] and bdata[group][bc][shRNA] > 3: 
            b += min(1, nLuc/meanLuc)
    return b
for shRNA in shRNAs:
    balance[shRNA] = get_balance_for(shRNA)

print('sumESS:', sumESS)
for shRNA in shRNAs: balance[shRNA] = -log(1 - balance[shRNA]/sumESS)
for shRNA in shRNAs: balance[shRNA] *= len(balance) / sum(list(balance.values()))
meanctrl = np.mean([balance[ctrl] for ctrl in controls])
for shRNA in shRNAs: balance[shRNA] /= meanctrl
len(controls), controls[:3], len(balance), meanctrl
short = ['Aco2.1603_Plate.1_A09', 'Aco2.571_Plate.1_A07', 'Aco2.928_Plate.1_A08', 'Actn3.1064_Plate.1_A14', 'Actn3.2563_Plate.1_A13', 'Actn3.551_Plate.1_A15', 'Arcn1.2201_Plate.1_B21', 'Arcn1.2325_Plate.1_B20', 'Arcn1.3214_Plate.1_B19', 'Atp2a1.1242_Plate.1_C17', 'Atp2a1.1620_Plate.1_C16', 'Atp2a1.1674_Plate.1_C18', 'Capn3.1149_Plate.1_D14', 'Capn3.1204_Plate.1_D15', 'Capn3.1794_Plate.1_D13', 'Cryba2.252_Plate.1_G18', 'Cryba2.676_Plate.1_G17', 'Cryba2.677_Plate.1_G16', 'Fam166a.31_Plate.1_K09', 'Fam166a.36_Plate.1_K08', 'Fam166a.761_Plate.1_K07', 'Gja1.1610_Plate.1_L21', 'Gja1.1980_Plate.1_L20', 'Gja1.6_Plate.1_L19', 'Klc2.1047_Plate.1_P18', 'Klc2.1762_Plate.1_P17', 'Klc2.2977_Plate.1_P16', 'Mybpc1.1347_Plate.2_D19', 'Mybpc1.1348_Plate.2_D21', 'Mybpc1.3244_Plate.2_D20', 'Myh13.2207_Plate.2_D24', 'Myh13.3324_Plate.2_D22', 'Myh13.3488_Plate.2_D23', 'Myh4.2708_Plate.2_E01', 'Myh4.3324_Plate.2_E02', 'Myh4.638_Plate.2_E03', 'Pitpnm1.105_Plate.2_G24', 'Pitpnm1.2892_Plate.2_G23', 'Pitpnm1.42_Plate.2_G22', 'Psmb7.198_Plate.2_I02', 'Psmb7.351_Plate.2_I01', 'Psmb7.352_Plate.2_I03', 'Slc12a4.270_Plate.2_L09', 'Slc12a4.3090_Plate.2_L08', 'Slc12a4.3312_Plate.2_L07', 'Slc16a3.1161_Plate.2_L12', 'Slc16a3.347_Plate.2_L10', 'Slc16a3.770_Plate.2_L11', 'Tecta.3794_Plate.2_N12', 'Tecta.4895_Plate.2_N10', 'Tecta.5399_Plate.2_N11', 'Tusc5.2581_Plate.3_B11', 'Tusc5.2582_Plate.3_B10', 'Tusc5.2586_Plate.3_B12', 'Xkr6.1309_Plate.3_E06', 'Xkr6.1454_Plate.3_E05', 'Xkr6.2101_Plate.3_E04']
long = [shRNA for shRNA in shRNAs if shRNA not in controls and shRNA not in short]

genes = []; geneshort = []; genelong = []; genectrl = []
for shRNA in shRNAs:
    g = shRNA.split('.')[0]
    if g not in genes: genes.append(g)
    if shRNA in short and g not in geneshort: geneshort.append(g)
    if shRNA in long  and g not in genelong:  genelong.append(g)
    if shRNA in controls and g not in genectrl: genectrl.append(g)

failed_samples = ['GGTAATGG', 'AGAAGCTG', 'ATATAATT', 'AACCGAGC', 'AGGTTAGG', 'ATAACACA', 'AGAGGTAT']
n_failed_Luciferase = [281,         71,        102,        161,      100,        1762,      122]
#cell type:            1stCMP       1stHSC     1stMPP   2ndB220low   2ndCMP,     2ndHSC     2ndMPP
# 1761 is likely cross contamination; others are barcode mutation

crh = CountReads('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/group_bc_target.dat')
cr200 = CountReads('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/ms_41hm_luc_200BJ.group_bc_target.dat')
cr120 = CountReads('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/ms_41hm_luc_120BJ_NoIndex_L008_R1_001.group_bc_target.dat')
cr500 = CountReads('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/hek500_group_bc_target.dat')
cr1 = CountReads('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/1st_group_bc_target.dat')
cr2 = CountReads('/Users/jiangnan/Documents/rudolphLab/mouseTransplant/2nd_group_bc_target.dat')
for cr in [crh, cr200, cr120, cr500, cr1, cr2]: 
    for shRNA in shRNAs:
        if shRNA not in cr.seqids: 
            cr.seqids.append(shRNA)
            cr.data.loc[shRNA] = 0
    cr.normalize_data(normalize_by='Luciferase', normcutoff=2, with_threshold=True)
group = 'HSC'
mmx = [min(cr1.presences[group].values)-1, max(cr1.presences[group].values)+1]
change_of_shRNAs = {}
for group in cr1.data.keys().levels[0]:
    change_of_shRNAs[group] = {}
    for transplant in ['1st', '2nd']:
        change_of_shRNAs[group][transplant] = {}
        for shRNA in cr1.seqids: 
            change_of_shRNAs[group][transplant][shRNA] = {}
            
pvalues = {}; zvalues = {}; ATresults = {}
groupmap = {'HSC':'HEK002500', 'MPP':'HEK008000', 'B220low':'HEK060000', 'CMP':'HEK1e5'}
for group in ['HSC']:
    figure
    res1 = count_analysis([crh.normdata[groupmap[group]], cr1.normdata[group]], transform='log', pre_neutralize=False,
                           controls=controls, LFDRthr0=0.5, minr0=0., data_name=group+' 1st', 
                           parallel=False, weighted=False, debug=True)
    subplot(121)
    text(-4, ylim()[1]*0.93, 'a', fontsize=18)
    subplot(122)
    text(0., ylim()[1]*0.93, 'b', fontsize=18)
    figure
    res2 = count_analysis([cr1.normdata[group], cr2.normdata[group]], transform='log', pre_neutralize=False,
                           controls=controls, LFDRthr0=0.5, minr0=0.00, 
                           parallel=False, data_name=group+' 2nd', weighted=False, debug=True)
    subplot(121)
    text(-4, ylim()[1]*0.93, 'c', fontsize=18)
    subplot(122)
    text(0., ylim()[1]*0.93, 'd', fontsize=18)
    res3 = count_analysis([crh.normdata[groupmap[group]], cr2.normdata[group]], transform='log', pre_neutralize=False,
                           controls=controls, LFDRthr0=0.5, minr0=0.00, nbins=50, 
                          parallel=False, data_name=group+' 0th 2nd', weighted=False, debug=True)

    pvalues[group], zvalues[group] = [res1.nps, res2.nps, res3.nps], [res1.nzs, res2.nzs, res3.nzs]
    ATresults[group] = [res1, res2, res3]
    
        

pvalues = {}; zvalues = {}; ATresults = {}
groupmap = {'HSC':'HEK002500', 'MPP':'HEK008000', 'B220low':'HEK060000', 'CMP':'HEK1e5'}
for group in ['HSC']:
    figure
    res1 = count_analysis([crh.normdata[groupmap[group]], cr1.normdata[group]], transform='log',
                           controls=controls, LFDRthr0=0.5, minr0=0., data_name=group+' 1st')
    subplot(121)
    text(-4, ylim()[1]*0.93, 'a', fontsize=18)
    subplot(122)
    text(0., ylim()[1]*0.93, 'b', fontsize=18)
    figure
    res2 = count_analysis([cr1.normdata[group], cr2.normdata[group]], transform='log',
                           controls=controls, LFDRthr0=0.5, minr0=0.00, data_name=group+' 2nd')
    subplot(121)
    text(-4, ylim()[1]*0.93, 'c', fontsize=18)
    subplot(122)
    text(0., ylim()[1]*0.93, 'd', fontsize=18)
    res3 = count_analysis([crh.normdata[groupmap[group]], cr2.normdata[group]], transform='log',
                           controls=controls, LFDRthr0=0.5, minr0=0.00, nbins=50, data_name=group+' 0th 2st')

    pvalues[group], zvalues[group] = [res1.nps, res2.nps, res3.nps], [res1.nzs, res2.nzs, res3.nzs]
    ATresults[group] = [res1, res2, res3]

N = 2000; ntop = int(N*0.05); pops = {}

# def create_simu_data(mdif, nrep, ntechrep=1, prop=0.05, gradual=False):
#     simudata = pd.DataFrame({(i,j):[] for i in range(2) for j in range(nrep)})
#     for i in range(2*nrep*ntechrep): 
#         group = i//(nrep*ntechrep); sample = i % (nrep*ntechrep)
#         if sample < nrep:
#             r = random.normal(0, 1, size=(N,))
#             if group==0: 
#                 if gradual: r[:int(N*prop)] += 2 * mdif * (1 + array(range(int(N*prop)))) / int(N*prop)
#                 else: r[:int(N*prop)] += mdif
#         else:
#             r = simudata[(group, sample % nrep)]
#         simudata[(group, sample)] = r
#     return simudata
# 
# 
# N = 10000
# data = create_simu_data(mdif=2, nrep=5, ntechrep=1, prop=0.8)
# _, ps8  = ss.ttest_ind(data[0].values, data[1].values, equal_var=True, axis=1)
# 
# data = create_simu_data(mdif=2, nrep=5, ntechrep=1, prop=0.8, gradual=True)
# _, ps8g = ss.ttest_ind(data[0].values, data[1].values, equal_var=True, axis=1)
# 
# figure(figsize=(12, 10))
# subplot(221)
# LFDRthr0 = 0.5
# neup(ps8, LFDRthr0=LFDRthr0, minr0=0, subplot_label='a', plotLFDR=False);
# text(0.2, 0.9, '80% different by 2σ')
# 
# subplot(222)
# neup(ps8g, LFDRthr0=LFDRthr0, minr0=0, subplot_label='b', plotLFDR=False);
# text(0.2, 0.9, '80% different by 2σ on average')

# ps = random.uniform(size=(10000,))
# ps[:100] **= 3
# ps **= 3
# ps[0] = 1e-299; ps[1] = 1e-199
# neup(ps, minr0=0.00)
# 
# wt = CountReads('~/Downloads/WT_raw.tsv',   sep='\t', has_groups=False, header=None)
# 
# for i in [21, 22, 25, 28, 34, 36]: del wt.data['all'][i]
# 
# wt.normalize_data(normalize_by='quantile5', normcutoff=10)
# wt_keys = sorted(wt.data['all'])
# 
# n = 3
# res = count_analysis([wt.subgroup('all', [44, 31, 13, 42, 46]),
#                       wt.subgroup('all', [ 4, 24, 35, 41, 23]),], span=0.5, method='',
#                         transform='log2.5', minmean=5, debug=True, with_plot=True, do_SVA=False,
#                         pre_neutralize=True, fine_tune=True, weighted=True)

def load_data(fn, normcutoff, **kwargs):
    cr = CountReads('~/Downloads/' + fn, **kwargs)
    cr.normalize_data(normalize_by='quantile5', normcutoff=normcutoff)
    print('Library size:', len(cr.seqids))
    groups = list(cr.data.keys())
    print('groups:', groups)
#     for n in cr.nRefs: print(n, ':', cr.nRefs[n])
    return cr

def analyze(sp, cr, ga, group1, group2, ref, grp='', **kwargs2):
    pops = {}; sig = {}
    for (do_SVA, ridge, name, method) in [(False, False, 'Ordinary', 'ascertained_ttest'), 
                                          (True, False, 'SVA', ''), 
                                          (True, True, 'SVA ridge', ''),
                                          (False, False, 'Ascertained', 'ascertained_ttest')
                                 ]:
        datalist = [cr.normdata[group1], cr.normdata[group2]]
        res = count_analysis(datalist, 
                        span=0.8, method=method, do_SVA=do_SVA, ridge = ridge, nSV=2, debug=True, with_plot=False,
                        LFDRthr0=0.5, minn=2, df_fit=4, minr0=0., data_name=name, parallel=False, **kwargs2)
        pops[name + ' power of p values'] = res.pop
        sig [name + ' LFDR < 0.3'] = sum(res.LFDR < 0.3)
        if do_SVA: pops[name + ' GLM power of p values'] = res.vbs['SVApop']
    n1 = cr.normdata[group1].shape[1]
    n2 = cr.normdata[group2].shape[1]
    compare = grp + group1 + '(' + str(n1) + ') ~ ' + group2 +'(' + str(n2) + ')'
    results[sp][(ga, compare)] = OrderedDict({**pops, **sig,
                              'Mean of n number': (n1 + n2) / 2,
                              '# of analyzed genes':res.results.shape[0],
                              'Citation': ref
                  })
    print(pops)
    return res


ga = 'GSE97976'
cr = load_data('GSE97976_featureCounts.mm10.e81.csv', 10, sep='\t', igroup=1)
ref = 'Schmidt K, Zhang Q, Tasdogan A, Petzold A et al. The H3K4 methyltransferase Setd1b is essential for hematopoietic stem and progenitor cell homeostasis in mice. Elife 2018 Jun 19;7. PMID: 29916805'
res = analyze('Mouse', cr, ga, 'HOM', 'HET', ref, transform='log5', minmean=10)
ga = 'GSE79442'
cr = load_data('GSE79442_FilteredTargetCounts.csv.gz', 10, sep=',', 
                igroup=1, group_sample_sep='-')
ref = 'Edgar JM, Robinson M, Willerth SM. Fibrin hydrogels induce mixed dorsal/ventral spinal neuron identities during differentiation of human induced pluripotent stem cells. Acta Biomater 2017 Mar 15;51:237-245. PMID: 28088670'
res = analyze('Human', cr, ga, 'D0', 'D5', ref, transform='log5', minmean=10)


cr = CountReads('~/Downloads/SrHt_expressionMatrix.txt', sep='\t', igroup=1, 
                group_sample_sep='ical', has_groups=True)
cr.normalize_data(normalize_by='quantile5', normcutoff=1)
print('Library size:', len(cr.seqids))
groups = list(cr.presences.keys())
print('groups:', groups)

del cr.normdata['class']['15']
del cr.normdata['nonclass']['06']

monocyte = count_analysis([cr.normdata['class'], 
                      cr.normdata['nonclass']], 
#                      sample_weights=[cr.nRefs_as_weights('class'), cr.nRefs_as_weights('nonclass')],
                     paired=True, pre_neutralize=True, span=0.6, parallel=False,
                    LFDRthr0=0.5, transform='log5', minmean=10, debug=True, minr0=0.0)
text(-0.03, 0.96, 'f', fontsize=24)
text(0.2, 0.85, 'monocyte study', fontsize=18)
subplot(121)
gca().set_visible(False)

wt = CountReads('~/Downloads/WT_raw.tsv',   sep='\t', has_groups=False, header=None)
ko = CountReads('~/Downloads/Snf2_raw.tsv', sep='\t', has_groups=False, header=None)

for i in [21, 22, 25, 28, 34, 36]: del wt.data['all'][i]
for i in [6, 13, 25, 35]: del ko.data['all'][i]

wt.normalize_data(normalize_by='quantile5', normcutoff=1)
ko.normalize_data(normalize_by='quantile5', normcutoff=1)
wt_keys = sorted(wt.data['all'])
ko_keys = sorted(ko.data['all'])

postclean = count_analysis([wt.subgroup('all', [27, 47, 45]),
                            wt.subgroup('all', [17, 26, 43])],
                            transform=0.25, minmean=10, debug=True, with_plot=True, span=0.8, #method='', 
                            LFDRthr0 = 0.5, pre_neutralize=True, fine_tune=False);

cr = CountReads('~/Downloads/GSE97976_featureCounts.mm10.e81.csv', sep='\t', igroup=1, 
               )
cr.normalize_data(normalize_by='quantile5', normcutoff=1)
print('Library size:', len(cr.seqids))
groups = list(cr.data.keys())
print('groups:', groups)
cr.nRefs

res = count_analysis([
                cr.normdata['HOM'], 
                cr.normdata['HET'],
                ], 
                paired=False, pre_neutralize=False, span=0.6, method='',
                do_SVA = False, debug=True, ridge=True,
                LFDRthr0=0.5, transform='log5', minmean=10, minn=2, df_fit=4, minr0=0., nSV=2, penaltyRatioToSV=1)

time.sleep(0.1)
