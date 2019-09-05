'''
This module implements ascertained t-test (Welch test) that uses the similarity of genes/constructs in the same 
experimental group to reduce the uncertainty of the variance of each gene/construct in an empirical Bayesian approach. 
The reduced uncertainty is converted back to a number for degrees of freedom to be compatible with conventional statistics. 

Created on Jun 13, 2017

@author: Jiang-Nan Yang
'''

import scipy.stats as ss, pandas as pd
from .lowess import loess_fit
from pylab import *
import numpy as np, warnings
from scipy.optimize import minimize, fsolve
from scipy.special import gammaln, polygamma, gammainc
from scipy.integrate import quad
import scipy.linalg
from multiprocessing.pool import Pool
import psutil
from _functools import partial
warnings.filterwarnings("error", category=RuntimeWarning)
from .neutralstats import neup
from . import integral as itg

__all__ = ['ascertained_ttest', 'mean_median_norm', 'segmented_mean_median_norm', 'ttest', 'ttest2', 'itg']

def partial_ridge_sva_residual(data, controls, nSV, penaltyRatioToSV=0., ridge=True, parallel=True):
    """
    Surrogate variable analysis:
    An improved and explicit surrogate variable analysis procedure by coefficient adjustment
        Seunggeun Lee Wei Sun Fred A. Wright Fei Zou
    Biometrika, Volume 104, Issue 2, 1 June 2017, Pages 303–316, https://doi.org/10.1093/biomet/asx018
    
    Ridge regression:
    Boonstra PS, Mukherjee B, and Taylor JMG (2015). 
    "A small-sample choice of the tuning parameter in ridge regression." Statistica Sinica 25, 1185-1206
    
    For partial ridge regression (ridge=True), we could apply a penalty only to the surrogate variables by  
    setting penaltyRatioToSV=0 to avoid their unnecessary removal of variations resulted from the design 
    variables and avoid problems of multicollinearity. 
    
    nSV: number of Surrogate Variables to use
    """
    dt = data.loc[controls].dropna().transpose()
#     i0 = 100
#     for i in range(i0, i0+200): 
#         dt.loc[0, 0][i] += (i-i0+50)/50.; dt.loc[0, 1][i] -= (i-i0+50)/50.; 
#         dt.loc[1, 1][i] -= (i-i0+50)/50.
    X0 = np.array([(k[0], 1.) for k in list(dt.index)])
    X = X0 + 0.
    for i in range(X.shape[1] - 1):
        X[:, i] -= np.mean(X[:, i])
        X[:, i] /= np.std (X[:, i])
    inv = np.linalg.inv
    XXX = inv(X.T @ X) @ X.T
    Bh = XXX @ dt
    M = X @ XXX
    R = (np.eye(M.shape[0]) - M) @ dt
    
    U, _, _ = scipy.linalg.svd(R)
    Uq = U[:, :nSV]
    F = inv(Uq.T @ Uq) @ Uq.T @ R
    J = np.ones((F.shape[1], 1))
    IMJFt = (np.eye(F.shape[1]) - J @ inv(J.T @ J) @ J.T) @ F.T
    S = Uq + X @ Bh @ IMJFt @ inv(F @ IMJFt) # S is the surrogate variables
#     print('Surrogate Variables:')
#     print(S)
    
    for i in range(S.shape[1]): # Standardize surrogate variables.
        S[:, i] -= np.mean(S[:, i])
        S[:, i] /= np.std (S[:, i])
    XS1 = np.concatenate((X0, S), axis=1)
    B = []
    ses = []; #Fs = []
    values = data.get_values()
    ys = (values - np.nanmean(values, axis=1)[:, np.newaxis]).tolist()
    nparam =  XS1.shape[1]
    dfes = (data.notnull().sum(axis=1).values - nparam).tolist()
    XS0 = np.concatenate((X[:, :-1], S), axis=1)
    XSXS0 = XS0.T @ XS0; I0 = np.eye(XS0.shape[0])
    XS1XS1 = XS1.T @ XS1
    
    ptrg = partial(partial_ridge, X=X, S=S, XS0=XS0, XSXS0=XSXS0, XS1=XS1, XS1XS1=XS1XS1, I0=I0, 
                        ridge=ridge, nSV=nSV, nparam=nparam, penaltyRatioToSV=penaltyRatioToSV)
    if parallel:
        pool = Pool(psutil.cpu_count(logical=False))
        results = pool.map(ptrg, zip(ys, values, dfes))
        pool.close()
    else: results = [ptrg(p) for p in zip(ys, values, dfes)]
    for (beta, se) in results:
        B.append(beta); ses.append(se)
    B = np.transpose(B)
    pvals = 2*(ss.t.cdf(-np.abs(B[0]/ses), dfes))
    _, _, pop = neup(pvals, with_plot=False, minr0=0, fine_tune=False)
    print('SVA pop of design variable'+' for ridge regression'*ridge+':', pop)
    sum_inv_counts = sum([1. / X0[:,0].tolist().count(x) for x in set(X0[:,0])])
    scale_var = (data.T - XS1 @ B).var(axis=0, ddof=2) * sum_inv_counts / (np.array(ses)**2)
    res = data.T - S @ B[X.shape[1]:, :]
    return res.T, pop, scale_var

def partial_ridge(params, X, S, XS0, XSXS0, XS1, XS1XS1, I0, ridge, nSV, nparam, 
                  penaltyRatioToSV):
    y, y0, dfe = params
    valid = np.logical_not(np.isnan(y))
    inv = np.linalg.inv
    if ridge:
        def ridge_GCVc(lmd): # return errors of the generalized cross validation.
            lmd = lmd[0]
            L = np.diag([lmd * penaltyRatioToSV] * (X.shape[1]-1) + [lmd] * nSV)
            P = XS @ inv(XSXS + L) @ XS.T
            IP2 = (I - P) @ (I - P)
            e = np.log(y @ IP2 @ y) - 2 * np.log(max(1e-99, 1 - np.trace(P)/n - 2./n))
            return e
        n = valid.sum()
        if n > nparam:
            y = np.array(y)[valid]
            if n == len(y0): 
                XS = XS0; XSXS = XSXS0; I = I0
            else: 
                XS = np.concatenate((X[valid, :-1], S[valid, :]), axis=1)
                XSXS = XS.T @ XS
                I = np.eye(n)
            lmda = minimize(ridge_GCVc, 1, bounds=[(0.0001, 200)]).x[0]
        else:
            lmda = np.nan
    else: lmda = 0
    L  = np.diag([0] * X.shape[1] + [lmda] * nSV)
    Lr = np.diag([lmda * penaltyRatioToSV / 4 * 0] * X.shape[1] + [lmda] * nSV)
    # /4 is used because std of the original X0 is only 0.5, beta would then be doubled
    if all(valid):
        D  = inv(XS1XS1 + L ) @ XS1.T
        Dr = inv(XS1XS1 + Lr) @ XS1.T
    else:
        XS1v = XS1[valid, :]
        D  = inv(XS1v.T @ XS1v + L ) @ XS1v.T
        Dr = inv(XS1v.T @ XS1v + Lr) @ XS1v.T
    beta  = D  @ y0[valid]
    betar = Dr @ y0[valid]
    yhat = XS1 @ beta
    yhat[~valid] = np.nan
    sse = np.nansum((y0-yhat)**2); #   % sum of squared errors
    mse = sse/dfe;
    if dfe>0: se = (D[0] @ D[0] * mse) ** 0.5  # http://web.as.uky.edu/statistics/users/pbreheny/764-F11/notes/9-1.pdf
    else: se = np.nan
    beta[1:] = betar[1:]
    return beta, se

def cross_group_normalization(data, idxes, controls):
    allms = []
    for i in idxes:       
        ms = data[i].loc[controls].mean(axis=1)
        allms.append(ms.mean())
    mall = mean(allms)
    for i in idxes:
        data[i] += mall - allms[i]
    return data

def normalize_by_least_change(data, equalV, min_pval, weights, paired):
    idxes = data.columns.levels[0].tolist()
    means = data.mean(axis=1).values
    vs = means.copy()
    vs.sort()
    included = (means > vs[len(vs)//3]) & (means < vs[len(vs)//4*3])
    mvns = calc_mvns(data, weights, paired)
    _, pvals, _ = calc_ttest_stats(mvns, paired, equalV, need_z=False)
    for idx in idxes[:-1]:
        ps = pvals[idx]
        nps, _, _ = neup(np.array(ps), with_plot=False, minr0=0, fine_tune=False)
        included = np.logical_and(nps > min_pval, included)
    for idx in idxes: # means of groups are calculated separately in case of sample size difference 
        in_sample_means = data[idx][included].mean()
        null = in_sample_means.isnull()
        in_sample_means[null] = data[idx][in_sample_means[null].index].mean()
        data[idx] += in_sample_means.mean() - in_sample_means # cross-sample normalization
    nd = cross_group_normalization(data, idxes, data[included].index)
    return nd, data[included].index

def calc_mvns(data, weights, paired, scale_var=1, repeat=1):
    idxes = data.columns.levels[0].tolist()
    mvns = {}
    for i in range(len(idxes) - paired):
        idx = idxes[i]
        if paired:
            md = (data[idx+1] + data[idx]) * 0.5
            dd = data[idx+1] - data[idx]
        else: 
            md = dd = data[idx]
        mvn = {}
        if weights is None:
            mdd = dd.mean(axis=1).values
            for _ in range(repeat):
                SE = (dd.values - mdd[:, None])**2
                SSE = np.nansum(SE, axis=1)
                valid = SSE > 0
                w = 1 / np.nanmean(SE[valid] / SSE[valid, None], axis=0)
                w = w * w.sum() / (w*w).sum() 
                # https://stats.stackexchange.com/questions/71081/degrees-of-freedom-for-a-weighted-average
                mvn['n'] = (dd.notnull() * w).sum(axis=1)
                m = mdd  = np.nansum(dd * w, axis=1) / mvn['n']
        else:
            if weights == False:
                w = np.ones((dd.shape[1],))
            else: 
                w = weights[idx]
                if paired: w = w * weights[idx+1]
            mvn['n'] = (dd.notnull() * w).sum(axis=1)
            m = mdd = np.nansum(dd * w, axis=1) / mvn['n']
        if paired or repeat==0: m = np.nansum(md * w, axis=1) / mvn['n']
        if paired:
            dx = data[idx+1] - data[idx]
            mdx = np.nansum(dx * w, axis=1) / mvn['n']
            SSE = np.nansum(((dx.values - mdx[:, None])**2), axis=1)
            mvn['mdx'] = pd.Series(mdx, index=dx.index)
        else:
            SSE = np.nansum(((dd.values - mdd[:, None])**2), axis=1)
        mvn['m'] = pd.Series(m, index=dd.index)
        # We do not use weighted variance due to its bias in the yeast WT data.
        # Low pops were resulted in when using weighted variance even when
        # V was adjusted by w.sum() - 1 instead of len(samples) - 1.
        # It is unclear whether variance of variances is biased.
        # But smaller ns would result in large expected sampling variance and 
        # much smaller prior variance and inflate significance.
        # For example that use weight only for calculating mean, see
        # http://www.analyticalgroup.com/download/WEIGHTED_MEAN.pdf
        mvn['ns'] = (dd.notnull()).sum(axis=1) # Number of Samples
        mvn['v' ] = SSE / (mvn['ns'] - 1) / scale_var
        mvn['En'] = mvn['ns']
        mvn['EV'] = mvn['v']
        
        mvn['m' ][mvn['ns']<2] = np.nan
        mvn['v' ][mvn['ns']<2] = np.nan
        mvn['ns'][mvn['ns']<2] = np.nan
        mvn['n' ][mvn['ns']<2] = np.nan
        mvns[idx] = mvn
    mvns['idxes'] = idxes[:len(idxes) - paired]
    return mvns

def calc_ttest_stats(mvns, paired, equalV, need_z = True):
    idxes = mvns['idxes']
    dxs = {}; pvals = {}; zs = {}
    if paired:
        for idx in idxes:
            mvn = mvns[idx]
            dx, v, n, En = mvn['mdx'], mvn['EV'], mvn['n'], mvn['En']
            t = dx / (v/n)**0.5
            p = ss.t.cdf(-abs(t), En-1) * 2
            pvals[idx] = pd.Series(p, index=t.index)
            if need_z:
                zs[idx] = pd.Series(np.sign(t) * -ss.norm.ppf(pvals[idx]/2), index=t.index)
            dxs[idx] = dx
    else: 
        for i in range(len(idxes[:-1])):
            idx = idxes[i];   idx1 = idxes[i+1]
            mvn0 = mvns[idx]; mvn1 = mvns[idx1]
            dx = mvn1['m'] - mvn0['m']
            v0, n0, En0 = mvn0['EV'], mvn0['n'], mvn0['En']
            v1, n1, En1 = mvn1['EV'], mvn1['n'], mvn1['En']
            p = ttest2(dx, v0, n0, En0, v1, n1, En1, equalV)
            pvals[idx] = pd.Series(p, index=dx.index)
            if need_z:
                zs[idx] = pd.Series(np.sign(dx) * -ss.norm.ppf(pvals[idx]/2), index=dx.index)
            dxs[idx] = dx
    return dxs, pvals, zs

def ttest(data, controls=None, paired=False, weights=None, equalV=False, 
          do_SVA=False, nSV=2, penaltyRatioToSV=0, normalize_min_pval=0.1, 
          ridge=True, parallel=True, need_z=True):
    idxes = data.columns.levels[0].tolist()
    if controls is None: 
        for _ in range(2): 
            data, controls = normalize_by_least_change(data, equalV=equalV, 
                                                       min_pval=normalize_min_pval, 
                                                       weights=weights, paired=paired)
    else:
#         mcs = data.loc[controls].mean(axis=1).dropna()
#         smc = sorted(mcs.values); l = len(smc)
#         mid_controls = mcs[(mcs >= smc[l//10]) & (mcs < smc[l//10*9])].index
        data = cross_group_normalization(data, idxes, controls)
    vbs = {}
    scale_var = 1
    if do_SVA: 
        data, SVApop, scale_var = partial_ridge_sva_residual(data, controls, nSV=nSV, penaltyRatioToSV=penaltyRatioToSV,
                                          ridge=ridge, parallel=parallel)
        vbs['SVApop'] = SVApop
    mvns = calc_mvns(data, weights, paired, scale_var)
    vbs['mvns'] = mvns
    dxs, pvals, zs = calc_ttest_stats(mvns, paired, equalV, need_z=need_z)
    if len(idxes)==2: 
        pvals = pvals[idxes[0]]; dxs = dxs[idxes[0]]
        if need_z: zs = zs[idxes[0]]
    vbs['dx'] = dxs
    return pvals, zs, vbs

def ascertained_ttest(data, controls=None, paired=False, weights=None, 
                      span=None, equalV=False, pre_neutralize=True, 
                      do_SVA=False, nSV=2, penaltyRatioToSV=0, normalize_min_pval=0.1, 
                      ridge=True, parallel=True):
    """
    data is a MultiIndex DataFrame with id (eg. shRNA, gene) as row index and different 
    experimental condition numbers (0, 1, 2, ...) and sample names in tuples as column index.
    eg.: 
>>> data.head()
                      0                             1            
                     S1        S2        S3        S6        S5        S7   
0610005C13Rik  1.058412       NaN  2.033172  1.283600  0.356335  1.330445   
0610007N19Rik       NaN  0.054081 -1.917866  2.215734  2.147257  2.367417   
0610007P14Rik  4.840528  4.656822  5.305128  5.117443  5.223823  4.544029   
0610008F07Rik  1.343370  1.813711  1.385292  1.711190  2.077108  1.935664   
0610009B14Rik  1.868023  2.189409  0.257601  0.774219  2.220796  2.242303   

    `data` has two experimental conditions (0 and 1), each with three replicates. 
    data should have already been normalized (with normal distribution; 
    Replicates should have the same mean or median).
     
    controls: the genes (or seqids) that are supposed to have neutral effects and 
    the mean of them will be used as controls.
    If controls is None, then the mean of all genes will be used as controls. 
    
    If paired is True, experimental conditions must have the same number of 
    replicates in matched order for a gene. 
    
    Return dicts:
    p values, z-scores
    """
    pops = {}
    idxes = data.columns.levels[0].tolist()
    tps, _, vbs = ttest(data, controls=controls, paired=paired, weights=weights, equalV=equalV, 
                        do_SVA=do_SVA, nSV=nSV, penaltyRatioToSV=penaltyRatioToSV,
                        normalize_min_pval = normalize_min_pval, need_z=False,
                        ridge=ridge, parallel=parallel)
    mvns = vbs['mvns']
    vbs['idxes'] = mvns['idxes']
    mns = [mvns[idx]['n'].mean() for idx in idxes[:len(idxes)-paired]]
    if len(idxes) == 2:
        _, _, pops[0] = neup(tps, with_plot=False, minr0=0, fine_tune=False)
    else:
        for i in range(len(idxes)-1): 
            _, _, pop = neup(tps, with_plot=False, minr0=0, fine_tune=False)
            pops[i] = pop
    if pre_neutralize: print(round(pops[0], 3), end='; ')
#         print("Student's t-test power of p-values:", list(pops.values()))
    base_pop = {}
    if span is None: span = 0.95 * min(1, 1 - 0.2 * (np.log10(len(mvns[0]['n']))-2))
    for i in range(len(mvns['idxes'])): 
        idx = idxes[i]
        i0 = max(0, i-1)
        pop = pops[i0]
        if 0 < i < len(pops)-1: pop = (pop + pops[i]) / 2
        mn = (mns[i0] + mns[i]) / 2
        if pre_neutralize == False: pop = 1
        elif len(idxes) == 2 and not paired:
            pop = pop ** (mns[i]/mn)
        pop = min(pop, 1)
        vbs[idx] = estimated_deviation(mvns[idx], pop, span=span, parallel=parallel)
        base_pop[idx] = (pop ** 0.5 * (mn - 1) + 1) / mn
        if paired: vbs[idx]['mdx'] = mvns[idx]['mdx']
    _, pvals, zs = calc_ttest_stats(vbs, paired, equalV)
#     for idx in pvals:
#         pvals[idx]  **= 1 / adj_pops[idx]
#         zs[idx] = np.sign(dxs[idx]) * -ss.norm.ppf(pvals[idx]/2)
            
    if len(idxes)==2: pvals = pvals[idxes[0]]; zs = zs[idxes[0]]
    vbs['mvns'] = mvns; vbs['ttest_pops'] = pops; vbs['base_pop'] = base_pop
    return pvals, zs, vbs
    
def ttest2(dx, v1, n1, En1, v2, n2, En2, equalV=False):
    if equalV:
        df = En1 + En2 - 2
        sp = np.sqrt( ( (En1-1)*v1 + (En2-1)*v2 ) / df )
        t = dx / sp / np.sqrt(1/n1 + 1/n2)
    else:
        vn1 = v1 / En1; vn2 = v2 / En2
        df = ((vn1 + vn2)**2) / ((vn1**2) / (En1 - 1) + (vn2**2) / (En2 - 1))
        t = dx / np.sqrt(v1/n1 + v2/n2)
    p = ss.t.cdf(-abs(t), df) * 2
    return p

def z_score_for_ttest_p_val(t, p):
    z = -itg.z_score(np.array(p)/2) * np.sign(t)
    return z

def cdfs2alV(s, n, V, a):
    return gammainc((n-1)/2, (n-1)*s**(1/a)/(2*V))

def means2a(n, V, a): 
    return np.exp(gammaln((n-1)/2+a) - gammaln((n-1)/2)  + a * np.log(2*V/(n-1)))

def vars2a(n, V, a, y0=0):
    G  = gammaln((n-1)/2)
    Gs = gammaln((n-1)/2+a)
    return (2*V/(n-1))**(2*a) * (np.exp(gammaln((n-1)/2+2*a)-G) - np.exp(2*(Gs-G))) - y0

# Oliphant, Travis E., "A Bayesian perspective on estimating mean, variance, and standard-deviation from data" (2006).All Faculty Publications. Paper 278.
# https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=1277&context=facpub
def dVsigma(n, r):
    return 1 - (n-3) / 2 * np.exp(gammaln(n/2-1) - gammaln((n-1)/2))**2 - r

def dfVls2 (V, s2, n, ElogV, VlogV0): 
    return (n-1) * s2 / 2 / V  -  (np.log(V) - ElogV) / VlogV0   -  (n+1)/2

def dfVls2a(V, s2, n, EVa, Vga, a, n0):  # @UnusedVariable
    return (n-1) * s2 / 2 / V  -  (V**a - EVa) * a * V**a / Vga  -  (n/2+1-a)

uselogV = False
fVls2 = itg.fVls2a; EVfV = itg.EVfVa
EsgmfV = itg.EsgmfVa; VsgmfV = itg.VsgmfVa
if uselogV:
    fVls2 = itg.fVls2; EVfV = itg.EVfV
    EsgmfV = itg.EsgmfV; VsgmfV = itg.VsgmfV

def calc_En_Es2(params):
    gene, n, V, EV0, VlogV0, a = params
    useNormal = False
    if uselogV: 
        s2 = max(V, EV0/100)
        args = [s2, n, np.log(EV0), VlogV0]
    else: 
        s2 = max(V, EV0/100)
        args = [s2, n, EV0**a, VlogV0, a, 0] # set n0 = 0
    En = n; EV = EV0
    if EV0 >= 0 and s2 >=0: # not NaN
        x0 = max(EV0, s2)
        while True:
            try: 
                if   uselogV  : Vmax = fsolve(dfVls2,  x0, args=(*args,))[0]
                elif useNormal: Vmax = fsolve(dfVls2a, x0, args=(*args,))[0]
                else:
                    rV = VlogV0 / EV0**(2*a)
                    n0 = 2*a**2/rV + 1
                    if rV >= 1e-3:
                        try: n0 = fsolve(vars2a, n0, args=(EV0, a, VlogV0))[0]
                        except:
                            try: n0 = fsolve(vars2a, 1+1e-11, args=(EV0, a, VlogV0))[0]
                            except: print('Cannot solve for n0. Use approximation', n0)
                    En = n0 + n - 1
                    EV = ((n0-1) * EV0 + (n-1) * s2) / (n0 + n - 2)
#                     return En, EV
                    Vmax = ((n0-n-2) + np.sqrt((n0-n-2)**2 + 4*(n0-1)/EV0*(n-1)*s2)) / (2 * (n0-1) / EV0)
                    args = [s2, n, EV0**a, EV0, a, n0]
                break
            except: 
                x0 *= 0.9
                if x0 < 0.5 * min(EV0, s2): 
                    print('x0={}; s2={}; n={}; ElogV={}; VlogV0={}'.format(
                        min(EV0, s2), s2, n, EV0, VlogV0)); raise
        pts = [EV0, s2, Vmax]
        args.append(Vmax)
        uplmt0 = upperlim = max(pts); lowerlim = min(pts)
        while fVls2(upperlim, *args) > 1e-20 and upperlim < 1000 * uplmt0: upperlim *= 1.1 + 10/n
        while fVls2(lowerlim, *args) > 1e-20: lowerlim /= 1.1 + 10/n
        try: 
            M   =  quad(fVls2, lowerlim, upperlim, args=(*args,), points=pts)[0]
        except: 
            print('EV0={}; s2={}; n={}; VlogV0={}; M={}; upperlim={}; lowerlim={}; pts={}'.format(
                    EV0, s2, n, VlogV0, M, upperlim, lowerlim, pts))
            raise
        if n < 100 or M > 0:
            EV   = quad(  EVfV, lowerlim, upperlim, args=(*args, M),       points=pts)[0]
            Esgm = quad(EsgmfV, lowerlim, upperlim, args=(*args, M),       points=pts)[0]
            Vsgm = quad(VsgmfV, lowerlim, upperlim, args=(*args, M, Esgm), points=pts)[0]
            if Vsgm / EV < 1e-3: En = np.e + 0.5 * EV / Vsgm
            else:
                try:  En = fsolve(dVsigma, 3.01, args=(Vsgm/EV,))[0]
                except: 
                    try: En = fsolve(dVsigma, 6, args=(Vsgm/EV,))[0]
                    except: 
                        print(gene, 'En Error: Vsgm={}; EV={}; EV0={}; s2={}; n={}; VlogV0={}; M={}; Vmax={}'.format(
                                    Vsgm, EV, EV0, s2, n, VlogV0, M, Vmax), params)
    return En, EV

def estimated_deviation(mvn, pop, span, parallel=True):
    Ens = {}; EVs = {}; vbs = {}
    genes = list(mvn['m'].index)
    
    n = mvn['ns']; v = mvn['v']; ms = np.array(mvn['m'])
    nn = 1 + (n - 1) * pop ** 0.5
    Vs = v * (n - 1) / n * nn / (nn-1)
    Ns = 1 + (mvn['n'] - 1) * pop ** 0.5
    Vars = Vs
    vbs['means'] = mvn['m']
    nn = np.array(nn); ns = np.array(mvn['ns'])
    Vs = np.array(Vs)
    x = np.array(nn)/np.nanmax(nn) + np.array(ms)/(np.nanmax(ms) - np.nanmin(ms))
    x[np.logical_or(Vs==0, np.isnan(Vs))] = np.NaN
    
    if uselogV: logVs = np.log(Vs); a = np.NaN
    else:
        a = segmented_mean_median_norm(x, Vs) 
        logVs = Vs ** a # logVs is actually Vs^a
    
    ElogVs  = loess_fit(x, logVs, span=span)
    for i in range(10):
        try:
            ElogVs2 = loess_fit(x, logVs, span=span*0.2*(i/2.+1))
            break
        except: pass
    ELoess  = loess_fit(x, (ElogVs - ElogVs2)**2, span=span)
    S2all   = loess_fit(x, (logVs  - ElogVs )**2, span=span)
    minS2all = np.nanmean(S2all) * 1e-3
    S2all[S2all<0] = 0
    S2all += minS2all
    meanns = loess_fit(x, ns, span=span) # Used to calculate V of sampling V (Vsmpl)
    # Variance of sampling variance should be the sun of sampling expression variance (VE) and 
    # sampling technical replicate variance (VT). Technical replicates would stabilize VT and thus
    # reduce Vsmpl. That is why we should not use nn to calculate Vsmpl. On the other hand,
    # technical replicates does not influence VE. But VE is difficult to estimate and could be small.

    """Variance normalization and correction. 
    Note: ONLY mean of s2 is an unbiased estimation of true variance; transformed s2 is not.
    https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
    """
    if uselogV: 
        EV0s = np.exp(ElogVs) / (means2a(meanns, 1, 1e-5)) ** (1/1e-5) # bias of log(s2) can be approximated by s2 with a small power
        Vsmpls = polygamma(1, (meanns-1)/2)
    else:
        EV0s = (ElogVs / means2a(meanns, 1, a)) ** (1/a) # Correct bias due to transformation
        Vsmpls = vars2a(meanns, EV0s, a)
    meanS2all = np.nanmean(S2all)
    rS2 = S2all / ElogVs**2; sortedrS2 = sorted(rS2[rS2>0])
    Vadj = Vsmpls * np.minimum(1, rS2 / mean(sortedrS2[len(sortedrS2)//10*9:]))
    meandif = np.nanmean(np.maximum(0, S2all - Vadj))
    VlogV0s = np.maximum(0, ELoess) + minS2all + meandif * S2all / meanS2all 
    if parallel:
        pool = Pool(psutil.cpu_count(logical=False))
        results = pool.map(calc_En_Es2, zip(genes, nn, Vs, EV0s, VlogV0s, [a]*len(nn)))
        pool.close()
    else:
        results = [calc_En_Es2(xi) for xi in zip(genes, nn, Vs, EV0s, VlogV0s, [a]*len(nn))]
    for i, (En, EV) in enumerate(results):
        gene = genes[i]
        Ens[gene] = En
        EVs[gene] = EV
    vbs['EV'] = EVs; vbs['En'] = Ens; vbs['n'] = Ns; vbs['m'] = mvn['m']
    
    vbs['S2all'] = S2all; vbs['Vs'] = Vars; vbs['mean'] = ms
    vbs['Vsmpls'] = Vsmpls; vbs['Vadj'] = Vadj; vbs['ElogVs'] = ElogVs
    vbs['EV0s'] = EV0s
    vbs['VlogV0s'] = VlogV0s; vbs['x'] = x; vbs['logVs'] = logVs
    for k in vbs: 
        if type(vbs[k]) is dict: vbs[k] = pd.Series(vbs[k])
    if not uselogV: vbs['a'] = a
    return vbs

def mean_median_norm(v, a0=0.3, only_positive=False):
    ''' 
    Normalize data by taking power a so that the difference between mean and median 
    (divided by std) is minimized. a0 is the initial guess of a. 
    If only_positive is True, then only positive (and non-zero) values are used 
    for the normalization. But all values are returned.
    
    Returns:
    tuple: Initial values raised by the optimal power a, and a itself.
    '''
    v0 = np.array(v); v = v0 + 0.
    if only_positive: v = v[v>0]
    else: v = v[v>=0]
    v.sort()
    v = v[len(v)//10 : len(v)//10*9]
    md = np.nanmedian(v)
    def loss(a): 
        va = v ** a
        return abs(np.nanmean(va) - md**a) / (np.nanstd(va) + 1e-19)
    a = minimize(loss, a0, bounds=[(0.001, 10)]).x[0]
    return v0**a, float(a)

def segmented_mean_median_norm(x, Vs):
    xs = np.array(x); xs.sort()
    es = []; l = len(x)
    steps = min(max(3, l//1000), 50)
    for i in range(steps):
        vi = Vs[(x>=xs[i*l//steps]) & (x<=xs[(i+1)*l//steps-1])]
        if len(vi)>1:
            _, e = mean_median_norm(vi, only_positive=True)
            es.append(e)
    a = np.median(es)
    return a