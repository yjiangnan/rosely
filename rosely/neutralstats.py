# -*- coding: utf-8 -*-
'''
This module implements the basic statistical tools for calculating neutrality-controlled
p-values and local false discovery rate (LFDR). 

Created on Aug 9, 2017

@author: Jiang-Nan Yang
'''


import time
import scipy.stats as ss, pandas as pd
from pylab import *
from scipy.optimize import curve_fit, minimize
# from scipy.integrate import quad # General purpose integration
from sklearn.isotonic import IsotonicRegression
# from scipy.interpolate import pchip_interpolate
from scipy.optimize import OptimizeWarning
from scipy.integrate import IntegrationWarning
warnings.filterwarnings("ignore", category=OptimizeWarning)
warnings.filterwarnings("ignore", category=IntegrationWarning)

__all__ = ['dicthist', 'false_discovery_rate', 'locfdr', 'neup', 'density', 'poisson_regression']

def dicthist(data, a=1, nbins=20):
    '''Plot a histogram for a dict.'''
    if type(data) is dict: v = array(list(data.values()))
    else: v = array(data)
    v = v[abs(v)>=0]
    return hist(v**a, nbins)

def false_discovery_rate(ps):
    ranks0 = ss.mstats.rankdata(np.ma.masked_invalid(ps))
    ranks = ranks0 + 0
    ranks[ranks==0] = np.nan
    FDR = ps / ranks * np.nanmax(ranks)
    FDR[FDR>1] = 1
    fdrs = np.array(FDR).tolist()
    rks = np.array(ps).argsort().argsort() # Two argsort produces rank
    zipped = sorted(zip(rks, range(len(rks))))[::-1]
    for j, (_, i) in enumerate(zipped[:-1]):
        i1 = zipped[j+1][1]
        if fdrs[i] < fdrs[i1]: fdrs[i1] = fdrs[i]
    FDR[:] = fdrs
    return FDR

def density(x, nbins=50):
    x = array(x); x = x[abs(x)>=0]
    bins = linspace(x.min(), x.max(), nbins)
    bin_nums = np.digitize(x, bins) - 1
    dx = bins[1] - bins[0]
    freqs = np.float64(np.bincount(bin_nums)); freqs /= sum(freqs) * dx
    bins += dx/2
    return bins, freqs

# @profile
def poisson_regression(z, bins, freqs, df_fit, symmetric=True, lowerfirst=True, 
                       method="SLSQP"): # SLSQP is fast and stable. Default method 'BFGS' is not deterministic.
    dx = bins[1] - bins[0]
    if symmetric: 
        x = bins[(bins)>=-dx/2]; y0 = freqs[(bins)>=-dx/2]
        ir = IsotonicRegression(increasing=False)
        y = ir.fit_transform(x, y0)
    else: x = bins[abs(bins) >= 0]; y = y0 = freqs[abs(bins) >=0]
    def powered(x):
        xs = [ones((len(x),))]
        for _ in range(df_fit-1): xs.append(xs[-1] * x)
        return array(xs)
    X = powered(x)
    wx = exp(-(x*0.333)**2)
    def combpoly(a):
        p1 = dot(a[:df], Xdf)
        p2 = dot(a[df:], Xdf)
        return (p1-p2) * wx + p2
    def fun_poisson(a, creg=1e-4): # https://en.wikipedia.org/wiki/Poisson_regression
        p = combpoly(a)
        reg = (abs(a) *  wreg).sum() / (2*df) * creg
        return -(p * y - exp(p)).sum()/len(x) + reg
    if lowerfirst:
        df = (df_fit+1)%2 + 1
        a0 = [0] * (df*2)
    else: a0 = [0] * (df_fit*2)
    while len(a0) <= df_fit * 2:
        df = len(a0)//2
        Xdf = X[:df, :]  # @UnusedVariable
        wreg = array(list(array(range(1, df+1))**2) * 2)  
        # regularize a, especially for coefficients of higher degrees
        res = minimize(fun_poisson, a0, method=method)
        a0 = list(res.x)[:df] + [0, 0] + list(res.x)[df:] + [0, 0]
    df = df_fit
    wx = exp(-(z*0.333)**2)
    Xdf = powered(z)
    Ey = exp(combpoly(res.x))
    return Ey

def locfdr(ps=None, zs=None, p0s=[0, 0.1], zthreshold=1.75, nbins = 50, df_fit=5, 
           data_name='data', with_plot=True, new_figure=True, ax0=121, axes=[], plot_density_only=False,
           lowerfirst=True, poisreg_method='SLSQP'):
    """ Calculate the local false discovery rate (LFDR) from a list or dict or pandas Series of pvalues ps or z-scores zs.
    p0s: The initial guess of gaussian parameters for neutral genes. The important part is not the values themselves,
         but the number of values. p0s=[0] would assume a standard normal distribution for neutral genes, and
         p0s=[0, 0.1] would fit the variance of neutral genes (important if ps and zs are not neutrally controlled). 
    zthreshold: absolute z scores smaller than zthreshold will be used as neutral gene values to fit the gaussian curve.
    nbins: the number of bins for the histogram of frequency distribution for z scores. 
           Use smaller values if the number of data is small (a few hundreds) and the curves do not fit well.
    df_fit: degree of freedom for the polynomial fit of density curve of all data.
            Use small values (e.g. 5) if number of data (and nbins) is small.
    data_name: the name displayed for the data dots in the density function plot.
    with_plot: whether or not to plot the data, curves and LFDRs.
    
    Return: the caculated LFDR (either an array or a pandas Series, depending on the input)
    """
    if ps is not None:
        if type(ps) is dict: 
            zvals0 = ss.norm.ppf([ps[sid]/2 for sid in sorted(ps)])
        elif type(ps) is list: zvals0 = ss.norm.ppf(array(ps)/2)
        elif type(ps) is np.ndarray: zvals0 = ss.norm.ppf(ps/2)
        elif type(ps) is pd.Series: zvals0 = pd.Series(ss.norm.ppf(ps/2), index=ps.index)
        else: raise Exception('data type for ps only supports dict, list, numpy.ndarray, or pandas.Series')
    elif zs is not None:
        ps = zs
        if type(zs) is dict: 
            zvals0 = array([zs[sid] for sid in sorted(zs)])
        elif type(zs) is list: zvals0 = array(zs) + 0.
        else: zvals0 = zs + 0.
    # zvals0 = array(list((zs.values())))
    MAXZ = 30
    zvals0[zvals0==Inf] = MAXZ+1; zvals0[zvals0==-Inf] = -MAXZ-1
    try: minz = nanmin(zvals0[zvals0>-MAXZ]); maxz = nanmax(zvals0[zvals0<MAXZ])
    except: print('Error in locfdr:', zvals0, ps); raise
    zvals0[zvals0> MAXZ] =  (zvals0[zvals0> MAXZ] - maxz)/10 + maxz
    zvals0[zvals0<-MAXZ] =  (zvals0[zvals0<-MAXZ] - minz)/10 + minz
    zvals = zvals0[logical_and(abs(zvals0)>=0, abs(zvals0)<MAXZ)]
    zvals = np.sort(np.append(-abs(zvals), abs(zvals)))
    nbins = int(nbins * (nanmax(abs(zvals)/2))**0.5 * (len(zvals) / 1000)**0.5)
#     print('nbins:', nbins, 'len(zvals):', len(zvals))
    bins, freqs = density( zvals, nbins )

    x = bins[abs(bins)<zthreshold]
    y = freqs[abs(bins)<zthreshold]
    def gaus(x, a, v=1):
        return exp(a/2 - x**2 / (4*v))
    para0 = curve_fit(gaus, x, y**0.5, p0=p0s)[0]
    if with_plot:
        if new_figure: figure(figsize=(12,5))
        if len(axes): sca(axes[0])
        else: subplot(ax0)
        plot(bins, (freqs), 'b.', label=data_name)
        plot(bins, gaus(bins, *para0)**2, 'r', label='fit normal')
    
    den_all = poisson_regression(append(abs(zvals0), abs(zvals)), bins, freqs, df_fit, 
                                 lowerfirst=lowerfirst, method=poisreg_method)
    LFDR = gaus(abs(zvals0), *para0)**2 / den_all[:len(zvals0)]
    LFDR[LFDR>1] = 1; zpeak = min(3, max(abs(zvals0[LFDR==nanmax(LFDR)])))
    LFDR[abs(zvals0)<=zpeak] = max(LFDR[abs(zvals0)<=zpeak])
    LFDR /= nanmax(LFDR)
    idx = LFDR >= 0; x = abs(zvals0)[idx]; y = LFDR[idx]
    ir = IsotonicRegression(increasing=False)
    LFDR[idx] = ir.fit_transform(x, y)
    retLFDR = LFDR
    if with_plot:
        plot(zvals, den_all[len(zvals0):], 'g', label='fit all')
#         plot(x, y, '.y', label='isotonic')
        xlabel('± z score', fontweight='bold'); ylabel('probability density', fontweight='bold'); #ylim([0, ylim()[1]])
        legend()
        if not plot_density_only:
            if len(axes): sca(axes[1])
            else: subplot(ax0+1)
            pvals = ss.norm.cdf(-abs(zvals0))*2
            idx = pvals.argsort()
            x = ss.rankdata(pvals[idx])/len(pvals)
            plot(x, pvals[idx], color='orange', label='p value')
            plot(x, LFDR[idx], 'r.', label='LFDR')
            plot(x, false_discovery_rate(pvals[idx]), label='FDR')
            xlabel('rank proportion', fontweight='bold'); ylabel('FDR or $p$ value', fontweight='bold')
            legend(); grid('on')
    if type(ps) is dict:
        retLFDR = {}
        for i, sid in enumerate(sorted(ps)): retLFDR[sid] = LFDR[i]
    return retLFDR
    
def fit_pvalues(x, y, pvals, aa=1):
    def fun(logx, y, *a):
        return (a[0]**(1. + a[1] * y ** a[2])) * logx
    logx = log(x); logx[logx<-741] = -741
    para0 = [aa, 0., 2]
    found = False
    while not found:
        logy = log(y)
        res = minimize(lambda b: sum((fun(logx, y, *b) - logy)**2)/len(y) * (1 + 0.2 * b[1]**2), 
                       para0, bounds=([0.0, 2], [-1, 1], [1.5, 5]), method="SLSQP")
        para = list(res.x)
#         for method in ['trf', 'dogbox']:
#             try: 
#                 para = curve_fit(lambda logx, *para: fun(logx, y, *para), logx, log(y), 
#                                   p0=tuple(para0), method=method, bounds=([0.0, -1, 1.5], [2, 1, 5]))[0]
        found = True; break
#             except KeyboardInterrupt: return
#             except: pass
        if not found: 
            para0 = [aa * (0.5 + np.random.uniform()), np.random.uniform(), np.random.uniform(1.5, 3)]
            print('neup: Solution to p value fit not found. Trying new para0:', para0)
    if para[1] * para[2] > 0 or 1: 
        rk = ss.mstats.rankdata(np.ma.masked_invalid(pvals))
        if aa > 1 and para[0] > 1: para[0] **= 0.5 
        psneu = exp(fun(log(pvals), rk/max(rk), *para))
    return psneu, para

def fp(p, a, pa):
    return p ** (a ** (1 + pa[1] * p + pa[2] * p*p))

def lsq_fit_pvalues(r, p0, pvals):
#     w = 1 - r
    logr = log(r)
    logp = polyval((polyfit(logr, log(p0), 5)), logr)
    logp -= logp.max() + abs(logp[-1] - logp[-2]) * -logr[-1]/(logr[-1] - logr[-2])
    p = exp(logp)
    y = log(logr / (logp-1e-16))
    X = transpose([ones(len(p0)), p, p*p])
    XTX = X.T @ X
    XTX[1, 1] *= 1.3
    XTX[2, 2] *= 1.1
    para = inv(XTX) @ X.T @ y
    a = exp(para[0])
    if isnan(a):
        raise
    para /= para[0]
    return fp(pvals, a, para), a, para
    

def neup(pvals0, getLFDR=True, LFDRthr0 = 0.5, minr0=0.0, nbins=50, df_fit=5, data_name='data', 
         with_plot=True, fine_tune=True, new_figure=True, plotLFDR=True, ax0=121, axes=[], global_changes=False,
         poisreg_method='SLSQP', base_pop=1, subplot_label=''):
    """ Calculate the neutrality-controlled p values.
    getLFDR: whether or not to calculate local false discovery rate (LFDR).
    LFDRthr0: the threshold of LFDR to consider as the boundary between neutral and non-neutral observations
             if most observations have LFDR of 1 and p-values are random. 
    Other parameters are passed to locfdr()
    
    Returns:
    neutralized p values, LFDRs, power of raw p values
    """
    if type(pvals0) is dict:
        pvals = pd.Series(pvals0)
    elif type(pvals0) is list: pvals = array(pvals0)
    else: pvals = pvals0
    ps0 = ps = array(pvals[pvals<=1])
    ps.sort()
    aa = 1; a0 = 1; rp = np.arange(1, len(ps) + 1) / len(ps); r00 = 0; SE0 = Inf; r0 = re = 1
    for i in range(1000):
        try: lfdr = locfdr(ps, p0s=[0], with_plot=False, df_fit=df_fit//2+1, 
                           zthreshold=1.5, lowerfirst=False, poisreg_method=poisreg_method)
        except: print(r0, re, ps0, pvals); raise
        LFDRthr = LFDRthr0 * mean(lfdr) * aa**((1. - aa)*0.5)
        if lfdr.min() < LFDRthr: r0 = min(r0, min(0.9, rp[(lfdr<LFDRthr)][-1]))
        else: r0 = 0.
        if r00 and global_changes: r0 = (r0*r00**9) ** 0.1
        r00 = r0
        r0 = max(minr0, r0)
        pe = LFDRthr0/2
        re = r0 + (1 - r0) * pe
        x = ps0[int(len(ps) * re):]
        y = linspace(pe, 1, len(x)) - 0.5*(1 - pe)/(len(x) - 1)
        
#         w = y * y
#         logX = transpose([log(x)*w]); logY = log(y)*w
#         try: a = linalg.lstsq(logX, logY, rcond=None)[0][0]
#         except: a = linalg.lstsq(logX, logY)[0][0]
#         ps, param = fit_pvalues(x, y, ps0, aa)
#         a = param[0]

        ps, a, para = lsq_fit_pvalues(y, x, ps0)
        
        if (a-1)*(a-a0) <= 0: break
        aa = a
        if abs(a-1) < 1e-9: break
        a0 = a
        if global_changes: 
            yhat = fp(x, a, para)
            SE = ( sum((yhat - y)**2 * y*y)/sum(y*y) )**0.5 / len(x) * aa
        else: SE = 0
        if with_plot: 
            print('iterations:{} r0:{:.3f}, re:{:.3f}, a:{:.5f} lsq para:{}'.format(
                i, r0, re, a, float32(para[1:])) + ' SE:{:.7f}'.format(SE)*global_changes)
        if (SE > SE0 and global_changes) or r0==0: break
        SE0 = SE
    if aa > 1: aa **= 0.5
    fine_tuned = False
    psneu = fp(pvals, aa, para)
    if getLFDR:
        if fine_tune:
            x = ps0[int(len(ps) * re):]
            psneu, para = fit_pvalues(x, y, pvals, aa)
            fine_tuned = True
            if with_plot: print('finetune para:', para, ' raw aa:', aa)
            aa = para[0]
        LFDR = locfdr(psneu, p0s=[0.], with_plot=with_plot and plotLFDR, 
                      new_figure=new_figure, plot_density_only=True, ax0=ax0, axes=axes,
                      data_name=data_name, nbins=nbins, df_fit=df_fit, poisreg_method=poisreg_method)
    else: LFDR = psneu
    aa *= base_pop
    if with_plot:
        print('total iterations:{} r0:{:.3f}, re:{:.3f}, last a:{:.5f}'.format(i, r0, re, a))
        print('power of p-values: fit: {:.4f}, base: {:.4f}, final:{:.4f}'.format(aa/base_pop, base_pop, aa))
        if aa > 1.5:
            warnings.warn('''
        The power of p values is too large (ideally the power should be 1), meaning that most raw p values
        are much larger than that expected by random chance and the controls failed to serve as random data. 
        This could happen either by statistical underestimation of degrees of freedom, 
        or by stable noise in most data that blurred the difference between groups,
        or by different samples not having the same data quality,
        or by problems in data normalization. The calculated 
        LFDRs (adjusted by the power) should be treated with caution and could underestimate true LFDRs.
            ''')
            time.sleep(0.1)
        sortedpsneu = array(sorted((psneu[psneu>=0])))
        x = (arange(len(sortedpsneu))/len(sortedpsneu) + 1./len(sortedpsneu))
        rawFDR = sortedpsneu / x; FDRplot = LFDR[LFDR<=1]
        xunif = [r0, 1]; yunif = [1./len(x), 1]
        if subplot_label and plotLFDR:
            if len(axes): sca(axes[0])
            else: subplot(ax0)
            text(-3.95, 0.39, subplot_label[0], fontsize=24)
        if plotLFDR: 
            if len(axes): sca(axes[-1])
            else: subplot(ax0+1)
        cla()
        i = 0
        if i: 
            x = log10(x); ps0 = log10(ps0); sortedpsneu = log10(sortedpsneu); 
            rawFDR = log10(rawFDR); FDRplot = log10(FDRplot); xunif = log10(xunif); yunif = log10(yunif)
        plot(xunif, yunif, label='neutral')
        plot(x, ps0**(1/base_pop), label='$p$-val')
        plot(x, sortedpsneu, label='~'*int(fine_tuned) + '$p$-val$ ^{{{:.3f}}} $'.format(aa))
#             plot(x, rawFDR, '.', label='rawFDR')
        plot(x, sorted(FDRplot), '.', label='LFDR')
        plot([r0, r0], [0, 1], ':')
        grid('on'); legend(loc='lower right'); 
        xlabel('log$_{10}$ '*i + 'rank proportion', fontweight='bold'); 
        ylabel('log$_{10}$ ('*i + '$p$ value or LFDR' + ')'*i, fontweight='bold')
        if subplot_label: text(-0.03, 0.96, subplot_label[int(plotLFDR)], fontsize=24)
    return psneu, LFDR, aa
