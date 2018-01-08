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
    for i, rk in enumerate(ranks0):
        fdr = FDR[i]
        if rk: FDR[logical_and(ranks < rk, FDR > fdr)] = fdr
    return FDR

def density(x, nbins=30):
    x = array(x); x = x[abs(x)>=0]
    bins = linspace(min(x), max(x), nbins)
    bin_nums = np.digitize(x, bins) - 1
    dx = bins[1] - bins[0]
    freqs = np.float64(np.bincount(bin_nums)); freqs /= sum(freqs) * dx
    bins += dx/2
    return bins, freqs

def poisson_regression(z, bins, freqs, df_fit, symmetric=True):
    dx = bins[1] - bins[0]
    if symmetric: 
        x = bins[(bins)>=-dx/2]; y0 = freqs[(bins)>=-dx/2]
        ir = IsotonicRegression(increasing=False)
        y = ir.fit_transform(x, y0)
    else: x = bins[abs(bins) >= 0]; y = y0 = freqs[abs(bins) >=0 ]
    xcenter = x[y.argmax()]
    def fun_poisson(a, x, y): # https://en.wikipedia.org/wiki/Poisson_regression
        v = polyval(a, x)
        return -sum( (v * y - exp(v)) * ((abs(x - xcenter) * (logical_and(y0>0, abs(x - xcenter) < 5) + 0.5) + 1))**2 )
    a0 = zeros((df_fit,))
    res = minimize(fun_poisson, a0, args=(x, y))
    return exp(polyval(res.x, z))

def locfdr(ps=None, zs=None, p0s=[0, 0.1], zthreshold=1.75, nbins = 30, df_fit=5, data_name='data', with_plot=True):
    """ Calculate the local false discovery rate (LFDR) from a list or dict of pvalues ps or z-scores zs.
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
    
    Return: the caculated LFDR (either an array or a dict, depending on the input)
    """
    if ps is not None:
        if type(ps) is dict: 
            zvals0 = ss.norm.ppf([ps[sid]/2 for sid in sorted(ps)])
        else: zvals0 = ss.norm.ppf(array(ps)/2)
    elif zs is not None:
        ps = zs
        if type(zs) is dict: 
            zvals0 = array([zs[sid] for sid in sorted(zs)])
        else: zvals0 = array(zs) + 0.
    # zvals0 = array(list((zs.values())))
    zvals0[zvals0==Inf] = 21; zvals0[zvals0==-Inf] = -21
    try: minz = nanmin(zvals0[zvals0>-20]); maxz = nanmax(zvals0[zvals0< 20])
    except: print(zvals0, ps); raise
    zvals0[zvals0> 20] =  (zvals0[zvals0> 20] - maxz)/40 + maxz
    zvals0[zvals0<-20] =  (zvals0[zvals0<-20] - minz)/40 + minz
    zvals = zvals0[logical_and(abs(zvals0)>=0, abs(zvals0)<20)]
    zvals = np.sort(np.append(-abs(zvals), abs(zvals)))
    bins, freqs = density( zvals, int(nbins * (nanmax(abs(zvals)) / 2)**0.5 * (len(zvals) / 1000)**0.5) )

    x = bins[abs(bins)<zthreshold]
    y = freqs[abs(bins)<zthreshold]
    x = x[y>0]; y = y[y>0]
    def gaus(x, a, v=1):
        return a - x**2 / (2*v)
    para0 = curve_fit(gaus, x, log(y), p0=p0s)[0]
    if with_plot:
        figure(figsize=(12,5))
        subplot(121)
        plot(bins, (freqs), 'b.', label=data_name)
        plot(bins, exp(gaus(bins, *para0)), 'r', label='fit normal')
    
    den_all = poisson_regression(append(abs(zvals0), abs(zvals)), bins, freqs, df_fit)
    LFDR = exp(gaus(abs(zvals0), *para0)) / den_all[:len(zvals0)]
    LFDR[LFDR>1] = 1; zpeak = max(abs(zvals0[LFDR==nanmax(LFDR)]))
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
        subplot(122)
        pvals = ss.norm.cdf(-abs(zvals0))*2
        idx = pvals.argsort()
        x = ss.rankdata(pvals[idx])/len(pvals)
        plot(x, pvals[idx], color='orange', label='p value')
        plot(x, LFDR[idx], 'r.', label='Local FDR')
        plot(x, false_discovery_rate(pvals[idx]), label='FDR')
        xlabel('rank proportion', fontweight='bold'); ylabel('FDR or $p$ value', fontweight='bold')
        legend(); grid('on')
    if type(ps) is dict:
        retLFDR = {}
        for i, sid in enumerate(sorted(ps)): retLFDR[sid] = LFDR[i]
    return retLFDR
    

def neup(pvals0, getLFDR=True, LFDRthr0 = 0.5, minr0=0.00001, nbins=30, df_fit=5, data_name='data', 
         with_plot=True, fine_tune=True, force_fine_tune=False):
    """ Calculate the neutrality-controlled p values.
    getLFDR: whether or not to calculate local false discovery rate (LFDR).
    LFDRthr0: the threshold of LFDR to consider as the boundary between neutral and non-neutral observations
             if most observations have FDR of 1 and p-values are random. 
    Other parameters are passed to locfdr()
    
    Returns:
    neutralized p values, LFDRs, power of raw p values
    """
    if type(pvals0) is dict:
        names = sorted(pvals0)
        pvals = array([pvals0[n] for n in names])
    else: pvals = array(pvals0)
    ps0 = ps = array(sorted((pvals[pvals<=1])))
    aa = 1; a0 = 1; rp = np.arange(1, len(ps) + 1) / len(ps); r00 = 0; SE0 = Inf; r0 = re = 0
    for i in range(1000):
        try: FDR = locfdr(ps, p0s=[0], with_plot=False, zthreshold=1.5)
        except: print(r0, re, ps0, pvals); raise
        FDRthr = LFDRthr0 * mean(FDR) #* aa**0.3
        if FDR.min() < FDRthr: r0 = min(0.8, rp[(FDR<FDRthr)][-1])
        else: r0 = 0.
        if r00 and minr0: r0 = (r0*r00**9) ** 0.1
        r00 = r0
        r0 = max(minr0, r0)
        pe = LFDRthr0/2
        re = r0 + (1 - r0) * pe
        x = ps[int(len(ps) * re):]
        y = linspace(pe, 1, len(x)) + 1./len(x)
        a = linalg.lstsq(transpose([log(x)*y*y]), log(y)*y*y)[0][0]
#         a = dot(log(y), linalg.pinv([log(x)]))[0]
        if (1-a)*(1-a0) < 0: aa *= a; break
        aa *= a
        ps = ps0**(aa)
        if abs(a-1) < 1e-9: break
        a0 = a
        SE = ( sum((x**a - y)**2 * y*y)/sum(y*y) )**0.5 / len(x)
        if with_plot: 
            print('iterations:{} r0:{:.3f}, re:{:.3f}, aa:{:.5f} a:{:.5f} SE:{:.7f}'.format(i, r0, re, aa, a, SE))
        if SE > SE0 and minr0: break
        SE0 = SE
    psneu = pvals ** (aa); fine_tuned = False
    if getLFDR:
        if (aa < 0.5 or force_fine_tune) and fine_tune:
            def fun(logx, y, *a):
                return (a[0]**(1. - a[1] * y ** a[2])) * logx
            x = ps0[int(len(ps) * re):]
            logx = log(x); logx[logx<-741] = -741
            para0 = [aa, 0.3, 1.3]
            found = False
            while not found:
                for method in ['lm', 'trf', 'dogbox']:
                    try: 
                        para = curve_fit(lambda logx, *para: fun(logx, y, *para), logx, log(y), 
                                          p0=tuple(para0), method=method, bounds=([0, 0, 0], [2*aa, 1, 5]))[0]
                        found = True; break
                    except KeyboardInterrupt: return
                    except: pass
                if not found: 
                    para0 = [aa * np.random.uniform() * 2, np.random.uniform(), np.random.uniform(0, 2)]
                    print('neup: Solution to p value fit not found. Trying new para0:', para0)
            if para[1] * para[2] > 0 or 1: 
                fine_tuned = True
                rk = ss.mstats.rankdata(np.ma.masked_invalid(pvals))
                psneu = exp(fun(log(pvals), rk/max(rk), *para))
#             aa = para[0]
        LFDR = locfdr(psneu, p0s=[0.], with_plot=with_plot, data_name=data_name, nbins=nbins, df_fit=df_fit)
    else: LFDR = psneu
    if with_plot:
        print('total iterations:{} r0:{:.3f}, re:{:.3f}, last a:{:.5f}'.format(i, r0, re, a))
        print('power of pval:', round(aa, 4))
        if fine_tuned: print('pvalue fit para:', para)
        if aa > 1.5:
            warnings.warn('''
        The power of p values is too large (ideally the power should be 1), meaning that most raw p values
        are much larger than that expected by random chance and the controls failed to serve as random data. 
        This could happen either by statistical underestimation of degrees of freedom, 
        or by stable noise in most data that blurred the difference between groups. The calculated 
        LFDRs (adjusted by the power) should be treated with caution and could underestimate true LFDRs.
            ''')
            time.sleep(0.1)
        sortedpsneu = array(sorted((psneu[psneu>=0])))
        x = (arange(len(sortedpsneu))/len(sortedpsneu) + 1./len(sortedpsneu))
        rawFDR = sortedpsneu / x; FDRplot = LFDR[LFDR<=1]
        xunif = [r0, 1]; yunif = [1./len(x), 1]
#         figure(figsize=(14,5))
        for i in range(1):
            subplot(122)
            cla()
            if i: 
                x = log10(x); ps0 = log10(ps0); sortedpsneu = log10(sortedpsneu); 
                rawFDR = log10(rawFDR); FDRplot = log10(FDRplot); xunif = log10(xunif); yunif = log10(yunif)
            plot(xunif, yunif, label='neutral')
            plot(x, ps0, label='$p$-val')
            plot(x, sortedpsneu, label='~'*int(fine_tuned) + '$p$-val$ ^{{{:.3f}}} $'.format(aa))
#             plot(x, rawFDR, '.', label='rawFDR')
            plot(x, sorted(FDRplot), '.', label='LFDR')
            grid('on'); legend(loc='lower right'); 
            xlabel('log$_{10}$ '*i + 'rank proportion'); ylabel('log$_{10}$ ('*i + '$p$ value or LFDR' + ')'*i)
    if type(pvals0) is dict:
        psdict = {}; FDRdict = {}
        for (i, n) in enumerate(names):
            psdict[n] = psneu[i]; FDRdict[n] = LFDR[i]
        psneu = psdict; LFDR = FDRdict
    return pd.Series(psneu), pd.Series(LFDR), aa
