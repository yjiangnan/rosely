# Rosely
### Robust and sensitive tools for high-throughput data analysis
Rosely implements tools in Python for data normalization, identification of significantly changed genes/transcripts/constructs in microarray, RNA-seq, or shRNA/CRISPR library data, and downstream Gene Ontology and KEGG pathway enrichment analysis. It puts central emphasis on robustness and gives you *p*-values that are tested and can be supported by the data. Even better, this does not come with the cost of reduced sensitivity. 

For usage examples, see [Rosely examples](https://github.com/yjiangnan/rosely/tree/master/examples).

# Table of Contents
1. [Background](#Background)
2. [Practical Guide](#Practical)
3. [Installation](#Installation)

## Background<a name="Background"></a>

The combination of high-throughput data and modern statistical tools has enabled sensitive discoveries of important genes/pathways involved in biological processes. However, important questions remained unanswered. How reliable are these discoveries? When a statistical tool give you numbers of false discovery rates (FDRs), how can you trust these numbers? Are there biases and inflated false significance in the results? The reproducibility crisis hotly discussed in recent years in biological research indicates that these questions are indeed realistic concerns. There are no reasons why the same researchers, the same biological samples, and the same experimental procedures would not bring the same reproducibility crisis just because of some fancy high-throughput techniques and statistical tools that are boosted mostly for sensitivity instead of robustness. 

The situation may even be worse in high-throughput data analysis because the sophistication of the methods of both the experimental part and the statistical analysis part usually results in the complete separation/isolation of the two parts into the hands of different people, with the biologists rarely fully understand the statistical consequences of some easily ignored but critical actions and the bioinformaticians rarely fully understand what was going on in the experiments to formulate realistic proper assumptions to carry out the analysis. People rarely appreciate the pressure, frustration, and even guilt of the bioinformaticians when it becomes difficult for them to present very excitingly possitive results to the biologists, which encourages the bioinformaticians to come up with clever ways to produce at least a decent number of discoveries and explains why most statistical tools focus on sensitivity instead of robustness. On the other hand, biologists usually only choose the most beautiful results (possibly from different analyses or even different bioinformaticians sometimes) to present in their high-impact papers. The papers may look nice, but the devil would eventually come back at the researchers and eat their life because the results point to the wrong direction and get them lost in the so many genes/pathways identified. The worst part is, the results cannot be called into question because of their cost, sophistication, and their often assumed unbiasedness. 

`Rosely` is developed with the emphasis of robustness in mind. It is developed for the purpose of parallel studies of many genes to give a reliable confidence level for each gene instead of just identifying some top candidates. So, there is no motivation for false significance but strong desire against it. I personally carried out most part of the years of tedious experiments (shRNA library production; lenti-virus production; stem cell isolation, culture, transduction, and transplatation; generation of libraries ready for sequencing), not just to follow standard protocols, but made substantial improvements in the experimental procedures. This allows me to deeply think and understand the many problems in biological processes/experiments and their statistical consequences. Therefore, `Rosely` does not simply give you a list of *p*-values and FDRs, it is optimized for reducing the long-term cost of researchers and at the same time providing reliable sensitive discoveries. 

A central component of `Rosely` is the development of neutrality-controlled *p*-values and local false discovery rates (LFDRs). Neutrality-controlled *p*-values combat against almost all types of false significance and do not allow inflated significance get reported without raising eyebrows and requiring strong justifications. How it works is that, with some simplifications, all *p*-values must be transformed by a power (*pop*) to get the final neutrality-controlled *p*-values so that sorted non-significant *p*-values must approach a straight line expected by random chance. LFDRs are calculated based neutrality-controlled *p*-values, with great improvements of stability from [Efron (2007)'s](https://projecteuclid.org/download/pdfview_1/euclid.aos/1188405614) inspiring method. Note: conventional false discovery rate (FDR or *q*-value) is a measurement of false discoveries for all the genes/constructs whose *p*-values are smaller than a gene/construct; it says nothing about the probability of false discovery for the gene/construct itself.

Ascertained *t*-test is develop to achieve conservative but sensitive discoveries. It is slightly more conservative than Welch's *t*-test for false discoveries, but only slightly less sensitive than *z*-test with known population variance for true discoveries in normally distributed simulation data. Its method is similar to, and inspired by, moderated *t*-test in `limma`, but with great emphasis on not overestimating statistical significance. It strictly and completely follows Baysian statistics instead of some simple approximations to accurately calculate *p*-values so that one can more easily focus on its assumptions and data quality in case of bias in *p*-values (e.g. low *pop*). 

`Rosely` also makes serious efforts for data normalization. It implements partial-ridge-regression-regularized surrogate variable analysis (SVA) to conservatively remove the effects of hidden variables. It separates within-group sample normalization for variance reduction and cross-group normalization to guarantee that a set of control genes do not change over conditions (which could easily arise as an artifact when the data distribution across experimental groups change). 

## Practical Guide<a name="Practical"></a>
For usage examples, see [Rosely examples](https://github.com/yjiangnan/rosely/tree/master/examples).

`Rosely` also implements the following tools:
### 1. Classes for reading and normalizing data and for results
`ReadCounts` is a class for reading excel and text files containing read counts for each gene/constructs in each sample. It implements data normalization (to bring different samples to the same level) by either a specific control gene/construct, or the median of all data, or a quantile (1-9) of all data around the median, or a list of control genes/constructs.

`Result` is a class to wrap the results of raw *p*-values, *z*-scores, their neutralized values, local false discovery rates (and a member function for recalculation from raw *p*-values with different parameters) and the changes of normalized values. It has a member `results` of type of `pandas DataFrame` to utilize `pandas` tools for further processing.

### 2. Ascertained *t*-test
`ascertained_ttest` calculates *p*-values (with high statistical significance) of the changes of genes/constructs in different experimental groups using the similarity of genes/constructs in the same experimental group to reduce the uncertainty of the variance of each gene/construct in a rigorous empirical Bayesian approach. It has careful consideration of NOT overestimating the degrees of freedom and statistical significance. It outperforms mainstream RNA-seq analysis tools edgeR and DESeq2 for controlling false discovery rate in comparison of WT vs. WT from 48 yeast replicates. It uses LOWESS regression (locally weighted scatterplot smoothing) to model mean-variance relationship and thus does not assume a specific version of variance distribution. Therefore, it can be equally used for other high-thoughput data such as for shRNA and CRISPR library analysis. The `rosely` tools were actually initially developed for analyzing small shRNA libraries used in *in vivo* stem cell screening that were statistically challenging due to the high stochasticity of the changes.

`count_analysis` is a function to more easily call `ascertained_ttest` with the option for data pre-normalization and subsequent neutrality control for *p*-values and returns a `Result` object. The data normalization options include `log` for log-normal data but excludes 0s in the data (only positive values are used), or simiarly `log1` to transorm *x* into log(1+*x*) and to include 0s, or a power to minize the difference between mean and median (divided by range) that can be calculated by function `mean_median_norm`.

### 3. Neutrality-controlled *p*-values and their local false discovery rates (LFDR)
`neup` is a function to calculate neutralized *p*-values by raising the raw *p*-values by a power to force the larger *p*-values to approximate the neutral distribution of *p*-values (uniform distribution between 0 and 1). 

`locfdr` calculates the local false discovery rates, which measures the probability for each gene/construct to be false positive discovery. 

### 4. Presence/absence analysis
When many replicates are available and each construct does not exist in all replicates, `presence_absence_analysis` offers a Bayesian method to accurately calculate the *p*-value for the presence/absence difference in different groups for each construct. It is useful when the read counts contain exponential noises and cannot accurately measure the biological effects of the construct. It will estimate a dilution factor between groups so that different groups do not have to have the same probability for presence. For example, during continuous cell culture, the shRNAs in the library could gradually be lost. For different time points, the average probability for presence in each sample for each shRNA will not be the same. Instead of directly comparing the frequency of presences in different groups, one needs to normalize the frquencies of different shRNAs to the mean of shRNAs around the median of all shRNAs or some control shRNAs. Such normalization is not simply multiplying the ratio of frequencies, but the ratio of coverages (i.e., how many average copies of shRNAs in each sample; shRNAs would not necessarily have a 100% presence in all samples even if each sample on average contains many copies). `presence_absence_analysis` has incorporated such considerations when calculating *p*-values.

### 5. Combining the *p*-values of different constructs for the same gene
`combine_by_gene` implements both Stouffer's z-score method and Fisher's *p*-value method to combine the results of different constructs for the same gene.

### 6. KEGG pathway enrichment analysis
Module `pathwayanalysis` includes tools for Gene Ontology (GO) and KEGG pathway enrichment analysis, gene id mapping, and for drawing the enriched pathways with nodes color-coded by provided values such as z-scores and fold changes. 


## Installation<a name="Installation"></a>
### Option 1
Download and install as a normal python package:

    sudo python setup.py install

or,

    sudo python3 setup.py install

### Option 2 
If you have pip3 and git installed:

    sudo pip3 install git+https://github.com/yjiangnan/rosely.git
    
### Option 3
Use it locally and add path to search directory. This is helpful if you want to modify it.

Simply copy the folder into any directory. When using, add the directory to the path. For example, I have `rosely` under `rudolphLab/scripts` folder, the following codes imports and reload all methods in `rosely` (including code changes):

    import sys, imp
    if 'rudolphLab/scripts' not in sys.path: sys.path.insert(0, '/Users/jiangnan/Documents/rudolphLab/scripts')
    import rosely
    imp.reload(rosely)
    from rosely import *

### Dependencies
Importantly, ascertained t-test in rosely depends on a specific version of `pyloess` https://github.com/yjiangnan/pyloess for python3 support. This pyloess requires `numpy` to have `Blas` (Basic Linear Algebra Subprograms) to run, otherwise errors would occur (e.g., "undefined symbol: dswap_"). To avoid this, install Blas (e.g., [OpenBLAS](https://www.openblas.net/)) with default `make` and `make install`, then download [numpy](https://github.com/numpy/numpy) and change path to its directory, copy the file `site.cfg.example` as `site.cfg`, uncomment the lines associated with the version of Blas you installed. For `OpenBLAS`, the lines are

    [openblas]
    libraries = openblas
    library_dirs = /opt/OpenBLAS/lib
    include_dirs = /opt/OpenBLAS/include
    runtime_library_dirs = /opt/OpenBLAS/lib

Then install numpy by `python3 setup.py install`. After that, also copy the `site.cfg` to `/python_install_path/site-packages/numpy/disutils/`. 

Then, if you have git and pip3, install `pyloess` by

    sudo pip3 install git+https://github.com/yjiangnan/pyloess.git

Additionally, `rosely` depends on some python packages.

If you have pip3 installed (for python3), install them by pip3 can be done as:

    sudo pip3 install scipy pandas matplotlib sklearn psutil mygene biopython Pillow
    
### iPython notebook and Jupyter
All the examples here are iPython notebooks and it is strongly recommended that you also use iPython notebooks to keep the code, documentation and results together for good scientific practice purposes and reproducibility. To instead Jupyter, go to https://jupyter.org/.

Then, in your terminal, you can go (`cd`) to the folder of your data/notebooks, and start by running
```
jupyter notebook --script
```
to automatically open the folder in the browser and then create new notebooks or open existing notebooks.

### Install without root privileges
If you do not have root privileges, you can try using virtual environment [JuNest](https://github.com/fsquillace/junest) and build `python3`, `python-pip`, `git`, `tk` and `gcc-fortran` inside it:

    junest -f
    pacman -Sy python3 python-pip git tk gcc-fortran wget
    
However, before that, make sure that `pip3` does not pre-exist in the paths inside the junest virtual environment. Then, install python packages associated with `rosely` and run analysis inside `junest` for local or short-time interactive analysis. Or, for long runs, if you just `ssh` logged into your account on a server (outside `junest`), run
    
    nohup junest -- python3 -u script.py > log &
