# rosely
### Robust and sensitive tools for high-throughput data analysis

`rosely` mainly implements the following tools:
### 1. Classes for reading and normalizing data and for results
`ReadCounts` is a class for reading excel and text files containing read counts for each gene/constructs in each sample. It implements data normalization (to bring different samples to the same level) by either a specific control gene/construct, or the median of all data, or a quantile (1-9) of all data around the median, or a list of control genes/constructs.

`Result` is a class to wrap the results of raw *p*-values, *z*-scores, their neutralized values, local false discovery rates (and a member function for recalculation from raw *p*-values with different parameters) and the changes of normalized values. It has a member `results` of type of `pandas DataFrame` to utilize `pandas` tools for further processing.

### 2. Ascertained *t*-test
`ascertained_ttest` calculates *p*-values (with high statistical significance) of the changes of genes/constructs in different experimental groups using the similarity of genes/constructs in the same experimental group to reduce the uncertainty of the variance of each gene/construct in a rigorous empirical Bayesian approach. It has careful consideration of NOT overestimating the degrees of freedom and statistical significance. It outperforms mainstream RNA-seq analysis tools edgeR and DESeq2 for controlling false discovery rate in comparison of WT vs. WT from 48 yeast replicates. It uses LOWESS regression (locally weighted scatterplot smoothing) to model mean-variance relationship and thus does not assume a specific version of variance distribution. Therefore, it can be equally used for other high-thoughput data such as for shRNA and CRISPR library analysis. The `rosely` tools were actually initially developed for analyzing small shRNA libraries used in *in vivo* stem cell screening that were statistically challenging due to the high stochasticity of the changes.

`count_analysis` is a function to more easily call `ascertained_ttest` with the option for data pre-normalization and subsequent neutrality control for *p*-values and returns a `Result` object. The data normalization options include `log` for log-normal data but excludes 0s in the data (only positive values are used), or simiarly `log1` to transorm *x* into log(1+*x*) and to include 0s, or a power to minize the difference between mean and median (divided by range) that can be calculated by function `mean_median_norm`.

### 3. Neutrality-controlled *p*-values and their local false discovery rates (LFDR)
`neup` is a function to calculate neutralized *p*-values by raising the raw *p*-values by a power to force the larger *p*-values to approximate the neutral distribution of *p*-values (uniform distribution between 0 and 1). 

`locfdr` calculates the local false discovery rates, which measures the probability for each gene/construct to be false positive discovery. Note: conventional false discovery rate (FDR or *q*-value) is a measurement for all the genes/constructs whose *p*-values are small; it says nothing about the probability of false discovery for the gene/construct itself.

### 4. Presence/absence analysis
When many replicates are available and each construct does not exist in all replicates, `presence_absence_analysis` offers a Bayesian method to accurately calculate the *p*-value for the presence/absence difference in different groups for each construct. It is useful when the read counts contain exponential noises and cannot accurately measure the biological effects of the construct. It will estimate a dilution factor between groups so that different groups do not have to have the same probability for presence. For example, during continuous cell culture, the shRNAs in the library could gradually be lost. For different time points, the average probability for presence in each sample for each shRNA will not be the same. Instead of directly comparing the frequency of presences in different groups, one needs to normalize the frquencies of different shRNAs to the mean of shRNAs around the median of all shRNAs or some control shRNAs. Such normalization is not simply multiplying the ratio of frequencies, but the ratio of coverages (i.e., how many average copies of shRNAs in each sample; shRNAs would not necessarily have a 100% presence in all samples even if each sample on average contains many copies). `presence_absence_analysis` has incorporated such considerations when calculating *p*-values.

### 5. Combining the *p*-values of different constructs for the same gene
`combine_by_gene` implements both Stouffer's z-score method and Fisher's *p*-value method to combine the results of different constructs for the same gene.


## Installation
### Option 1
Install as a normal python package:

    sudo python setup.py install

or,

    sudo python3 setup.py install
    
### Option 2
Use it locally and add path to search directory. This is helpful if you want to modify it.

Simply copy the folder into any directory. When using, add the directory to the path. For example, I have `rosely` under `rudolphLab/scripts` folder, the following codes imports and reload all methods in `rosely` (including code changes):

    import sys, imp
    if 'rudolphLab/scripts' not in sys.path: sys.path.insert(0, '/Users/jiangnan/Documents/rudolphLab/scripts')
    import rosely
    imp.reload(rosely)
    from rosely import *

### Dependencies
Importantly, ascertained t-test in rosely depends on a specific version of pyloess https://github.com/jcrotinger/pyloess for python3 support. You need to download and install it to use ascertained t-test.

Additionally, `rosely` depends on the following python packages:

`numpy scipy pandas matplotlib sklearn`

If you have pip3 installed (for python3), install them by pip3 can be done as:

    sudo pip3 install numpy scipy pandas matplotlib sklearn
