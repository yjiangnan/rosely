This folder contains three iPython Notebooks for reproducing the results of manuscript entitled 
"**High-throughput data analysis by Rosely: a genomic approach to the reproducibility crisis**".  The first two notebooks can be directly run 
to produce almost identical results in the manuscript except for some small fluctuations due to the use of random numbers. 

The first notebook [test_rosely_with_random_data](https://github.com/yjiangnan/rosely/examples/manuscript_code/test_rosely_with_random_data.ipynb)
contains all tests on the correctness of the functions/equations used by `Rosely` and the statistical properties of `Rosely` on random data.
It also shows various senarios where *p*-values get heavily biased.

The second notebook [rosely_test_yeast_48rep_and_monocyte_RNAseq.ipynb](rosely_test_yeast_48rep_and_monocyte_RNAseq.ipynb) contains tests
for the statistical properties of `Rosely` when analyzing real RNA-seq datasets, mainly highly repeated analysis of WT vs. WT yeast biological
replicates and an analysis on classical and nonclassical monocytes.


The third notebook [GEO_power_of_p-value_analysis](GEO_power_of_p-value_analysis.ipynb)  analyzing the power of *p*-values in 
43 RNA-seq studies cannot be directly run due to the dependency on many heterogenious RNA-seq 
data files downloaded from [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/gds). However, GEO Accession symbols and filenames 
are provided in the notebook so that one can individually download them if necessary. Additionally, some unzipped data files may have a 
format difficult to be directly processed by Rosely. So, some slight changes to match the loaded format shown in the notebook may be required
to successfully run the code. Some code in the analysis depends on [GEOparse](https://pypi.org/project/GEOparse/) and GEO application 
programming interface (API) so it could also fail due to possible future change/breakdown of these dependencies. 
