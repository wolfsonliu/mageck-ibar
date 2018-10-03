# mageck-ibar #

## Description ##

mageck-ibar is the analysis software for CRISPR/Cas9 with iBAR system applied library screening.


## Installation ##

The software is developed for Linux system using [Python 3](https://www.python.org) language.

### Dependency ###

The software takes great advantages of following Python packages.

#### Operation System ####

For the operation system, **Linux** are the applying platform. And any of the major maintaining distribution is sufficient, centos, fedora, debian, ubuntu for instance.

#### Python and Information ####

1. [Python](https://www.python.org/) (require version > 3.x)
1. [Pandas](http://pandas.pydata.org) (require version > 0.18) is an open source data structures and data analysis tools for Python.
2. [NumPy](http://www.numpy.org/) (require version > 1.10) is the fundamental Python package for scientific computing.
3. [SciPy](https://www.scipy.org) (require version > 0.17) is the python ecosystem for mathematics, science, and engineering. Pandas and NumPy are also the core packages of SciPy.
4. [RRA](https://sourceforge.net/projects/mageck/) in MAGeCk is integrated in the ibar software. 

### From Source ###

1. Download the souce code.
2. enter the directory.
3. `python3 setup.py install` or local install with `python3 setup.py install --user`

## Usage ##

```{shell}
usage: mageck-ibar [-h] -i INPUT [-b] [-n] [--col-gene COL_GENE]
                   [--col-guide COL_GUIDE] [--col-barcode COL_BARCODE]
                   -c COL_CONTROL [COL_CONTROL ...] -t COL_TREAT [COL_TREAT ...]
                   [-o OUTPREFIX] [--largerthan LARGERTHAN] [--test {norm}]
                   [--gene-test-fdr-threshold GENE_TEST_FDR_THRESHOLD]
                   [--RRApath RRAPATH] [-p {DEBUG,INFO,WARNING,ERROR}]

Analysis CRISPR/Cas9 screening data, capable of analysis data with or without
barcode integrated.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT Count table, should include <gene> <guide> <barcode> <control> <treatment>.
  -b, --with-barcode    Whether the data contain barcode.
  -n, --two-rra         Using two cycles RRA for barcode analysis.
  --col-gene COL_GENE   The column name of gene column in input file.
  --col-guide COL_GUIDE The column name of guide column in input file.
  --col-barcode COL_BARCODE The column name of barcode column in input file.
  -c COL_CONTROL [COL_CONTROL ...], --col-control COL_CONTROL [COL_CONTROL ...] The column name of control column in input file.
  -t COL_TREAT [COL_TREAT ...], --col-treat COL_TREAT [COL_TREAT ...] The column name of treatment column in input file.
  -o OUTPREFIX, --outprefix OUTPREFIX Output file prefix.
  --largerthan LARGERTHAN Normalized count should be larger than the threshold gaven, default is 10.
  --test {norm}         The test method used in analysis.
  --gene-test-fdr-threshold GENE_TEST_FDR_THRESHOLD p value threshold for alpha value of RRA in gene test (RRA -p)
  --RRApath RRAPATH     The Robust Rank Aggregation program path.
  -p {DEBUG,INFO,WARNING,ERROR}, --print-level {DEBUG,INFO,WARNING,ERROR} The information print level of the running program.
```


## Demo ##

For typical library screening data, the run time can be 5 minutes (1E6 barcodes with two replicates) or more, depending on the data size.

Sample data is provied in the ./sample directory.

### Input sample ###

161124 row data with two reference and two experiment data columns.

```{shell}
gene,guide,barcode,D0R1,D0R2,PSR1,PSR2
EMP1,AAAAAAAAGAGCCAACATGT,ATGCCA,205,225,6,6
EMP1,AAAAAAAAGAGCCAACATGT,CTAGTA,166,138,7,3
EMP1,AAAAAAAAGAGCCAACATGT,GTCGCG,32,20,1,0
EMP1,AAAAAAAAGAGCCAACATGT,TTCTTC,42,64,2,5
```

### Running code ###

Under linux ternimal:

```{shell}
mageck-ibar -i sample.csv -c D0R1 D0R2 -t PSR1 PSR2 -o ./sample_result
```

### Output file ###

#### sample_result.sgrna.txt ####

File contains the sgRNA or barcode information.

#### sample_result.phigh.txt ####

Analysis file for the up regulated sgRNAs.

#### sample_result.plow.txt ####

Analysis file for the down regulated sgRNAs.

#### sample_result.gene.high.txt ####

Result of ibar analysis for the enriched genes.

```{shell}
group_id, items_in_group, lo_value, p, FDR, goodsgrna
HPRT1, 12, 3.1175e-45, 2.5936e-07, 0.000134, 12
ITGB1, 8, 4.8659e-17, 2.5936e-07, 0.000134, 8
SRGAP2, 7, 5.8243e-15, 2.5936e-07, 0.000134, 7
```

#### sample_result.gene.low.txt ####

Result of ibar analysis for the drop-out genes.

```{shell}
group_id, items_in_group, lo_value, p, FDR, goodsgrna
THAP11, 10, 1.5596e-08, 5.1873e-07, 0.004455, 10
MAK16, 12, 4.4624e-08, 1.5562e-06, 0.004455, 12
CPSF6, 7, 4.6853e-08, 1.8156e-06, 0.004455, 7
```

## License ##

Released under GNU General Public License v3

Author: Zhiheng Liu <zhiheng.liu@pku.edu.cn>
Version: 0.1.1
