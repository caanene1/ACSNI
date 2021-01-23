# ACSNI
System biology information extraction for genomics.

Characterisation of molecular mechanisms underlying cellular functions allows us to develop novel therapeutic approaches. However, the current methods for pathway network inference from transcriptomic data involve single variable activity values, which leads to a significant loss of system information. We report ACSNI, a computational framework that uses artificial neural network pattern extraction and conventional machine learning approaches to infer pathway circuitry from empirical RNA-Seq evidence. ACSNI leverages prior knowledge of signalling pathways to provide context-specific interactors, thereby helping researchers better characterise the molecular basis of a biological system.

This tool is built in python3 with keras and tensorflow backend.

# To clone the source repository
git clone https://github.com/caanene1/ACSNI

# Installing the requirements
The best way to get all the dependencies for this tool is to install the "requirements.txt" included in the root directory.

On  Linux and Mac OS::
pip install -r requirements.txt

On Windows OS::
python -m pip install -U pip setuptools

# Input

Expression Matrix - The expression matrix file (expression.csv), specified by ```-i```, should be in the csv format where columns are samples and rows are genes. The expression values should be normalised (eg. TPM, CPM, RSEM). Make sure the column name of the 1st column is "gene". You can also supply gene IDs instead of gene symbols.

| gene  | Sample1 | Sample2 | Sample3 |
| --- | --- | --- | --- |
| Foxp1 | 123.2 | 274.1 | 852.6 |
| PD1  |  324.2  | 494.1  |  452.6  |
| CD8  |  523.6  | 624.1  |  252.6 |

This input should NOT be transformed in any way (e.g. log, z-scale)

Prior Matrix - The prior matrix (geneSets.csv) file, specified by ```-t```, should be in the csv format where rows are genes and column is a binary pathway membership. Where "1" means that a gene is in the pathway and "0" means that the gene is not a representative of the given pathway.
The standard prior looks like below. Make sure the column name of the 1st column is "gene". You can also supply gene IDs instead of gene symbols.

| gene  | Pathway |
| --- | --- |
| Foxp1 | 0 |
| PD1  |  0  |
| CD8  |  1  |


The tool can handle multiple pathway columns in the ```-t``` csv file as below.

| gene  | Pathway1 | Pathway2 | Pathway3 |
| --- | --- | --- | --- |
| Foxp1 | 0 | 0 | 0 |
| PD1  |  0  | 1  |  0  |
| CD8  |  1  | 0  |  1 |

Note: Each pathway above is analysed independently and the outputs have no in-built relationship.
The tool is designed to get a granular view of a single pathway.   


# Running the tool
To see input parameters use:

python ACSNI.py -help

```ACSNI.py [-h] [-m MAD] [-b BOOT] [-c ALPHA] [-p LP] [-f FULL] -i INPUT -t PRIOR [-w WEIGHT] [-s SEED]```

```optional arguments:
  -h,       --help        show this help message and exit
  -m MAD,   --mad MAD     Minimum median absolute deviance for geneSets
  -b BOOT,  --boot BOOT   Number of ensemble models to run
  -c ALPHA, --alpha ALPHA Alpha threshold to make prediction calls
                          Do not change this except you know what you are doing
  -p LP,    --lp LP       Percentage of gene set for model layers
  -p LP,    --lp LP       Percentage of gene set for model layers
  -f FULL,  --full FULL   Run tool in 1=full 0=sub (error only) mode
  -i INPUT, --input INPUT Input expression data
  -t PRIOR, --prior PRIOR Prior matrix, binary
  -w TF, --tf TF          Use weights for the genes
  -s SEED, --seed SEED    Set seed for reproducibility
```

```python ACSNI.py -i expression.csv -t GeneSets.csv```

```# Examples
Example files are included in the "example data" folder.
-i exp_derive.csv
-f biotype.csv
```


# ACSNI-derive
De-novo generation of gene sets (ACSNI-derive)

The requirement for a pre-annotated gene set is a challenge because many biological processes remain completely unknown. For example, few gene sets are available for the pathways of non-coding RNAs. To address this problem, we develop the ACSNI-derive tool for de novo creation of a gene set given a single gene. First, the expression of the gene of interest is first correlated with the expression of other genes (can also be set for a specific biotype only - eg., lncRNA, miRNAs or protein coding). Then the top most correlated targets are then selected to construct a search space and negative controls. We search this space using multiple iterations of the main model and extract genes predicted in more than 60% of the iterations in the presence of the target as the functional gene set.

# Input

Expression Matrix - The expression matrix file (exp_derive.csv), specified by ```-i```, should be in the csv format where columns are samples and rows are genes/genes_id. The expression values should be normalised (eg. TPM, CPM, RSEM). Make sure the column name of the 1st column is "gene". You can also supply gene IDs instead of gene symbols.

| gene  | Sample1 | Sample2 | Sample3 |
| --- | --- | --- | --- |
| Foxp1 | 123.2 | 274.1 | 852.6 |
| PD1  |  324.2  | 494.1  |  452.6  |
| CD8  |  523.6  | 624.1  |  252.6 |

This input should NOT be transformed in any way (e.g. log, z-scale)

PLEASE NOTE - We highly recommend to remove any constitutive/or un-desirable genes (eg. MT, RPL, Receptor genes) from the expression matrix prior to running ACSNI-derive as they usually interfere during initial prior matrix generation steps.

Biotype file (Optional) - The biotype file specified by ```-f```, should be in the csv format. This file should only be given if the generation of geneset should be based on a specific biotype (eg., "lncRNA"). If this file is provided, then the initial prior matrix will be made based on the biotype of interest. Parameter specified by ```-b``` which takes the desired biotype of interest (eg. ```-b "lncRNA"```) must be provided as well. Make sure the column names of the 1st and the 2nd column are "gene" and "biotype", respectively. You can also supply gene IDs instead of gene symbols.

| gene | biotype |
| --- | --- |
| Foxp1 | protein_coding |
| PD1  |  protein_coding  |
| MALAT1  |  lncRNA  |
| SNHG12  |  lncRNA  |
| RNU1-114P  |  snRNA  |

Correlation file (Optional) - The correlation file specified by ```-u```, should be in the csv format. This file should only be given if the user wishes to replace "some" specific genes with other genes to be used as a prior for the 1st iteration run of ACSNI. This is usually done if there is an  availability of a confident set of genes based on experimental knowledge for the tissue whose expression matrix is provided. We only recommend this if you are confident in the genes you will be testing. Otherwise don't provide this file and let the program choose the top-correlated genes from your expression matrix. Make sure the column names of the 1st and the 2nd column are "gene" and "cor", respectively. You can also supply gene IDs instead of gene symbols. Correlation values in the "cor" column of this file comes from the correlation analysis for the expression matrix provided as input with parameter ```-i```.

| gene | cor |
| --- | --- |
| Foxp1 | 0.9 |
| PD1  |  0.89  |
| MALAT1  |  0.85  |
| SNHG12  |  0.80  |
| RNU1-114P  |  0.72  |

# Running the tool

```
To set up your executable you can define a variable in your bash profile and source it:

acsni_derive=$"path/to/ACSNI-derive.py"
```

To see input parameters use:

```
python $acsni_derive -h
```

```ACSNI-derive.py [-h] -i INPUT -g GENE [-f BIO_FILE] [-b BIO_TYPE] [-m MAD] [-p LP] [-c ALPHA] [-t CT] [-z PC] [-u CORR_FILE]```

```optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input expression data (default: None)
  -g GENE, --gene GENE  Gene ID/symbol of the gene whose geneset needs to be
                        generated. (default: None)
  -f BIO_FILE, --bio_file BIO_FILE
                        Gene Bio_type table in .csv file (default: None)
  -b BIO_TYPE, --bio_type BIO_TYPE
                        Gene Bio_type of interest (default: None)
  -m MAD, --mad MAD     Minimum median absolute deviance for geneSets
                        (default: 1.2)
  -p LP, --lp LP        Percentage of gene set for model layers (default: 16)
  -c ALPHA, --alpha ALPHA
                        Alpha threshold to make prediction calls. Do not
                        change this except you know what you are doing
                        (default: 0.01)
  -t CT, --ct CT        Threshold to use for correlation (default: 0.6)
  -z PC, --pc PC        Number of prior matrix columns (default: 5)
  -u CORR_FILE, --corr_file CORR_FILE
                        File containing the top correlated genes to be used as
                        a prior for the 1st iteration of ACSNI (default: None)
```

# Example ACSNI-derive run
```
#If Biotype file is not provided
Usage: python $acsni_derive -i exp_derive.csv -m 1.2 -g "MTOR" --ct 0.6 --pc 5

#If Biotype file is provided
Usage: python $acsni_derive -f biotype.csv -b "lncRNA" -i exp_derive.csv -m 1.2 -g "SNHG12" --ct 0.6 --pc 5

```

# Main output files from ACSNI-derive are stored in the directory "main_results"
```
1) <gene>_Top_cor.csv = List of genes that were used as initial prior for ACSNI predictions.
2) <gene>_All_cor.csv = Full results of correlation for your gene of interest.
3) <gene>_fulldataset.csv = Entire dataset containing predicted and non-predicted genes.
4) <gene>_finalpredict.csv = De-novo generated genesets for your gene of interest. This file contains the list of genes that were predicted and can be used for further validation.

```
