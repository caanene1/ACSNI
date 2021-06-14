# OVERVIEW

Determining context-specific circuit of biological pathways is a fundamental goal of molecular cell biology.
ACSNI combines prior knowledge of biological processes (gene set) with a deep neural network to decompose gene 
expression profiles (GEP) into pathway activities and identify unknown pathway components, see Anene et al., 2021.

# Required inputs

1. Gene expression matrix with genes in rows and samples in columns (format .csv).
2. Gene set membership file representing prior knowledge of gene functions (format .csv)
   For a single gene analysis, the second input is a gene name provided at run time.
3. Optional- weight file (as integer values) for the genes in the second input. 


# CASES

For this tutorial, we will use the cases reported in the manuscript (Anene et al., 2021) (herein referred to as MTOR, ATF2 and HOTAIRM1 cases) to demonstrate how to set up, run ACSNI, interpret the results and navigate the extended database. 

# INSTALLATION

You can install ACSNI with the PIP command, which automates installing the required packages (e.g. Tensorflow). 

Ensure you have python version 3.8 installed before running the code below; if not, see https://www.python.org/downloads/


```python
!pip3 install ACSNI
# or 
!pip install ACSNI
```

The above should install the latest version of ACSNI. 

In addition to specifying a version during installation (pythonic way), 
you can also install directly from the .wheel file provided at https://github.com/caanene1/ACSNI or
compile the code yourself using setup.py.

ACSNI has three entry commands, including:
- ACSNI-run : multiple genes prior
- ACSNI-derive : single gene prior
- ACSNI-get : phenotype linking

Run the code with option -h to check installation and parameters. 
Further, you can use ACSNI functions can be used in regular python imports and calls.

Below are the arguments for the ACSNI-run entrance.


```python
!ACSNI-run -h
```

    usage: ACSNI-run [-h] [-m MAD] [-b BOOT] [-c ALPHA] [-p LP] [-f FULL] -i INPUT
                     -t PRIOR [-w WEIGHT] [-s SEED]
    
    System biology information extraction for genomics.
    
    optional arguments:
      -h, --help            show this help message and exit
      -m MAD, --mad MAD     Minimum median absolute deviance for geneSets
      -b BOOT, --boot BOOT  Number of ensemble models to run
      -c ALPHA, --alpha ALPHA
                            Alpha threshold to make prediction calls
      -p LP, --lp LP        Dimension of the pathway layer. It is also half of the
                            subprocess,set to 0 or default for automatic
                            estimation
      -f FULL, --full FULL  Run tool in 1=full 0=sub (error only) mode
      -i INPUT, --input INPUT
                            Input expression data (.csv)
      -t PRIOR, --prior PRIOR
                            Prior matrix, binary
      -w WEIGHT, --weight WEIGHT
                            Use weights for the genes
      -s SEED, --seed SEED  Set seed for reproducibility


Except for the two required arguments -i and -t, the rest of the arguments have well-tested defaults. 
You can tune these parameters to your specific needs. Caution!!

# mTOR case
The first case infers the extended mTOR signalling network in clear cell renal cell carcinoma (ccRCC).

Input Files (included): 

    1. TCGA_.csv - gene expression matrix (source-TCGA)
    2. mTOR.csv - mTOR gene set from Pathway interaction database
    3. sample_info.csv - sample phenotype for the first input

Here, we set the:

    -f to 1 for the full run and output
    -p to 0 for automatic estimation of layers
    -b to 5 for five models and 
    -m to 2.5 for minimum absolute deviation


```python
!ACSNI-run -i TCGA_.csv -t mTOR.csv -f 1 -p 0 -b 5 -m 2.5
```

    Running for PID_MTOR_4PATHWAY
    Results will be saved to TCGA_PID_MTOR_4PATHWAY-704AHXJ
    Geneset with 67 genes in the expression
    146 samples
    2021-04-02 20:13:24.916282: I tensorflow/compiler/jit/xla_cpu_device.cc:41] Not creating XLA devices, tf_xla_enable_xla_devices not set
    2021-04-02 20:13:24.919062: I tensorflow/core/platform/cpu_feature_guard.cc:142] This TensorFlow binary is optimized with oneAPI Deep Neural Network Library (oneDNN) to use the following CPU instructions in performance-critical operations:  AVX2 FMA
    To enable them in other operations, rebuild TensorFlow with the appropriate compiler flags.
    2021-04-02 20:13:25.041742: I tensorflow/compiler/mlir/mlir_graph_optimization_pass.cc:116] None of the MLIR optimization passes are enabled (registered 2)
    The optimal number of dimension is 11
    Geneset with 67 genes in the expression
    146 samples
    Geneset with 67 genes in the expression
    146 samples
    Geneset with 67 genes in the expression
    146 samples
    Geneset with 67 genes in the expression
    146 samples


The command will print the run progress and location of the output.
As shown above, the results are in the folder "TCGA_PID_MTOR_4PATHWAY-704AHXJ".
We can visualise the contents of the output folder below.


```python
!ls TCGA_PID_MTOR_4PATHWAY-704AHXJ
```

    NULL_TCGA_.csv    Network_TCGA_.csv dbsTCGA_.ptl


The output of ACSNI-run can be identified by the prefix, specifically:

**NULL** - randomly shuffled expression matrix derived from -i. You can use it to check randomness in the predictions. For this, rerun the model with the argument -i set to this file.

**dbs** - is the database of intermediate files and run details.

**Network** - the inferred network and components. Most users only need this output.

The network file has three columns including:

    1. name - gene name, 
    2. sub - name of the inferred sub-network 
    3. Direction - strength of the interaction in the sub-network


```python
!head TCGA_PID_MTOR_4PATHWAY-704AHXJ/Network_TCGA_.csv
```

    name,sub,Direction
    ANAPC13,AE_SE8IA2_0,0.0067188637331128
    ANKRD46,AE_SE8IA2_0,0.0033135202247649
    ATP9A,AE_SE8IA2_0,0.0136904744431376
    BPGM,AE_SE8IA2_0,0.0053974902257323
    C18orf32,AE_SE8IA2_0,0.0077944286167621
    C6orf226,AE_SE8IA2_0,0.0003798712568823
    C8orf80,AE_SE8IA2_0,-0.0013907115207985
    CCDC56,AE_SE8IA2_0,0.0074006947688758
    COL24A1,AE_SE8IA2_0,0.0001676924148341


# Advanced usage
You can further inspect the dbs output of ACSNI-run (as well as ACSNI-derive) using the python pickle package. Internally, the database is a python class with methods and member variables.

This database also contains predictions made through linear decomposition approaches; see Anene et al., 2021.


```python
# Import the os and pickle packages
import os
import pickle

# Load the database and visualise the information
database = pickle.load(open("TCGA_PID_MTOR_4PATHWAY-704AHXJ/dbsTCGA_.ptl", "rb"))
database
```




    
    ACSNI result with 14 modules over 5 bootstraps.




```python
# Extract run information
database.get_run_info()

# Extract predicted network and save to file
output = database.get_p()
output.to_csv("TCGA_PID_MTOR_4PATHWAY-704AHXJ/results.csv")

# This file has extended output, including results from decomposing with PCA, NMF and median
output.head(5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AE_SE8IA2_0</th>
      <th>AE_SE8IA2_1</th>
      <th>AE_SE8IA2_2</th>
      <th>AE_SE8IA2_3</th>
      <th>AE_SE8IA2_4</th>
      <th>AE_SE8IA2_5</th>
      <th>AE_SE8IA2_6</th>
      <th>AE_SE8IA2_7</th>
      <th>AE_SE8IA2_8</th>
      <th>AE_SE8IA2_9</th>
      <th>...</th>
      <th>AE_QD6XXF_8</th>
      <th>AE_QD6XXF_9</th>
      <th>AE_QD6XXF_10</th>
      <th>PCA_QD6XXF_0</th>
      <th>NMF_QD6XXF_0</th>
      <th>NMF_QD6XXF_1</th>
      <th>MEDIAN_QD6XXF_0</th>
      <th>Predicted</th>
      <th>Sum_stat</th>
      <th>Boot_Count</th>
    </tr>
    <tr>
      <th>name</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A1BG</th>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>...</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
    </tr>
    <tr>
      <th>A1CF</th>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>...</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>0</td>
      <td>1</td>
      <td>5</td>
    </tr>
    <tr>
      <th>A2LD1</th>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>...</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>0</td>
      <td>1</td>
      <td>5</td>
    </tr>
    <tr>
      <th>A2M</th>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>...</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
    </tr>
    <tr>
      <th>A4GALT</th>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>...</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>B</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 82 columns</p>
</div>



Above the column prefix allows for additional assessment "AE_" (autoencoder) against linear methods decomposition methods, including
> NMF_ non-negative factorisation
> PCA_ principle component analysis 
> MEDIAN_ simple median expression 

# Linking to phenotype
Further, you can link the predicted sub-processes to clinical or biological information. 
The included "sample_info.csv" contains a grouping variable 1.normal and 2.tumour.

We can check the file


```python
!head sample_info.csv
```

Then, we can run the ACSNI-get command with the database and the sample information file.


```python
!ACSNI-get -r TCGA_PID_MTOR_4PATHWAY-704AHXJ/dbsTCGA_.ptl -v sample_info.csv -c character
!head group_to\ subprocess\ associations.csv
```

    Statistics of variations in subprocesses explained by group
    q25 0.5068493150684932 
     q75 0.5856164383561644 
     mean 0.5620174346201743 
     std 0.09858342604975479
    sub,0,Association
    AE_SE8IA2_0,0.5342465753424658,Weak
    AE_SE8IA2_1,0.5068493150684932,Weak
    AE_SE8IA2_2,0.4863013698630137,Weak
    AE_SE8IA2_3,0.6712328767123288,Strong
    AE_SE8IA2_4,0.589041095890411,Strong
    AE_SE8IA2_5,0.5068493150684932,Weak
    AE_SE8IA2_6,0.7054794520547946,Strong
    AE_SE8IA2_7,0.5958904109589042,Strong
    AE_SE8IA2_8,0.4863013698630137,Weak


 

# ATF2 case
The second case investigates ATF2-dependent bzip transcriptional output in healthy artery aorta.

Input Files (included): 

    1. AA_.csv - gene expression matrix (source-GTEX)
    2. ATF2.csv - ATF2 gene set from Pathway interaction database

Here, we set the:

    -f to 1 for a full run and output
    -p to 0 for automatic estimation of layers
    -b to 10 for ten models and 
    -m to 1.2 for minimum absolute deviation


```python
!ACSNI-run -i AA_.csv -t ATF2.csv -f 1 -p 0 -b 10 -m 1.2
```

    Running for PID_ATF2_PATHWAY
    Results will be saved to AA_PID_ATF2_PATHWAY-51GJJUE
    Geneset with 37 genes in the expression
    432 samples
    2021-04-02 20:34:22.983361: I tensorflow/compiler/jit/xla_cpu_device.cc:41] Not creating XLA devices, tf_xla_enable_xla_devices not set
    2021-04-02 20:34:22.983526: I tensorflow/core/platform/cpu_feature_guard.cc:142] This TensorFlow binary is optimized with oneAPI Deep Neural Network Library (oneDNN) to use the following CPU instructions in performance-critical operations:  AVX2 FMA
    To enable them in other operations, rebuild TensorFlow with the appropriate compiler flags.
    2021-04-02 20:34:23.045616: I tensorflow/compiler/mlir/mlir_graph_optimization_pass.cc:116] None of the MLIR optimization passes are enabled (registered 2)
    The optimal number of dimension is 5
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples
    Geneset with 37 genes in the expression
    432 samples


The generated outputs can be interpreted in a similar manner to the outputs of the first case above.

# HOTAIRM1 case
Third case explores the regulatory network of the lncRNA HOTAIRM1 in healthy kidney. 
It is a case of a single gene, thus we need the ACSNI-derive command.

Input Files (included):

    1. KID_.csv - gene expression matrix (source-GTEX)
    2. ENSG00000233429 - gene ID for HOTAIRM1 in the first input
    3. biotype.csv - gene biotype information for the first input (grch38)
    4. exclude.csv - exclude biotype file for the third input 
    
Before runing the command, we can check the arguments as before.


```python
!ACSNI-derive -h
```

    usage: ACSNI-derive [-h] -i INPUT -g GENE [-f BIO_FILE] [-b BIO_TYPE] [-m MAD]
                        [-p LP] [-c ALPHA] [-ex EXCLUDE] [-t CT] [-z PC]
                        [-u CORR_FILE]
    
    De-Novo generation of gene sets
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Input expression data (.csv)
      -g GENE, --gene GENE  Gene ID/symbol to analyse
      -f BIO_FILE, --bio_file BIO_FILE
                            Gene Bio_type table (.csv)
      -b BIO_TYPE, --bio_type BIO_TYPE
                            Gene Bio_type of interest
      -m MAD, --mad MAD     Minimum median absolute deviation
      -p LP, --lp LP        Percentage of gene_set for model layers
      -c ALPHA, --alpha ALPHA
                            Alpha threshold to make prediction calls
      -ex EXCLUDE, --exclude EXCLUDE
                            Name of bio_types to exclude in csv format
      -t CT, --ct CT        Threshold to use for correlation
      -z PC, --pc PC        Number of prior matrix columns
      -u CORR_FILE, --corr_file CORR_FILE
                            Pre-computed top correlated genes (.csv)


Here, we set the:
                 
    -b to "lncRNA" to restrict the de novo prior to only lncRNAs 
         (optional see Anene et al., 2021). In most cases, you do not need to specify this.
    -t to 0.8 for correlation filtering
    -z to 5 for five models and 


Run the command to get the network of the gene of interest using automatic estimation method (**-p 0**):


```python
!ACSNI-derive -f biotype.csv -b "lncRNA" -i KID_.csv -m 1.2 -g "ENSG00000233429" -t 0.80 -z 5 -ex exclude.csv
```

or using a user defined fixed estimation method (example **-p 16**):


```python
ACSNI-derive -f biotype.csv -b "lncRNA" -i KID_.csv -m 1.2 -g "ENSG00000233429" -t 0.80 -z 5 -ex exclude.csv -p 16
```

Both will return three outputs, which can be identifed by their prefix, including;
**Predicted** which is a list of predicted genes that are under the network of the gene of interest (-g). Most users of this command only need this file.


The remaining outputs are **NULL** data for randomness assessment and "dbs" database for intermediate files checks and mining, as described above. Please, note that the **dbs** file here cannot be used with ACSNI-get command as it is not meaningful. Here, you can correlate the expression of the single gene to your phenotype directly to obtain similar results as ACSNI-get.

# Conclusion
The tutorial shows how to use the different ACSNI commands. The predicted genes are the components of the analysed pathways or network; see the original manuscript. Ultimately, you should apply orthogonal validation to refine the predictions; please see the original manuscript for an extended discussion on such approaches.

The included R-scripts further demonstrates some downstream analysis and other approaches are available depending on your intended use case (hypothesis generation or conclusive insight).
