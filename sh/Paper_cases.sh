# Case 1: Figure 2A-G
ACSNI-run -i TCGA_.csv -t mTOR.csv -f 1 -p 10 -b 5 -m 2.5
ACSNI-run -i TCGA_.csv -t mTOR.csv -f 1 -p 0 -b 5 -m 2.5
ACSNI-get -r dbsTCGA_.ptl -v sample_info_kidney_norm_cancer.csv -c character

# Case 2: Figure 2H-K
ACSNI-run -i AA_.csv -t ATF2.csv -f 1 -p 5 -b 10 -m 1.2
ACSNI-run -i AA_.csv -t ATF2.csv -f 1 -p 0 -b 10 -m 1.2

# Case 3: Figure 2L-N
ACSNI-derive -f grch38_biotype2.csv -b "lncRNA" -i GTEx_kidney_exp_norerplmt_eIDs.csv -m 1.2 -g "ENSG00000233429" --ct 0.80 --pc 5 --ex exclude.csv

ACSNI-derive -f biotype.csv -b "lncRNA" -i expression.csv -m 1.2 -g "ENSG00000233429" --ct 0.80 --pc 5 --ex exclude.csv
