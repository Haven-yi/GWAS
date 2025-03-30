# Genome-Wide Association Study (GWAS) 
This tutorial provides a guide to conducting GWAS using SNP and INDEL data with EMMAX. The workflow includes VCF-to-TPED conversion, phenotype data preparation, kinship matrix computation, and GWAS execution. The final output can be used for visualization through Manhattan plots to identify significant genomic associations.

Required Tools and Files  
Software  
•	EMMAX: GWAS analysis software  
•	PLINK: Genotype data processing tool  
Required Files  
•	SNP VCF File: genotype_snp.vcf  
•	INDEL VCF File: genotype_indel.vcf  
•	Phenotype Data File: phenotype_data.txt  
•	Population Structure (Q Matrix): population_structure.Q  
# 1. GWAS Analysis with SNP Data
## Step 1: Convert VCF to TPED Format (Required by EMMAX)
```
plink --vcf genotype_snp.vcf \
      --maf 0.05 --geno 0.1 \  # Filter out SNPs with MAF < 0.05 and missing rate > 10%
      --recode 12 transpose \  # Convert genotype to 12 format (TPED format)
      --output-missing-genotype 0 \  # Represent missing data as 0
      --out snp_filtered \  # Output file prefix
      --allow-extra-chr \
      --keep-allele-order
```
Output Files  
snp_filtered.log  
snp_filtered.nosex  
snp_filtered.tfam  
snp_filtered.tped  
## Step 2: Sort Phenotype Data
```
perl sort_phenotype.pl snp_filtered.tfam phenotype_data.txt > phenotype_sorted.txt
```
## Step 3: Compute Kinship Matrix
```
emmax-kin-intel64 -v -d 10 -o kinship_matrix snp_filtered
```
## Step 4: Perform GWAS Analysis
```
emmax-intel64 -v -d 10 \
              -t snp_filtered \
              -p phenotype_sorted.txt \
              -k kinship_matrix \
              -o gwas_results
```
## Step 5: Process Results for Visualization
```
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"}' snp_filtered.tped | \
paste - gwas_results.ps | \
awk 'BEGIN{print "SNP\tCHR\tBP\tP"}{if($2==$5){print $2"\t"$1"\t"$4"\t"$NF}}' > gwas_results_processed.txt
```
## Step 6: SNP GWAS with Q Matrix (Population Structure as Covariate)
```
awk '{print $1}' snp_filtered.tfam | paste - population_structure.Q | awk '{print $1" "$1" 1 "$2" "$3}' > Q_matrix_formatted.txt
emmax-intel64 -v -d 10 \
              -t snp_filtered \
              -p phenotype_sorted.txt \
              -k kinship_matrix \
              -c Q_matrix_formatted.txt \
              -o gwas_results_Q
```
# 2. GWAS Analysis with INDEL Data
## Step 1: Convert VCF to TPED Format
```
plink --vcf genotype_indel.vcf \
      --maf 0.01 --geno 0.5 \
      --recode 12 transpose \
      --out indel_filtered \
      --allow-extra-chr \
      --keep-allele-order
```
## Step 2: Sort Phenotype Data
```
perl sort_phenotype.pl indel_filtered.tfam phenotype_data.txt > phenotype_sorted.txt
```
## Step 3: Compute Kinship Matrix
```
emmax-kin-intel64 -v -d 10 -o kinship_matrix_indel indel_filtered
```
## Step 4: Run GWAS Analysis
```
emmax-intel64 -v -d 10 \
              -t indel_filtered \
              -p phenotype_sorted.txt \
              -k kinship_matrix_indel \
              -o gwas_results_indel
```
## Step 5: Process Results for Visualization
```
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"}' indel_filtered.tped | \
paste - gwas_results_indel.ps | \
awk 'BEGIN{print "SNP\tCHR\tBP\tP"}{if($2==$5){print $2"\t"$1"\t"$4"\t"$NF}}' > gwas_results_indel_processed.txt
```

