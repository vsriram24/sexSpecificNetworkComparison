#!/bin/bash
#$ -N vivek_male_rg
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=250:00:00
#$ -t 1-500

# Function: Run Linkage Disequilibrium Score Regression on all pairs of phenotypes
#           to calculate genetic correlations between each pair of disease nodes.

# Get an input list of phenotypes between which to calculate correlations
FILENAME="phenotypesToCalculateCorrelations.txt"
# Get phenotypes from the input list
PHE=$(cat $FILENAME)
# Initialize a list of phenotype pairs to be processed
declare -a phenotypeComboArray=()

# Nested for loop to process all pairs of phenotypes
for PHECODE1 in $PHE
do
        for PHECODE2 in $PHE
        do
                phecodeCombo1=${PHECODE1}_${PHECODE2}
                phecodeCombo2=${PHECODE2}_${PHECODE1}

                # If we haven't seen this pair of traits before, then we add it to our list and run LDSC
                if [[ ! " ${phenotypeComboArray[*]} " =~ " ${phecodeCombo1} " ]]; then
                        if [[ ! " ${phenotypeComboArray[*]} " =~ " ${phecodeCombo2} " ]]; then
                                phenotypeComboArray=(${phenotypeComboArray[@]} ${phecodeCombo1} ${phecodeCombo2})

                                # Run LDSC to calculate genetic correlation between the two phenotypes
                                source activate ldsc

                                ldsc/ldsc.py \
                                --rg rawNealeLabVCFFiles/male/${PHECODE1}.gwas.imputed_v3.male.vcf,rawNealeLabVCFFiles/male/${PHECODE1}.gwas.imputed_v3.male.vcf \
                                --ref-ld-chr /cbica/projects/kimlab/software_dir/ldsc/eur_w_ld_chr/ \
                                --w-ld-chr /cbica/projects/kimlab/software_dir/ldsc/eur_w_ld_chr/ \
                                --out rgResults/male/result_phe1_${PHECODE1}_phe2_${PHECODE2}
                        fi
                fi
             
        done
done