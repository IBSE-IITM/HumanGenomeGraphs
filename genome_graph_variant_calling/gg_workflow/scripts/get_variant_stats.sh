# Script to extract variant counts and QC metrics from vcf files created by "rule variant_filtering" in the gg_workflow
# Usage: bash get_variant_stats.sh [vcf files of samples]
#
# Author: Venkatesh Kamaraj
# Contact: ic11570@imail.iitm.ac.in

# Command below is for error handling in bash. Stops execution of the script once an error is faced.
set -euo pipefail

for sample_vcf in "$@"
do
    total_variants=$(bcftools stats $sample_vcf | sed -n '24p' | awk {'print $6'});
    SNPs=$(bcftools stats $sample_vcf | sed -n '26p' | awk {'print $6'});
    Indels=$(bcftools stats $sample_vcf | sed -n '28p' | awk {'print $6'});
    ts=$(bcftools stats $sample_vcf | sed -n '34p' | awk {'print $3'});
    tv=$(bcftools stats $sample_vcf | sed -n '34p' | awk {'print $4'});
    het_indels=$(bcftools view --genotype het --types indels $sample_vcf | bcftools stats | sed -n '24p' | awk {'print $6'});
    hom_indels=$(bcftools view --genotype hom --types indels $sample_vcf | bcftools stats | sed -n '24p' | awk {'print $6'});
    het_snps=$(bcftools view --genotype het --types snps $sample_vcf | bcftools stats | sed -n '24p' | awk {'print $6'});
    hom_snps=$(bcftools view --genotype hom --types snps $sample_vcf | bcftools stats | sed -n '24p' | awk {'print $6'});

    sample_name=$(echo $sample_vcf | \
    awk '{gsub("_filtered_variants.recode.vcf", "")}1' | \
    awk '{gsub("data/variant_calling/filtered_variants/", "")}1')

    echo "$sample_name, Graph, $total_variants, $SNPs, $Indels, $ts, $tv, $het_indels, $hom_indels, $het_snps, $hom_snps" >> data/variant_calling/summary_gg_pipeline_variants.csv
done
