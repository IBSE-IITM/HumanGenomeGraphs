## Steps to construct genome graph using VG Toolkit

**Inputs required:**

1) Reference genome: A fasta file that represent the haploid DNA sequence. This input is called in the commands below as ```Homo_sapiens_assembly38.fasta```

2) Variants Callset: A VCF file with variants in GI samples w.r.t. the above reference genome. This input is called in the commands below as ```hard_filtered.vcf```

**VCF Processing**

*Step 1: Filter only PASS variants from the hard filtered VCF file*

```
bcftools view -i "%FILTER='PASS'" hard_filtered.vcf  --threads 10 -Oz -o PASS_filtered.vcf.gz
```

*Step 2: Remove all annotation except for GT in FORMAT column*

```
bcftools annotate -x ^FORMAT/GT PASS_filtered.vcf.gz -Oz -o simplified_variants.vcf.gz
```

*Step 3: Split the MNP sites into primitive alleles*

```
vcfallelicprimitives -kg simplified_variants.vcf > processed_variants.vcf
```

*Step 4: Break multiallelic variant sites*

```
bcftools norm -m- processed_variants.vcf.gz --threads 10 -Oz -o normalized_variants.vcf.gz
```

*Step 5: Remove non-alts by selecting only SNPs, INDELs, MNPs*

```
bcftools view --types snps,indels,mnps normalized_variants.vcf.gz --threads 10 -Oz -o alt_removed_variants.vcf.gz
```


*Step 6: Process chrY variants separately*

Step 6.1: Filter chrY variants

```
bcftools view -r chrY alt_removed_variants.vcf.gz --threads 10 -Oz -o chrY_all.vcf.gz
```

Step 6.2: Extract only male samples from chrY

```
bcftools view -S male_samples.txt chrY_all.vcf.gz --min-ac 1 --threads 10 -Oz -o chrY_male.vcf.gz
```

Step 6.3: Apply MAF filter for males samples' chrY variants

```
bcftools view -q 0.05 chrY_male.vcf.gz --threads 10 | sed 's/\.\/\./0\/0/g' | sed 's/\.|\./0|0/g'  > ../final_vcf_files/chrY.vcf
```

*Step 7: Apply MAF cutoff*

```
zcat alt_removed_variants.vcf.gz | sed 's/\.\/\./0\/0/g' | sed 's/\.|\./0|0/g' | bcftools view -q 0.05  --threads 10 -Oz -o maf_filtered_variants.vcf.gz
```

*Step 8: Make one VCF file for each chromosome except chrY*

```
for i in $(seq 1 22) X ; 
do 
bcftools view -r chr$i maf_filtered_variants.vcf.gz > ../final_vcf_files/chr$i.vcf 
done
```

*Step 9: Index the individual VCF files*

```
for i in $(seq 1 22) X Y ;
do 
bgzip chr$i.vcf 
tabix -p vcf chr$i.vcf.gz
done
```


**Use VG Toolkit to create the genome graph**

*Step 10: Order VCF files as 1-22, X, Y to match reference FASTA*

```
VCFS=()
for CHROM in {1..22} X Y ; do
  VCFS+=("../final_vcf_files/chr${CHROM}.vcf.gz")
done
```


*Step 11: Create genome graph*

```
vg autoindex \
--workflow giraffe \
--ref-fasta ../reference_fasta_file/Homo_sapiens_assembly38.fasta \
--vcf "${VCFS[@]}" \
--prefix genome_graph \
--threads 10 \
--request XG
``` 