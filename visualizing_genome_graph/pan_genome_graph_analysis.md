> conda activate gg_workflow

**Download VCF files**

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/phase3.chr21.GRCh38.GT.crossmap.vcf.gz


wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/phase3.chr{1..21}.GRCh38.GT.crossmap.vcf.gz.tbi


# VCF processing

**Rename files for convenience**
mv ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi 1KGP_chrY.vcf.gz.tbi
mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz 1KGP_chrX.vcf.gz

for CHROM in {1..22} ; 
do
mv ALL.chr${CHROM}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi 1KGP_chr${CHROM}.vcf.gz.tbi
done


**Add chr to all chromosomes to match with reference genome fasta**
for vcf_file in $(ls *.vcf.gz)
do
vcf_name=$(echo $vcf_file | awk '{gsub(".GRCh38.GT.crossmap.vcf.gz", "")}1')
zcat $vcf_file | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}'  > with_chr_${vcf_name}.vcf
bgzip with_chr_${vcf_name}.vcf
tabix -p vcf with_chr_${vcf_name}.vcf.gz
done




**Apply MAF filter > 5%**
*chrY has h=genomtypes for only males*
for vcf_file in $(ls ../1KGP_phase3_hg38_full_vcf_files/phase3.chr*.vcf.gz)
do
vcf_name=$(echo $vcf_file | awk '{gsub("../1KGP_phase3_hg38_full_vcf_files/phase3.", "")}1' | awk '{gsub(".GT.crossmap.vcf.gz", "")}1')
echo $vcf_name
bcftools view -q 0.05 $vcf_file --threads 20 -Oz -o ${vcf_name}_maf_filtered.vcf.gz
done

for i in $(seq 1 22) X Y 
do 
echo $i
tabix -p vcf phase3.chr$i.GRCh38.GT.crossmap.vcf.gz
done


for i in $(seq 1 22) X Y 
do 
echo $i
mv 1KGP_complete_chr$i.vg 1KGP_maf_filtered_chr$i.vg
done


# VG manual construction

**VG files construction**
> cd vg_construct


> ((seq 1 22; echo X; echo Y) | \
parallel -j 24 \
"time \
vg construct --alt-paths \
--region-is-chrom --region chr{} \
--reference /cn4/data4/venkatesh_data4/1KGP_genome_graphs/reference_genome/Homo_sapiens_assembly38.fasta \
--vcf /cn4/data4/venkatesh_data4/1KGP_genome_graphs/1KGP_vcf_files/liftover_files/1KGP_phase3_hg38_maf_filtered_vcf_files/chr{}.GRCh38_maf_filtered.vcf.gz \
--threads 1 --flat-alts --handle-sv > 1KGP_maf_filtered_chr{}.vg") 2>> vg_construct.log

**Node ID coordination**
> vg ids -j -m mapping $(for i in $(seq 1 22; echo X; echo Y); do echo 1KGP_maf_filtered_chr$i.vg; done)

> cp mapping mapping.backup

**GFA creation**
> (seq 1 22; echo X; echo Y) | \
parallel -j 24 \
"time \
vg view -v /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/vg_files/1KGP_complete_chr{}.vg \
-g > /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr{}.gfa"

**Indexing with XG with alt paths visible**
> vg index -L -x 1KGP_maf_filtered_graph.xg \
$(for i in $(seq 22; echo X; echo Y); do echo /cn4/data4/venkatesh_data4/1KGP_genome_graphs/maf_filtered_1KGP_graph/vg_manual_construct/vg_files/1KGP_maf_filtered_chr$i.vg; done) \
--threads 20

## AutoIndexing

**VCF order must be 1-22, X, Y to match FASTA**
> VCFS=()
> for CHROM in {1..22} X Y ; do
  VCFS+=("/cn4/data4/venkatesh_data4/1KGP_genome_graphs/1KGP_vcf_files/liftover_files/1KGP_phase3_hg38_maf_filtered_vcf_files/chr${CHROM}.GRCh38_maf_filtered.vcf.gz")
done

**Create graph and index it**
> vg autoindex \
--workflow giraffe \
--ref-fasta /cn4/data4/venkatesh_data4/1KGP_genome_graphs/reference_genome/Homo_sapiens_assembly38.fasta \
--vcf "${VCFS[@]}" \
--prefix gi_200_graph \
--threads 10 \
--request XG \
--request GBWT

*non-alt sites are being skipped*


# Variant calling with genome graphs

**Run the genome graph pipeline for GIAB samples**
*Dry run*
snakemake -npr

*Full run*
snakemake --jobs -p -r --use-conda --cluster-config config/cluster.json --cluster "qsub -S /bin/bash -l ncpus={cluster.ppn}"


# Visualizing genome graph

*Activate the conda environment*
> conda activate viz_vcf



**Create 10MB subgraphs GFA**
*size of chrs in hg38 from ucsc*
> https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes

*Get sizes of chr 1-22,x,y only*
> head hg38.chrom.sizes -n 24 > hg38_subset.chrom.sizes

*make a bed file for a preffered window size*
> bedtools makewindows -g hg38_subset.chrom.sizes -w 10000000 > 10M_windows.bed

> cd /cn4/data4/venkatesh_data4/1KGP_genome_graphs/maf_filtered_1KGP_graph/vg_manual_construct/gfa_files/10MB_subgraphs/

> qsub create_gfa.sh


**Find complexity of genome graph in each bin**
> python genome_graph_complexity_finder.py


**Make circos plot of genome graph**
*Edit the genome_graph_complexity_finder.py output to make compatibility with circos*
> awk '{if($0 !~ /^#/) print "hs"$0; else print $0}' bin_level_complete_1kgp_graph_complexity.txt > bin_level_complete_1kgp_graph_complexity.dat

> awk '{if($0 !~ /^#/) print "hs"$0; else print $0}' bin_level_maf_1kgp_graph_highly_complexity.txt > bin_level_maf_1kgp_graph_highly_complexity.dat

*Get max count from file and edit circos config files*
> sort -k4nr bin_level_maf_1kgp_graph_complexity.dat | head -n 1 
> sort -k4nr bin_level_maf_1kgp_graph_highly_complexity.dat | head -n 1 

> circos  -outputdir ./   -outputfile  graph_complexity_circos_plot -conf ./config_files/graph_complexity_circos.conf


**Plots for hypervariable regions**
> RScript plot_hypervariable_distribution.R


**plotting paths in the most complex region of the graph**
*Get node ids for the complex nodes from chr1. Most complex nodes are identified from plot_hypervariable_distribution.R*
> cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/10MB_subgraphs/gfa_files/chr1_200000001_210000001.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 12'

*19298453 and 19298461 have degree 12*

> vg find -n 19298457 -c 100 -x /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/vg_files/1KGP_complete_chr1.vg | vg view -v - -d | dot -Tsvg -o complete_1kgp_gg_most_complex_region.svg


> vg find -n 19298457 -c 7 -x /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/vg_files/1KGP_complete_chr1.vg | vg view -v - -d | dot -Tsvg -o complete_1kgp_gg_most_complex_region_trimmed.svg

vg find -n 291240874 -c 10 -x /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/vg_files/1KGP_complete_chrX.vg | vg view -v - -d | dot -Tsvg -o complete_1kgp_x_chr_hypvar.svg


**Find invariable regions for each chromosome**
> python genome_graph_invariability_finder.py


**Plots for invariable regions**
> RScript plot_invariant_distribution.R



# HAP.PY analysis of GIAB variants
****
> for i in $(ls *vcf);
do
bcftools sort $i -T ./tmp -Oz -o ${i}.gz
tabix -p vcf ${i}.gz
done




**Filter pass variants and get variant stats for linear pipeline**


> bash get_variant_stats.sh /cn4/data4/venkatesh_data4/1KGP_genome_graphs/giab_happy_analysis/linear_vcf_files/vcf_files_from_gatk/*.vcf.gz


*Single file comparison*
hap.py /cn4/data4/venkatesh_data4/1KGP_genome_graphs/giab_happy_analysis/GIAB_high_confidence_variants/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /cn4/data4/venkatesh_data4/1KGP_genome_graphs/giab_happy_analysis/complete_1kgp_vcf_files/vcf_files_from_gg/HG001_filtered_variants.recode.vcf.gz \
-f /cn4/data4/venkatesh_data4/1KGP_genome_graphs/giab_happy_analysis/GIAB_high_confidence_variants/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
-r /cn4/data4/venkatesh_data4/1KGP_genome_graphs/reference_genome/Homo_sapiens_assembly38.fasta \
-o complete_1kgp_gg_HG001 \
--threads 20


*All files in a loop*
for i in {1..7};
do
echo "HG00${i}: hap.py analysis started"
hap.py /cn4/data4/venkatesh_data4/1KGP_genome_graphs/giab_happy_analysis/GIAB_high_confidence_variants/HG00${i}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /cn4/data4/venkatesh_data4/1KGP_genome_graphs/giab_happy_analysis/maf_filtered_1kgp_vcf_files/vcf_files_from_gg/HG00${i}_filtered_variants.recode.vcf.gz \
-f /cn4/data4/venkatesh_data4/1KGP_genome_graphs/giab_happy_analysis/GIAB_high_confidence_variants/HG00${i}_GRCh38_1_22_v4.2.1_benchmark.bed \
-r /cn4/data4/venkatesh_data4/1KGP_genome_graphs/reference_genome/Homo_sapiens_assembly38.fasta \
-o HG00${i}/HG00${i} \
--threads 20
echo "HG00${i} process started"
done


**Collate summaries for all (3 pipelines*7 samples) into a single file**
echo "File,$(cat maf_filtered_1kgp_vcf_files/happy_analysis/HG001/HG001.summary.csv | head -n 1)" > happy_results_merged_summary.csv


for summary_file in $(ls */happy_analysis/HG00*/*summary.csv);
do
echo "${summary_file},  $(cat ${summary_file} | head -n 2 | tail -n 1)" >> happy_results_merged_summary.csv
echo "${summary_file},  $(cat ${summary_file} | head -n 3 | tail -n 1)" >> happy_results_merged_summary.csv
echo "${summary_file},  $(cat ${summary_file} | head -n 4 | tail -n 1)" >> happy_results_merged_summary.csv
echo "${summary_file},  $(cat ${summary_file} | head -n 5 | tail -n 1)" >> happy_results_merged_summary.csv
done



























# Saturation analysis for variant count

**Filter variants for different stratified sample sets**
> for i in $(seq 1 22) X Y 
do 
echo $i
bcftools view \
--samples-file ../../startified_1252_samples.txt \
/cn4/data4/venkatesh_data4/1KGP_genome_graphs/1KGP_vcf_files/liftover_files/1KGP_phase3_hg38_full_vcf_files/phase3.chr${i}.GRCh38.GT.crossmap.vcf.gz \
--threads 15 \
--min-ac 1 \
--force-samples \
-o 1kgp_1252_samples_chr${i}_variants.vcf

bgzip 1kgp_1252_samples_chr${i}_variants.vcf

tabix -p vcf 1kgp_1252_samples_chr${i}_variants.vcf.gz
done


> for i in $(seq 1 22) X Y 
do 
echo $i
bgzip 1kgp_1252_samples_chr${i}_variants.vcf
tabix -p vcf 1kgp_1252_samples_chr${i}_variants.vcf.gz
done

**FIlter those different samples vcf files for different maf cutoff**
> for i in $(seq 1 22) X Y 
do 
echo $i
bcftools view -q 0.05 ../all_variants/1kgp_1252_samples_chr${i}_variants.vcf.gz --threads 15 -Oz -o 1kgp_1252_samples_chr${i}_maf_filtered.vcf.gz
done


> bash /cn4/data4/venkatesh_data4/1KGP_genome_graphs/1KGP_vcf_files/liftover_files/variant_count_saturation_analysis/get_variant_stats.sh $(ls *vcf.gz)



**hypervariable region analysis**


for i in $(seq 1 22) X Y
do
echo "${i} is being processed"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 5' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "5 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 6' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "6 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 7' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "7 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 8' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "8 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 9' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "9 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 10' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "10 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 11' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "11 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/1KGP_complete_chr${i}.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 12' >> /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/recombination_analysis/chr${i}_hypervariable_node_ids.txt
echo "12 degree nodes captured"

echo "${i} is processed for hypervariable nodes"
done




**Create LD blocks using LDKit**


java -jar /apps/venkatesh/tools/LDkit/LDkit.jar --infile /cn4/data4/venkatesh_data4/1KGP_genome_graphs/1KGP_vcf_files/liftover_files/1KGP_phase3_hg38_full_vcf_files/phase3.chr21.GRCh38.GT.crossmap.vcf.gz \
--ws 100 \
--maf 0.01 \
--threads 20 \
--Intermediate yes


/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/hyper_variable_node_analysis/ld_blocks/ldkit_blocks/temp.txt


cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/maf_filtered_1KGP_graph/vg_manual_construct/gfa_files/10MB_subgraphs/gfa_files/chr6_30000001_40000001.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 5' 
echo "5 degree nodes captured"

cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/maf_filtered_1KGP_graph/vg_manual_construct/gfa_files/10MB_subgraphs/gfa_files/chr6_30000001_40000001.gfa | grep "L" | cut -f 2 | sort -h | uniq -c | awk '$1 == 6' 
echo "6 degree nodes captured"

**Hyp_var nodes id in common**
41783204
41435822
41435824


**Hyp_var nodes (pos,id) in common**
31340514,41435824 - 6
31340587,41435844 - 4
32492191,41493258 - 3
32589513,41507828 - 4





cat /cn4/data4/venkatesh_data4/1KGP_genome_graphs/maf_filtered_1KGP_graph/vg_manual_construct/gfa_files/10MB_subgraphs/gfa_files/chr6_30000001_40000001.gfa | grep "L" | cut -f 2 | grep "41435824"










for i in $(seq 1 22) X Y; 
do 
echo $i ; 
awk -F "," '{if($9 == 1) print $3}' chr${i}_gene_ann_check.csv | sort | uniq -c; 
done



for i in $(seq 1 22) X Y; 
do 
echo $i ; 
awk -F "," '{if($10 == 1) print $4}' chr${i}_gene_ann_check.csv | sort | uniq -c; 
done



for i in $(ls /data/shared_GI/GI_test_set_215);
do
ln -s /data/shared_GI/GI_test_set_215/${i} ${i}
done