# Script to align reads (extracted from sorted bam files) to an existing whole genome graph and call variants
# Usage: bash vg_giraffe_freebayes_pipeline.sh [bam_filenames]


# Command below is for error handling in bash
set -euo pipefail

# Set the number of threads to be used by the pipelne
n_threads=25

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate free_sam_bcf

# or Run "conda activate free_sam_bcf" before running this script

# Set path to the existing graph indices
gbz_index=/apps/venkatesh/genome_graph/vg_giraffe_pipeline/vg_giraffe_graph_files/hg38-pangenome.giraffe.gbz
xg_index=/apps/venkatesh/genome_graph/vg_giraffe_pipeline/vg_giraffe_graph_files/hg38-pangenome.xg
min_index=/apps/venkatesh/genome_graph/vg_giraffe_pipeline/vg_giraffe_graph_files/hg38-pangenome.min
dist_index=/apps/venkatesh/genome_graph/vg_giraffe_pipeline/vg_giraffe_graph_files/hg38-pangenome.dist

# Set path to the existing reference genomes
reference_fasta_index=/apps/venkatesh/genome_graph/vg_giraffe_pipeline/GRCH38_contigs/fasta_files_renamed_from_fna/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fa.fai
reference_fasta=/apps/venkatesh/genome_graph/vg_giraffe_pipeline/GRCH38_contigs/fasta_files_renamed_from_fna/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fa

# perform the following tasks for all bam files
for bam_file in "$@"
do 
    echo "$(date "+%d/%m/%y    %H:%M:%S"): Working with $bam_file"
    
    # Extract sample name from the bam file name
    sample_name=$(echo $bam_file |awk '{gsub(".bam", "")}1' | awk '{gsub("/apps/venkatesh/CCMB_samples/", "")}1')
    
    # Create a folder for each sample
    mkdir $sample_name
    cd $sample_name

    # create a log file for each sample
    touch $sample_name.log

    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Pipeline started for $sample_name with $n_threads threads"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Pipeline started for $sample_name with $n_threads threads in $(hostname)" >> $sample_name.log

    # make a subfolder to store files outputted in each task
    mkdir alignment
    mkdir bam_processing
    mkdir variant_calling
    
    # Map reads to whole genome graph with vg giraffe
    samtools bam2fq $bam_file | \
    vg giraffe -if - -Z $gbz_index \
    -m $min_index \
    -d $dist_index \
    -t $n_threads > ./alignment/vg_giraffe_mapping_$sample_name.gam
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Mapped to genome graph for $sample_name with VG GIRAFFE"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Mapped to genome graph for $sample_name with VG GIRAFFE" >> $sample_name.log
    
    # surject gam file to bam file for calling variants 
    vg surject -x $xg_index \
    -b ./alignment/vg_giraffe_mapping_$sample_name.gam \
    -t $n_threads > ./alignment/vg_giraffe_mapping_surjected_$sample_name.bam
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  GAM is surjected to BAM for $sample_name with VG"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  GAM is surjected to BAM for $sample_name with VG" >> $sample_name.log

    # Add Read group names to the BAM file
    samtools addreplacerg ./alignment/vg_giraffe_mapping_surjected_$sample_name.bam \
    --threads $n_threads \
    -r ID:normal -r SM:$sample_name \
    -o ./bam_processing/rg_added_wgg_mapping_$sample_name.bam
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  RGIDs added to BAM for $sample_name with SAMTOOLS"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  RGIDs added to BAM for $sample_name with SAMTOOLS" >> $sample_name.log

    # Sort the BAM file
    samtools sort ./bam_processing/rg_added_wgg_mapping_$sample_name.bam \
    --threads $n_threads \
    -o ./bam_processing/sorted_wgg_mapping_$sample_name.bam
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  BAM file is sorted for $sample_name with SAMTOOLS"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  BAM file is sorted for $sample_name with SAMTOOLS" >> $sample_name.log

    # Mark Duplicates
    samtools markdup ./bam_processing/sorted_wgg_mapping_$sample_name.bam \
    ./bam_processing/markdup_sorted_wgg_mapping_$sample_name.bam -@ $n_threads
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Marked Duplicates for $sample_name with SAMTOOLS"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Marked Duplicates for $sample_name with SAMTOOLS" >> $sample_name.log

    # Index BAM file
    samtools index ./bam_processing/markdup_sorted_wgg_mapping_$sample_name.bam -@ $n_threads
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Index BAM for $sample_name with SAMTOOLS"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Index BAM for $sample_name with SAMTOOLS" >> $sample_name.log

    # Use FreeBayes-Parallel to call variants
    freebayes-parallel <(fasta_generate_regions.py $reference_fasta_index 100000) $n_threads \
    -f $reference_fasta \
    ./bam_processing/markdup_sorted_wgg_mapping_$sample_name.bam \
    --min-mapping-quality 10 \
    --min-base-quality 10 \
    --min-coverage 4 > ./variant_calling/wgg_graph_variants_raw_parallel_freebayes_$sample_name.vcf
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Variants called for $sample_name with FreeBayes"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Variants called for $sample_name with FreeBayes" >> $sample_name.log

    # Convert MNPs to SNPs
    vcfallelicprimitives \
    -kg ./variant_calling/wgg_graph_variants_raw_parallel_freebayes_$sample_name.vcf \
    > ./variant_calling/wgg_graph_variants_processed_parallel_freebayes_$sample_name.vcf
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Variants processed for $sample_name with vcflib"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Variants processed for $sample_name with vcflib" >> $sample_name.log

    # Filter variants for quality
    vcftools --vcf ./variant_calling/wgg_graph_variants_processed_parallel_freebayes_$sample_name.vcf \
    --minQ 10 --recode \
    --out ./variant_calling/qual_10_variants_$sample_name.vcf
    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Variants filtered for $sample_name with vcftools"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Variants filtered for $sample_name with vcftools" >> $sample_name.log

    #echo "$(date "+%d/%m/%y    %H:%M:%S"):  Pipeline for $sample_name ran Successfully"
    echo "$(date "+%d/%m/%y    %H:%M:%S"):  Pipeline for $sample_name ran Successfully in $(hostname)" >> $sample_name.log
    
    # delete intermediete files
    rm ./variant_calling/wgg_graph_variants_processed_parallel_freebayes_$sample_name.vcf
    rm ./variant_calling/wgg_graph_variants_raw_parallel_freebayes_$sample_name.vcf
    rm ./bam_processing/sorted_wgg_mapping_$sample_name.bam
    rm ./bam_processing/rg_added_wgg_mapping_$sample_name.bam
    rm ./alignment/vg_giraffe_mapping_surjected_$sample_name.bam
    
    # get out of the folder containing the sample files
    cd ..
done

