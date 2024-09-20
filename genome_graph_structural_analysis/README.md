<h1>Structural analysis of Genome Graph</h1>

This repo contains codes to construct genome graphs and analyse its structure in-depth. <br>
Kindly execute the steps linearly and conserve the directory format. <br>
*Tested on Linux/Unix machines*

<br>

**Creating and activating the conda environnment**

The tools required for performing genome graph related analyses can be installed using conda [Install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). <br>
To install all the tools in one go, create a conda environment using the yml file provied with this repo at ```envs/structure_gg.yaml```.

1) Create the conda virtual environment: 
```
conda env create -f structure_gg.yaml
```

To verify if above step ran smoothly, run ```conda info --envs``` and check if ```structure_gg``` gets listed. 

2) Activate the environment created: 
```
conda activate structure_gg
```

<br>
<br> 

**Structural Analysis of genome graphs**

The directory contains a portable Snakemake based workflow that can be used for the structural analysis of genome graphs.

1) Move the processed bgzipped and indexed VCF files to ```datasets/input_files/processed_vcf_files/```. 
    * The variants in each chromosome should be present in individual VCF files. VCF files should be names as CHROMOSOME_NAME.vcf.gz
    * CHROMOSOME_NAME should match the headers in the FASTA file as well as the CHROM column in the VCF file
    * The corresponding index file should be named CHROMOSOME_NAME.vcf.gz.tbi

2) Edit the filepaths of the linear reference genome and the karyotype file in ```config/config.yaml```.

3) Move the karyotype file describing the chromosomes into ```datasets/input_files/```. 
    * This repo when cloned will contain the karyotype file for human reference genome GRCH38
    * Make sure the keys for individual chromosomes in the karyotype file match the CHROMOSOME_NAME in the VCF filename

4) Optionally the parameters of the analysis (like the bin_width for subgraph creation, node_cutoff for invariance etc.) can also be altered in ```config/config.yaml```.

5) Activate the ```structure_gg``` conda environment.
```
conda activate structure_gg
```

6) Perform a Snakemake dry run to ensure that the workflow works.
```
snakemake -n
``` 


7) Once the dry run has run successfully, the structural analysis can be started with the command below:
```
snakemake --cores [integer]--jobs [integer] -p 
```

> Note: Edit the parameters accordingly <br>
> --cores : Use at most N CPU cores/jobs in parallel
> --jobs : max number of jobs to be performed parallely. <br>
> --use-conda: to enable conda through job scheduler <br>
> *Optional parameters:* <br>
> -p & -r: print the command and the reason to run the command
> --forceall: Mandatorily run all the rules in the Snakefile <br>
> --rerun-incomplete: If the analysis stops before completion, this parameter can be added to recreate any incomplete files. <br>


Kindly look into [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/) for more details.


*Note: As circos plots are used to create a single figure depiction of the entire genome graph with all the chromosomes, the bin_width has to be altered such that there is less than 360 bins in total, across all the chromosomes. The visualization of the circos plot can be neglected to avoid any resolution limitation of this analysis, allowing the bin_width to be as small as possible.*

8) For increased resolution of analysis, the workflow can be forced to omit the creation of single figure depiction of the genome graph with the following command:
```
snakemake -R --until visualise_genome_graph_metrics collate_files
```
 <br><br>

**Outputs Generated**
 <br> <br>

1) ```datasets/vg_files/{chromosome}.vg``` : A genome graph for each chromosome in .vg format (binary file) <br> <br>

2) ```datasets/gfa_files/{chromosome}.gfa``` : A genome graph for each chromosome in .gfa format (human readable file) <br> <br>

3) ```datasets/node_ref_lookup/{chromosome}.csv``` : (Optional) A reference table for each chromosome that relates the node ids from its genome graph to genomic location <br> <br>

4) ```datasets/chromosome_level_structural_analysis/all_nodes_summary/{chromosome}.csv``` : For every chromosome, count of nodes with each out-degree till 30 in each bin <br> <br>

5) ```datasets/chromosome_level_structural_analysis/variable_nodes_summary/{chromosome}.tsv``` : For every chromosome, count of variable nodes in each bin <br> <br>

6) ```datasets/chromosome_level_structural_analysis/hypervariable_nodes_summary/{chromosome}.tsv``` : For every chromosome, count of hypervariable nodes in each bin <br> <br>

7) ```datasets/chromosome_level_structural_analysis/invariable_nodes_summary/{chromosome}.csv``` : For every chromosome, information of invariable regions <br> <br>

8) ```datasets/collated_results/genome_graph_variability.tsv``` : Count of variability in each bin, collated for all the chromosomes <br> <br>

9) ```datasets/collated_results/genome_graph_hypervariability.tsv``` : Count of hypervariability in each bin, collated for all the chromosomes <br> <br>

10) ```datasets/collated_results/genome_graph_all_degrees_summary.csv``` : Count of all nodes with each out-degree till 30 in each bin, collated for all the chromosomes <br> <br>

11) ```datasets/collated_results/genome_graph_invariability.csv``` : Information about invariability collated for all the chromosomes <br> <br>

12) ```datasets/genome_graph_visualisation/nodes_degree_distribution.png``` : Plot depicting the distribution of nodes w.r.t. out-degree for all chromosomes <br> <br>

13) ```datasets/genome_graph_visualisation/invariable_zones_number.png``` : Plot depicting number of invariable zones in each chromosome <br> <br>

14) ```datasets/genome_graph_visualisation/invariable_zones_max_length.png``` : Plot depicting length of the largest invariable zone in each chromosome <br> <br>

15) ```datasets/genome_graph_visualisation/invariable_zones_median_length.png``` : Plot depicting median length of invariable zones in each chromosome <br> <br>

16) ```datasets/genome_graph_visualisation/invariable_zones_normalised_length.png``` : Plot depicting aggregate invariability in each chromosome <br> <br>

17) ```datasets/genome_graph_visualisation/complete_genome_graph_visualisation.png``` : A single figure depicting the entire genome graph with all the chromosomes.

<br> <br>
