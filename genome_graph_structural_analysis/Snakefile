'''
Snakefile to run the following workflow: 
   1) Create genome graphs from linear genome and set of variants using VG toolkit
   2) Extract GFA files from genome graphs
   3) Break complete genome graphs into [user chosen] bp binned subgraphs
   4) Measure the degree-based complexity of the subgraph 
   5) Visualize the genome graph using Circos

Required Inputs: 
   1) Linear reference genome in FASTA format. 
   2) Processed genetic variants in .vcf.gz format. Location: /datasets/input_files/processed_vcf_files
        VCF files should be names as CONTIG_NAME.vcf.gz. CONTIG_NAME should match the headers in the FASTA file
   3) Karyotype file for the chromosomes that the graph is made for.
Author: Venkatesh Kamaraj
Contact: ic11570@imail.iitm.ac.in
'''

configfile: "config/config.yaml"

def get_chromosome_names():
    '''
    Function that gets the names of chromosomes for which genetic variants are present in /datasets/input_files/processed_vcf_files folder
    '''
    import glob
    chromosome_names = []
    for file in glob.glob('datasets/input_files/processed_vcf_files/*.vcf.gz'):
        chromosome_names.append(file[41:].split(".")[0]) 
        # "datasets/input_files/processed_vcf_files/" = 41 characters && the chromosome name should be the name of the vcf file.
    return list(set(chromosome_names))

#CHROMS = ["chr1", "chr2"] # uncomment to run per sample by mentioning list of chromosomes names
CHROMS = get_chromosome_names()



#This is the rule that snakemake automatically tries to run. Will run all rules required to produce the input files.
rule all:
    input: 
        expand("datasets/genome_graph_visualisation/{plot}.png", plot = ["nodes_degree_distribution", "invariable_zones_number", "invariable_zones_max_length", "invariable_zones_median_length", "invariable_zones_normalised_length", "complete_genome_graph_visualisation"]),

        

# Create genome graphs for the input reference and variants using vg tolkit. Get output .vg files. 
rule construct_vg_files:
    input: 
        linear_reference=config["linear_reference_genome"],
        vcf_file="datasets/input_files/processed_vcf_files/{chromosome}.vcf.gz"
    output:
        "datasets/vg_files/{chromosome}.vg"
    conda:
        "envs/structure_gg.yaml"
    log:
        "logs/vg_construct/{chromosome}.log"
    threads: 25
    benchmark:
        "benchmarks/rule_1.{chromosome}.vg_construct.benchmark.txt"
    shell:
        "echo $(date '+%d/%m/%y    %H:%M:%S'): 'VG CONSTRUCT starts' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
        "(vg construct --alt-paths --region-is-chrom --region {wildcards.chromosome} --reference {input.linear_reference} --vcf {input.vcf_file} --threads {threads} --handle-sv > {output}) "         # The shell command
        "2>> {log}"                                                                                                  # write CLI log to log file

        " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'Created genome graph for {wildcards.chromosome} using VG CONSTRUCT' >> {log} "
        " && echo '-' >> {log}"                                                                                      # Add end time to log file


# Create the chromosome level genome graphs in GFA format. 
rule construct_gfa_files:
    input: 
        vg_file="datasets/vg_files/{chromosome}.vg"
    output:
        gfa_file="datasets/gfa_files/{chromosome}.gfa"
    conda:
        "envs/structure_gg.yaml"
    log:
        "logs/gfa_construct/{chromosome}.log"
    threads: 25
    benchmark:
        "benchmarks/rule_2.{chromosome}.gfa_construct.benchmark.txt"
    shell:
        "echo $(date '+%d/%m/%y    %H:%M:%S'): 'VG CONSTRUCT starts' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
        "(vg view -v {input.vg_file} -g > {output.gfa_file}) "         # The shell command
        "2>> {log}"                                                                                                  # write CLI log to log file

        " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'Created genome graph for {wildcards.chromosome} using VG CONSTRUCT' >> {log} "
        " && echo '-' >> {log}"                                                                                      # Add end time to log file


# perform structural analysis of chromosome-level genome graphs 
rule gg_structural_analysis:
    input: 
        gfa_file="datasets/gfa_files/{chromosome}.gfa",
    params:
        bin_width=config["Bin_Width"],
        save_lookup=config["Save_Lookup"],
        invar_cutoff=config["Invariance_cutoff"]
    output:
        lookup_output="datasets/node_ref_lookup/{chromosome}.csv",
        degree_output="datasets/chromosome_level_structural_analysis/all_nodes_summary/{chromosome}.csv",
        var_output="datasets/chromosome_level_structural_analysis/variable_nodes_summary/{chromosome}.tsv",
        hypvar_output="datasets/chromosome_level_structural_analysis/hypervariable_nodes_summary/{chromosome}.tsv",
        invar_output="datasets/chromosome_level_structural_analysis/invariable_nodes_summary/{chromosome}.csv"
    conda:
        "envs/structure_gg.yaml"
    log:
        "logs/gg_structural_analysis/{chromosome}.log"
    threads: 25
    benchmark:
        "benchmarks/rule_3.{chromosome}.gg_structural_analysis.benchmark.txt"
    shell:
        "echo $(date '+%d/%m/%y    %H:%M:%S'): 'Structural analysis starts for {wildcards.chromosome}' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
        "(python scripts/genome_graph_structural_analysis.py --gfa_file {input.gfa_file} --ref_path {wildcards.chromosome} --bin_width {params.bin_width} --invariance_cutoff {params.invar_cutoff} --save_lookup {params.save_lookup} --lookup_output {output.lookup_output} --degree_output {output.degree_output} --var_output {output.var_output} --hypvar_output {output.hypvar_output} --invar_output {output.invar_output}) "         # The shell command
        "2>> {log}"                                                                                                  # write CLI log to log file

        " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'Structural analysis ends for {wildcards.chromosome}' >> {log} "
        " && echo '-' >> {log}"                                                                                      # Add end time to log file



# Collate the output files from the previous rule
rule collate_files:
    input: 
        chromosome_variablity_files=expand("datasets/chromosome_level_structural_analysis/variable_nodes_summary/{chromosome}.tsv", chromosome = CHROMS),
        chromosome_hypvar_files=expand("datasets/chromosome_level_structural_analysis/hypervariable_nodes_summary/{chromosome}.tsv", chromosome = CHROMS),
        karyotype=config["karyotype_file"]
    output:
        variability_output="datasets/collated_results/genome_graph_variability.tsv",
        hypvar_output="datasets/collated_results/genome_graph_hypervariability.tsv",
        ideogram_conf="datasets/collated_results/circos_config_files/ideogram.conf",
        ideogram_label_conf="datasets/collated_results/circos_config_files/ideogram.label.conf",
        ticks_conf="datasets/collated_results/circos_config_files/ticks.conf",
        main_conf="datasets/collated_results/circos_config_files/circos_main.conf"
    params:
        ticks_multiplier=config['Ticks_Multiplier'],
        variability_color=config['Variability_Color'],
        hypvar_color=config['Hypervariability_Colour']
    conda:
        "envs/structure_gg.yaml"
    log:
        "logs/collate_files.log"
    threads: 25
    benchmark:
        "benchmarks/rule_4.collate_files.benchmark.txt"
    shell:
        "echo $(date '+%d/%m/%y    %H:%M:%S'): 'File collation starts' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
        "(Rscript scripts/collate_files.R --contig_variablity_files {input.chromosome_variablity_files} --contig_hypvar_files {input.chromosome_hypvar_files} --variability_output {output.variability_output} --hyp_var_output {output.hypvar_output} --karyotype {input.karyotype} --ideogram_conf {output.ideogram_conf} --ideogram_label_conf {output.ideogram_label_conf} --ticks_conf {output.ticks_conf} --main_conf {output.main_conf} --ticks_multiplier {params.ticks_multiplier} --variability_color {params.variability_color} --hypvar_color {params.hypvar_color}) "         # The shell command
        "2>> {log}"                                                                                                  # write CLI log to log file

        " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'File collation ends' >> {log} "
        " && echo '-' >> {log}"        




# Plot hypervariability plot, invariability plots
rule visualise_genome_graph_metrics:
    input: 
        chromosome_degree_files=expand("datasets/chromosome_level_structural_analysis/all_nodes_summary/{chromosome}.csv", chromosome = CHROMS),
        chromosome_invariablity_files=expand("datasets/chromosome_level_structural_analysis/invariable_nodes_summary/{chromosome}.csv", chromosome = CHROMS),
        karyotype=config["karyotype_file"],
    output:
        degree_output="datasets/collated_results/genome_graph_all_degrees_summary.csv",
        invariability_output="datasets/collated_results/genome_graph_invariability.csv",
        degree_plot="datasets/genome_graph_visualisation/nodes_degree_distribution.png",
        invar_number_plot="datasets/genome_graph_visualisation/invariable_zones_number.png",
        invar_length_plot="datasets/genome_graph_visualisation/invariable_zones_max_length.png",
        invar_median_plot="datasets/genome_graph_visualisation/invariable_zones_median_length.png",
        invariability_plot="datasets/genome_graph_visualisation/invariable_zones_normalised_length.png"
    conda:
        "envs/structure_gg.yaml"
    log:
        "logs/visualise_genome_graph_metrics.log"
    threads: 25
    benchmark:
        "benchmarks/rule_5.visualise_genome_graph_metrics.benchmark.txt"
    shell:
        "echo $(date '+%d/%m/%y    %H:%M:%S'): 'Visualisation of Genome Graphs Metrics Starts' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
        "(Rscript scripts/plot_summary.R --contig_degree_files {input.chromosome_degree_files} --contig_invariablity_files {input.chromosome_invariablity_files} --karyotype {input.karyotype} --degree_output {output.degree_output} --invariability_output {output.invariability_output} --degree_plot {output.degree_plot} --invar_number_plot {output.invar_number_plot} --invar_length_plot {output.invar_length_plot} --invar_median_plot {output.invar_median_plot} --invariability_plot {output.invariability_plot}) "         # The shell command
        "2>> {log}"                                                                                                  # write CLI log to log file

        " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'Visualisation of Genome Graph Metrics Ends' >> {log} "
        " && echo '-' >> {log}"                                                                                      # Add end time to log file
 

# Plot the entire genome graph in a circos plot
rule visualise_genome_graph_circos:
    input: 
        circos_conf="datasets/collated_results/circos_config_files/circos_main.conf",
    output:
        "datasets/genome_graph_visualisation/complete_genome_graph_visualisation.png"
    params:
        circos_plot="datasets/genome_graph_visualisation/complete_genome_graph_visualisation"
    conda:
        "envs/structure_gg.yaml"
    log:
        "logs/visualise_genome_graph_circos.log"
    threads: 25
    benchmark:
        "benchmarks/rule_6.visualise_genome_graph_circos.benchmark.txt"
    shell:
        "echo $(date '+%d/%m/%y    %H:%M:%S'): 'Visualisation of Overall Genome Graphs Starts' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
        "(circos  -outputdir ./   -outputfile  {params.circos_plot} -conf {input.circos_conf})" # shell command
        "2>> {log}"                                                                                                  # write CLI log to log file

        " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'Visualisation of Overall Genome Graph Ends' >> {log} "
        " && echo '-' >> {log}"                                                                                      # Add end time to log file
 

# # Index the chromosome level genome graphs in XG format. 
# rule index_xg_files:
#     input: 
#         vg_file="datasets/vg_files/{chromosome}.vg",
#         gfa_file="datasets/gfa_files/{chromosome}.gfa"
#     output:
#         xg_file="datasets/xg_files/{chromosome}.xg"
#     #conda:
#         #"envs/structure_gg.yaml"
#     log:
#         "logs/xg_index/{chromosome}.log"
#     threads: 25
#     benchmark:
#         "benchmarks/rule_3.{chromosome}.xg_index.benchmark.txt"
#     shell:
#         "echo $(date '+%d/%m/%y    %H:%M:%S'): 'VG CONSTRUCT starts' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
#         "(vg index -L -x {output.xg_file} {input.vg_file} --threads {threads}) "         # The shell command
#         "2>> {log}"                                                                                                  # write CLI log to log file

#         " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'Created genome graph for {wildcards.chromosome} using VG CONSTRUCT' >> {log} "
#         " && echo '-' >> {log}"                                                                                      # Add end time to log file


# # coordinate the node IDs od the chromosome level graphs 
# rule node_id_coordination:
#     input: 
#         vg_files=expand("datasets/vg_files/{chromosome}.vg", chromosome = CHROMS),
#     output:
#         mapping="datasets/vg_files/mapping",
#         mapping_backup="datasets/vg_files/mapping.backup"
#     #conda:
#         #"envs/structure_gg.yaml"
#     log:
#         "logs/vg_construct/node_id_coordination.log"
#     threads: 25
#     benchmark:
#         "benchmarks/rule_2.node_id_coordination.benchmark.txt"
#     shell:
#         "echo $(date '+%d/%m/%y    %H:%M:%S'): 'VG Node ID coordination started' >> {log} && echo '-' >> {log} &&"                 # Add start time to log file
        
#         "(vg ids -j -m {output.mapping} {input.vg_files}) "         # The shell command
#         "2>> {log}"

#         "&& cp {output.mapping} {output.mapping_backup}"                                                                                                  # write CLI log to log file

#         " && echo '-' >> {log} && echo $(date '+%d/%m/%y    %H:%M:%S'):  'VG Node ID coordination ended' >> {log} "
#         " && echo '-' >> {log}"                                                                                      # Add end time to log file


