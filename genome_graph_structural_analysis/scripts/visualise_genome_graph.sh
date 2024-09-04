# Script to visualise the genome graph as a circos and plot summary visualisations
# 
# Author: Venkatesh Kamaraj
# Contact: ic11570@imail.iitm.ac.in

#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 
    --contig_degree_files CONTIG_DEGREE_FILES [...] list of files describing contig level degree-bin summary 
    --contig_invariablity_files CONTIG_INVARIABLITY_FILES [...] list of files describing contig level invariability 
    --karyotype KARYOTYPE Karyotype filepath 
    --degree_output DEGREE_OUTPUT Output filepath for the combined degree-bin summary file 
    --invariability_output INVARIABILITY_OUTPUT Output filepath for the combined invariability file
    --degree_plot DEGREE_PLOT Filepath for the combined degree-bin summary plot 
    --invar_number_plot INVAR_NUMBER_PLOT Filepath for the number of invariable zones plot
    --invar_length_plot INVAR_LENGTH_PLOT Filepath for the Maximum Length of invariable zones plot
    --invar_median_plot INVAR_MEDIAN_PLOT Filepath for the Median Length of invariable zones plot
    --invariability_plot INVARIABILITY_PLOT Filepath for the Aggregate Invariability plot
    --circos_conf CIRCOS_CONF Circos main Configuration filepath
    --circos_plot CIRCOS_PLOT" Filepath for circos plots depicting the entire genome graph
    exit 1
}

# Parse command-line arguments
contig_degree_files=()
contig_invariablity_files=()
while [[ "$#" -gt 0 ]]; do
    echo "outer loop $1"
    case $1 in
        --contig_degree_files) 
            shift
            while [[ "$#" -gt 0 && ! "$1" =~ ^-- ]]; do
                echo "inner loop $1"
                contig_degree_files+=("$1")
                shift
            done
            ;;
        --contig_invariablity_files) 
            echo "invar_files_encountered"
            shift
            while [[ "$#" -gt 0 && ! "$1" =~ ^-- ]]; do
                contig_invariablity_files+=("$1")
                shift
            done
            ;;
        --karyotype) karyotype="$2"; shift ;;
        --degree_output) degree_output="$2"; shift ;;
        --invariability_output) invariability_output="$2"; shift ;;
        --degree_plot) degree_plot="$2"; shift ;;
        --invar_number_plot) invar_number_plot="$2"; shift ;;
        --invar_length_plot) invar_length_plot="$2"; shift ;;
        --invar_median_plot) invar_median_plot="$2"; shift ;;
        --invariability_plot) invariability_plot="$2"; shift ;;
        --circos_conf) circos_conf="$2"; shift ;;
        --circos_plot) circos_plot="$2"; shift ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done


# Check if all required arguments are provided
if [ ${#contig_degree_files[@]} -eq 0 ] || [ ${#contig_invariablity_files[@]} -eq 0 ] || [ -z "$karyotype" ] || [ -z "$degree_output" ] || [ -z "$invariability_output" ] || [ -z "$degree_plot" ] || [ -z "$invar_number_plot" ] || [ -z "$invar_length_plot" ] || [ -z "$invar_median_plot" ] || [ -z "$invariability_plot" ] || [ -z "$circos_conf" ]|| [ -z "$circos_plot" ]; then
    echo "Error: Missing arguments"
    usage
fi

# Call the R script that plots the hypervariability and invariability plots
Rscript scripts/plot_summary.R --contig_degree_files "$contig_degree_files" --contig_invariablity_files "$contig_invariablity_files" --karyotype "$karyotype" --degree_output "$degree_output" --invariability_output "$invariability_output" --degree_plot "$degree_plot" --invar_number_plot "$invar_number_plot" --invar_length_plot "$invar_length_plot" --invar_median_plot "$invar_median_plot" --invariability_plot "$invariability_plot"


# Call Circos to plot the entire genome graph in a single figure
circos  -outputdir ./   -outputfile  "$circos_plot" -conf "$circos_conf"




#bash scripts/visualise_genome_graph.sh --contig_degree_files datasets/contig_level_structural_analysis/all_nodes_summary/chr22.csv datasets/contig_level_structural_analysis/all_nodes_summary/chr21.csv --karyotype /cn4/data4/venkatesh_data4/genome_graph_structural_analysis/datasets/input_files/karyotype.human.hg38.txt --contig_invariablity_files datasets/contig_level_structural_analysis/invariable_nodes_summary/chr22.csv datasets/contig_level_structural_analysis/invariable_nodes_summary/chr21.csv --degree_output datasets/collated_results/genome_graph_all_degrees_summary.csv --invariability_output datasets/collated_results/genome_graph_invariability.csv --degree_plot datasets/genome_graph_visualisation/genome_graph_degree_distribution.png --invar_number_plot datasets/genome_graph_visualisation/invariable_zones_number.png --invar_length_plot datasets/genome_graph_visualisation/invariable_zones_max_length.png --invar_median_plot datasets/genome_graph_visualisation/invariable_zones_median_length.png --invariability_plot datasets/genome_graph_visualisation/invariable_zones_normalised_length.png --circos_conf datasets/collated_results/circos_config_files/circos_main.conf --circos_plot datasets/genome_graph_visualisation/genome_graph_visualisation