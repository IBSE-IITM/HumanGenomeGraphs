# R script to plot the folloeing summary of genome graph structural analysis: 
#    1) Degree distribution of hypervariable nodes
#    2) Numberof invariable regions, Median Length of invariable regions, and Average invariability of contigs 
# 
# Author: Venkatesh Kamaraj
# Contact: ic11570@imail.iitm.ac.in


library(argparse)
library(ggplot2)
library(dplyr)
library(scales)
library(reshape2)

# Function to read CSV files
combine_csv_files <- function(files) {
  combined_data <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  return(combined_data)
}


# Argument parser setup
parser <- ArgumentParser(description = "Prepare the files for genome graph visualisation")

# Adding argument for list of CSV files
parser$add_argument("--contig_degree_files", nargs = "+", required = TRUE, help = "list of files describing contig level degree-bin summary")
parser$add_argument("--contig_invariablity_files", nargs = "+", required = TRUE, help = "list of files describing contig level invariability")
parser$add_argument("--karyotype", required = TRUE, help = "Circos Files: karyotype filepath")
parser$add_argument("--degree_output", required = TRUE, help = "Output filepath for the combined degree-bin summary file")
parser$add_argument("--invariability_output", required = TRUE, help = "Output filepath for the combined invariability file")
parser$add_argument("--degree_plot", required = TRUE, help = "Filepath for the combined degree-bin summary plot")
parser$add_argument("--invar_number_plot", required = TRUE, help = "Filepath for the number of invariable zones plot")
parser$add_argument("--invar_length_plot", required = TRUE, help = "Filepath for the Maximum Length of invariable zones plot")
parser$add_argument("--invar_median_plot", required = TRUE, help = "Filepath for the Median Length of invariable zones plot")
parser$add_argument("--invariability_plot", required = TRUE, help = "Filepath for the Aggregate Invariability plot")

# Parse arguments
args <- parser$parse_args()

# # First combine the contig level summary files to one combined file.

# Combine the CSV files
combined_degree_data <- combine_csv_files(args$contig_degree_files)
combined_invariability_data <- combine_csv_files(args$contig_invariablity_files)


# Save the combined files to the specified output file
write.csv(combined_degree_data, args$degree_output, row.names = FALSE)
write.csv(combined_invariability_data, args$invariability_output, row.names = FALSE)

cat( Sys.time(), " DONE: Invariability and Degree Distribution of genome graph collated", "\n")

# Plot hypervariable nodes' degree distribution

# Read chromosome names from karyotype file
karyotype <- read.table(args$karyotype)
chr_order <- karyotype$V7


# chr_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
# "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
# "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY" )
# Get chromosome level summary
# chr_summary <- combined_degree_data %>%
#                 group_by(chromosome) %>%
#                 summarize(total_deg_5 = sum(deg_5_count),  total_deg_6 = sum(deg_6_count),  total_deg_7 = sum(deg_7_count),  
#                 total_deg_8 = sum(deg_8_count),  total_deg_9 = sum(deg_9_count),total_deg_10 = sum(deg_10_count),
#                 total_deg_11 = sum(deg_11_count),total_deg_12 = sum(deg_12_count),total_deg_13 = sum(deg_13_count),
#                 total_deg_14 = sum(deg_14_count),total_deg_15 = sum(deg_15_count),total_deg_16 = sum(deg_16_count),
#                 total_deg_17 = sum(deg_17_count),total_deg_18 = sum(deg_18_count),total_deg_19 = sum(deg_19_count),
#                 total_deg_20 = sum(deg_20_count),total_deg_21 = sum(deg_21_count),total_deg_22 = sum(deg_22_count),
#                 total_deg_23 = sum(deg_23_count),total_deg_24 = sum(deg_24_count),total_deg_25 = sum(deg_25_count),
#                 total_deg_26 = sum(deg_26_count),total_deg_27 = sum(deg_27_count),total_deg_28 = sum(deg_28_count),
#                 total_deg_29 = sum(deg_29_count),total_deg_30 = sum(deg_30_count))


high_degrees <- as.character(seq(5,30,1))
high_complexity_count <- c(sum(combined_degree_data$deg_5_count), sum(combined_degree_data$deg_6_count), 
sum(combined_degree_data$deg_7_count), sum(combined_degree_data$deg_8_count), sum(combined_degree_data$deg_9_count), 
sum(combined_degree_data$deg_10_count), sum(combined_degree_data$deg_11_count), sum(combined_degree_data$deg_12_count), 
sum(combined_degree_data$deg_13_count), sum(combined_degree_data$deg_14_count), sum(combined_degree_data$deg_15_count), 
sum(combined_degree_data$deg_16_count), sum(combined_degree_data$deg_17_count), sum(combined_degree_data$deg_18_count), 
sum(combined_degree_data$deg_19_count), sum(combined_degree_data$deg_20_count), sum(combined_degree_data$deg_21_count), 
sum(combined_degree_data$deg_22_count), sum(combined_degree_data$deg_23_count), sum(combined_degree_data$deg_24_count), 
sum(combined_degree_data$deg_25_count), sum(combined_degree_data$deg_26_count), sum(combined_degree_data$deg_27_count), 
sum(combined_degree_data$deg_28_count), sum(combined_degree_data$deg_29_count), sum(combined_degree_data$deg_30_count))


combined_results <- data.frame(Degree = high_degrees, Count = high_complexity_count)
combined_results$Degree <- factor(combined_results$Degree, levels = high_degrees)

hypervariability_plot <- ggplot(combined_results, aes(x=Degree,y=Count))+
                            geom_col(fill = "#56B4E9", color = "black")+
                            geom_text(aes(label = Count), vjust = -0.8, size = 5.5)+
                            #scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
                            ggtitle("Distribution of Hypervariable nodes", subtitle = "Nodes with out_degree >= 5")+
                            theme_classic()+
                            theme(axis.text.x = element_text(size = 18),
                                  axis.text.y = element_text(size = 18),
                                  axis.title.x = element_text(size = 18),
                                  axis.title.y = element_text(size = 18),
                                  plot.title = element_text(size = 20, face = "bold"),
                                  plot.subtitle = element_text(size = 18),
                                  legend.title = element_blank()) 

ggsave(args$degree_plot, plot = hypervariability_plot, dpi = 600, width = 16, height = 8, units = "in")

cat( Sys.time(), " DONE: Degree Distribution of genome graph visualised", "\n")

# Plot the invariable regions in the genome graph
combined_invariability_data$Chr <- factor(combined_invariability_data$Chr, levels = chr_order)


invar_number <- combined_invariability_data %>%
    count(Chr) %>%
    ggplot(aes(x = Chr, y = n)) +
    geom_col(fill = "lightgreen", color = "black") +
    xlab("Chromosome") +
    ylab("Zone Count") +
    ggtitle("Distribution of Invariable Zones", subtitle = "No. of Invariant zones")+
    theme_classic()+
  theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())

ggsave(args$invar_number_plot, plot = invar_number, dpi = 600, width = 8, height = 8, units = "in")




invar_length <- combined_invariability_data %>%
    group_by(Chr) %>%
    summarise(max_zone_length = max(Zone_length_in_nodes)*32) %>%
    ggplot(aes(x = Chr, y = max_zone_length)) +
    geom_col(fill = "lightgreen", color = "black") +
    xlab("Chromosome") +
    ylab("Zone Length") +
    ggtitle("Distribution of Invariable Zones", subtitle = "Length of the largest Invariant zones")+
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+ # millions
    theme_classic()+
  theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())

ggsave(args$invar_length_plot, plot = invar_length, dpi = 600, width = 8, height = 8, units = "in")




invar_median <- combined_invariability_data %>%
    group_by(Chr) %>%
    summarise(median_zone_length = median(Zone_length_in_nodes)*32) %>%
    ggplot(aes(x = Chr, y = median_zone_length)) +
    geom_col(fill = "lightgreen", color = "black") +
    xlab("Chromosome") +
    ylab("Zone Length") +
    ggtitle("Distribution of Invariable Zones", subtitle = "Median length Invariant zones")+
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+ # millions
    theme_classic()+
  theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())


ggsave(args$invar_median_plot, plot = invar_median, dpi = 600, width = 8, height = 8, units = "in")






# Get total length of invariant zones normalized by chr length

total_invariant_lengths <- combined_invariability_data %>%
                            group_by(Chr) %>%
                            summarise(all_inv_zones_length = sum(Zone_length_in_nodes)*32)
    
total_invariant_lengths <- merge(total_invariant_lengths, karyotype, by.x = "Chr", by.y = "V7", all = FALSE)  


invariability <- total_invariant_lengths %>%
mutate(invariant_length_normalized = all_inv_zones_length/V6) %>%
ggplot(aes(x = Chr, y = invariant_length_normalized)) +
geom_col(fill = "lightgreen", color = "black") +
xlab("Chromosome") +
ylab("Zone Length") +
ggtitle("Distribution of Invariable Zones", subtitle = "Normalized length of Invariant zones")+
theme_classic()+
theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 18),
      legend.title = element_blank())

ggsave(args$invariability_plot, plot = invariability, dpi = 600, width = 8, height = 8, units = "in")

cat( Sys.time(), " DONE: Invariability of genome graph visualised", "\n")