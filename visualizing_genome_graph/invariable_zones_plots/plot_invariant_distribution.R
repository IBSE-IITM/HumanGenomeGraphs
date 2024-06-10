library(ggplot2)
library(dplyr)
library(scales)
library(reshape2)



setwd("/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/invariable_zones_plots")


# Load the CSV files
invariant_zones_complete_1kgp <- read.csv("/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/gg_metrics/complete_1KGP_graph_invariant_zones.csv")

invariant_zones_maf_1kgp <- read.csv("/cn4/data4/venkatesh_data4/1KGP_genome_graphs/maf_filtered_1KGP_graph/vg_manual_construct/visualizing_genome_graphs/gg_metrics/maf_filtered_1KGP_graph_invariant_zones.csv")


# Factrize the chromosome column
chr_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
"chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY" )

invariant_zones_complete_1kgp$Chr <- factor(invariant_zones_complete_1kgp$Chr, levels = chr_order)

invariant_zones_maf_1kgp$Chr <- factor(invariant_zones_maf_1kgp$Chr, levels = chr_order)



# Add reference genome graphs correspondingly
invariant_zones_maf_1kgp$Genome_Graph <- "1KGP_Common_Genome_Graph"
invariant_zones_complete_1kgp$Genome_Graph <- "1KGP_Complete_Genome_Graph"


# Combine the two table sto create dodge bar plots

combined_invariant_zones <- rbind(invariant_zones_complete_1kgp, invariant_zones_maf_1kgp)


# plot the total number of invariant zones in the pan-genomes
combined_invariant_zones %>%
,   count(Chr, Genome_Graph) %>%
    ggplot(aes(x = Chr, y = n, fill = Genome_Graph)) +
    geom_col(position = "dodge", color = "black") +
    xlab("Chromosome") +
    ylab("Zone Count") +
    ggtitle("Number of Invariant zones", subtitle = "In human pangenomes")+
    theme_classic()+
    scale_fill_manual(values = c("#0072B2", "#D55E00"))+
    # geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.2)
  theme(#legend.position="none",
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("number_of_invariant_zones.png", dpi = 600, width = 8, height = 8, units = "in")







# plot the maximum length of  invariant zones in the pan-genomes
combined_invariant_zones %>%
    group_by(Chr, Genome_Graph)  %>%
    summarise(max_zone_length = max(Zone_length_in_nodes)*32) %>%
    ggplot(aes(x = Chr, y = max_zone_length, fill = Genome_Graph)) +
    geom_col(position = "dodge", color = "black") +
    xlab("Chromosome") +
    ylab("Zone Length") +
    ggtitle("Maximum size of Invariant zones", subtitle = "In human pangenomes")+
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+ # millions
    theme_classic()+
    scale_fill_manual(values = c("#0072B2", "#D55E00"))+
    # geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.2)
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("max_size_of_invariant_zones.png", dpi = 600, width = 10, height = 8, units = "in")



# invariant_zones %>%
#     group_by(Chr) %>%
#     summarise(mean_length = mean(Zone_length_in_nodes)*32, median_length = median(Zone_length_in_nodes)*32) %>%
#     melt(id=c("Chr")) %>%
#     ggplot(aes(x = Chr, y = value, fill = variable)) +
#     geom_col(stat="identity", position=position_dodge()) +
#     xlab("Chromosome") +
#     ylab("Zone Length") +
#     ggtitle("Central value of Invariant zones", subtitle = "In each Chromosome")+
#     scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+ # millions
#     theme_classic()+
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.title = element_text(size = 20, face = "bold"),
#         plot.subtitle = element_text(size = 18),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 15))


combined_invariant_zones %>%
    group_by(Chr, Genome_Graph)  %>%
    summarise(median_zone_length = median(Zone_length_in_nodes)*32) %>%
    ggplot(aes(x = Chr, y = median_zone_length, fill = Genome_Graph)) +
    geom_col(position = "dodge", color = "black") +
    xlab("Chromosome") +
    ylab("Zone Length") +
    ggtitle("Median length of Invariant zones", subtitle = "In human pangenomes")+
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))+ # millions
    theme_classic()+
    scale_fill_manual(values = c("#0072B2", "#D55E00"))+
    # geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.2)
  theme(legend.position="none",
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("median_length_of_invariant_zones_without_legends.png", dpi = 600, width = 8, height = 8, units = "in")





# Get total length of invariant zones normalized by chr length
chr_lengths <- read.table("hg38.chrom.sizes")
colnames(chr_lengths) <- c("Chr", "Contig_Length")


total_invariant_lengths <- combined_invariant_zones %>%
                           group_by(Chr, Genome_Graph)  %>%
                           summarise(all_inv_zones_length = sum(Zone_length_in_nodes)*32)
    
    
    
total_invariant_lengths <- merge(total_invariant_lengths, chr_lengths, by = intersect(names(total_invariant_lengths), names(chr_lengths)), all = FALSE)  


total_invariant_lengths %>%
mutate(invariant_length_normalized = all_inv_zones_length/Contig_Length) %>%
ggplot(aes(x = Chr, y = invariant_length_normalized, fill = Genome_Graph)) +
geom_col(position = "dodge", color = "black") +
xlab("Chromosome") +
ylab("Ratio") +
ggtitle("Aggregate Invariability", subtitle = "Of each chromosome in the human pangenomes")+
theme_classic()+
    scale_fill_manual(values = c("#0072B2", "#D55E00"))+
    # geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.2)
  theme(legend.position="none",
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("normalized_len_of_invariant_zones_without_legends.png", dpi = 600, width = 8, height = 8, units = "in")




# Chromosomes that have same no. of invariant zomes in both human pangenomes

no_of_invariant_zones <- combined_invariant_zones %>%
                              count(Chr, Genome_Graph) 

common_numbers <- no_of_invariant_zones[no_of_invariant_zones$Genome_Graph == "1KGP_Common_Genome_Graph", ]
complete_numbers <- no_of_invariant_zones[no_of_invariant_zones$Genome_Graph == "1KGP_Complete_Genome_Graph", ]

sum(common_numbers$n == complete_numbers$n)

common_numbers$Chr[common_numbers$n == complete_numbers$n]
# chr3  chr4  chr6  chr11 chr12 chr13 chr14 chr17 chr18 chr20 chr21 chr22




# Chromosomes that have same length of invariant zomes in both human pangenomes

common_lengths <- total_invariant_lengths[total_invariant_lengths$Genome_Graph == "1KGP_Common_Genome_Graph", ]
complete_lengths <- total_invariant_lengths[total_invariant_lengths$Genome_Graph == "1KGP_Complete_Genome_Graph", ]

sum(common_lengths$n == complete_lengths$n)

common_lengths$Chr[common_lengths$n != complete_lengths$n]


