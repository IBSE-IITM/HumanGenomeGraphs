library(ggplot2)
library(dplyr)

bin_summary <- read.csv("/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/gg_metrics/10mb_summary_of_complete_1kgp_gg.csv")

# chr_summary <- bin_summary %>%
#                 group_by(chromosome) %>%
#                 summarize(total_deg_6 = sum(deg_6_count),
#                           total_deg_7 = sum(deg_7_count),
#                           total_deg_8 = sum(deg_8_count),
#                           total_deg_9 = sum(deg_9_count),
#                           total_deg_10 = sum(deg_10_count),
#                           total_deg_11 = sum(deg_11_count),
#                           total_deg_12 = sum(deg_12_count),
#                           total_deg_13 = sum(deg_13_count))




bin_summary[bin_summary$deg_12_count >= 1, ]
# chr1 has the most complex nodes: 2 nodes with out_deg = 12


high_degrees <- c("5", "6", "7", "8", "9", "10", "11","12")
high_complexity_count <- c(sum(bin_summary$deg_5_count), sum(bin_summary$deg_6_count), 
sum(bin_summary$deg_7_count), sum(bin_summary$deg_8_count), sum(bin_summary$deg_9_count), 
sum(bin_summary$deg_10_count), sum(bin_summary$deg_11_count), sum(bin_summary$deg_12_count))

combined_results <- data.frame(Degree = high_degrees, Count = high_complexity_count)
combined_results$Degree <- factor(combined_results$Degree, levels = high_degrees)

ggplot(combined_results, aes(x=Degree,y=Count))+
geom_col(fill = "#009E73", color = "black")+
geom_text(aes(label = Count), vjust = -1)+
ggtitle("1KGP Complete Genome Graph", subtitle = "Count of Hypervariable nodes")+
theme_classic()+
    theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())

ggsave("/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/hyper_variable_region_plots/complete_1kgp_hypervariable_nodes_dist.png", dpi = 600, width = 8, height = 8, units = "in")
