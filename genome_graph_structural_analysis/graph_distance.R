library(ggplot2)
library(dplyr)
library(purrr)
library(reshape2)


distance_calculator <- function(metric_count_df, pop_1, pop_2){
    # Function to calculate the relative difference in the metric (variability/hypervariability) for two populations 
    # metric_count_df: A dataframe containing the metrics from both the graphs
    # pop_1: Colname in the dataframe corresponding to graph_1. Should be a character 
    # pop_2: Colname in the dataframe corresponding to graph_2. Should be a character  

    ratio_pop_1 <- metric_count_df[ , pop_1]/sum(metric_count_df[ , pop_1]) # normalise each bin count with the (count for all bins)

    ratio_pop_2 <- metric_count_df[ , pop_2]/sum(metric_count_df[ , pop_2])
    
    abs_difference <- abs(ratio_pop_2 - ratio_pop_1)

    return(mean(abs_difference, na.rm = TRUE))
    
}


correlation_calculator <- function(metric_count_df, pop_1, pop_2){
    # Function to calculate the relative difference in the metric (variability/hypervariability) for two populations 
    # metric_count_df: A dataframe containing the metrics from both the graphs
    # pop_1: Colname in the dataframe corresponding to graph_1. Should be a character 
    # pop_2: Colname in the dataframe corresponding to graph_2. Should be a character  

    return(cor(metric_count_df[ , pop_1], metric_count_df[ , pop_2]))
}


# read all the population-specific 
EUR_variability <- read.csv("./population_specific_genome_graphs/EUR_genome_graph/datasets/collated_results/genome_graph_variability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
AMR_variability <- read.csv("./population_specific_genome_graphs/AMR_genome_graph/datasets/collated_results/genome_graph_variability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
AFR_variability <- read.csv("./population_specific_genome_graphs/AFR_genome_graph/datasets/collated_results/genome_graph_variability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
SAS_variability <- read.csv("./population_specific_genome_graphs/SAS_genome_graph/datasets/collated_results/genome_graph_variability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
EAS_variability <- read.csv("./population_specific_genome_graphs/EAS_genome_graph/datasets/collated_results/genome_graph_variability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))

# merge the variability values based on chr, start, end
total_variabiliy <- list(EUR_variability, AMR_variability, AFR_variability, SAS_variability, EAS_variability) %>% reduce(inner_join, by=c("Chr", "Start", "End"))
colnames(total_variabiliy) <- c("Chr", "Start", "End", "EUR", "AMR", "AFR", "SAS", "EAS")



# define a vector with all possible combinations of populations
population_combinations <- c("AFR-AMR", "AFR-EAS", "AFR-EUR", "AFR-SAS", 
                            "AMR-EAS", "AMR-EUR", "AMR-SAS", 
                            "EAS-EUR", "EAS-SAS",
                            "EUR-SAS")

# Get a list of all chromosomes in the human genome
chromosomes <- paste("chr", c(1:22, "X", "Y"), sep = "")
chromosomes <- c("All", chromosomes)

# Assing a 0 matrix for capturing the variability distance between all combinations of population specific genome graphs
variability_distance_matrix <- matrix(rep(0,240), 
                                      nrow = length(chromosomes), 
                                      ncol = length(population_combinations), 
                                      dimnames = list(chromosomes, population_combinations))


# Assing a 0 matrix for capturing the variability correlation between all combinations of population specific genome graphs
variability_correlation_matrix <- matrix(rep(0,240), 
                                      nrow = length(chromosomes), 
                                      ncol = length(population_combinations), 
                                      dimnames = list(chromosomes, population_combinations))


for(i in chromosomes){
    for(j in population_combinations){
        pop1 <- strsplit(j, "-")[[1]][1] # get EAS seperated from "EAS-SAS"
        pop2 <- strsplit(j, "-")[[1]][2] # get SAS seperated from "EAS-SAS"
        
        if(i == "All"){ # the metrics should be calculated for all chr
            variability_distance_matrix[i,j] <- distance_calculator(total_variabiliy, pop1, pop2)
            variability_correlation_matrix[i,j] <- correlation_calculator(total_variabiliy, pop1, pop2)
        }
        else{ # the metrics should be calculated for that particular chr
            variability_distance_matrix[i,j] <- distance_calculator(total_variabiliy[total_variabiliy$Chr == i, ], pop1, pop2)
            variability_correlation_matrix[i,j] <- correlation_calculator(total_variabiliy[total_variabiliy$Chr == i, ], pop1, pop2)
        }

        
    }
}



# Melt the matrix into long format
variability_distance_matrix_melted <- melt(variability_distance_matrix)
colnames(variability_distance_matrix_melted) <- c("Chromosome", "Population", "Distance")

# Create the heatmap
variability_distance_plot <- ggplot(variability_distance_matrix_melted, aes(Population, Chromosome, fill = Distance)) +
  geom_tile(color = "darkgrey") +  # Add white borders to the tiles
  geom_text(aes(label = round(Distance, 4)), vjust = 1) + # Add the values as text
  ggtitle("Distance in Variability", subtitle = "between poputation-specific human genome graphs") +
  xlab("Combination of Populations") +
  ylab("Chromosome(s)") +
  scale_fill_gradient(low = "white", high = "#00a6ff") + # Define color scale
  theme_minimal() + # Use a minimal theme
  coord_fixed() + # Ensure tiles are square
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 18, vjust = -1),
        axis.title.y = element_text(size = 18, vjust = 3),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())

ggsave("./population_specific_genome_graphs/comparatory_analysis/graph_distance_analysis/variability_distance_plot.png", 
        plot = variability_distance_plot, dpi = 600, width = 8, height = 20, units = "in")


# Melt the matrix into long format
variability_correlation_matrix_melted <- melt(variability_correlation_matrix)
colnames(variability_correlation_matrix_melted) <- c("Chromosome", "Population", "Correlation")

# Create the heatmap
variability_correlation_plot <- ggplot(variability_correlation_matrix_melted, aes(Population, Chromosome, fill = Correlation)) +
  geom_tile(color = "white") +  # Add white borders to the tiles
  ggtitle("Correlation in Variability", subtitle = "between poputation-specific human genome graphs") +
  xlab("Combination of Populations") +
  ylab("Chromosome(s)") +
  geom_text(aes(label = round(Correlation, 4)), vjust = 1, size = 2.5) + # Add the values as text
  scale_fill_gradient(high = "white", low = "#00a6ff") + # Define color scale
  theme_minimal() + # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  coord_fixed() + # Ensure tiles are square 
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 18, vjust = -1),
        axis.title.y = element_text(size = 18, vjust = 3),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())


ggsave("./population_specific_genome_graphs/comparatory_analysis/graph_distance_analysis/variability_correlation_plot.png", 
        plot = variability_correlation_plot, dpi = 600, width = 8, height = 20, units = "in")







# Now repeat the analysis for hypervariability




# read all the population-specific 
EUR_hypervariability <- read.csv("./population_specific_genome_graphs/EUR_genome_graph/datasets/collated_results/genome_graph_hypervariability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
AMR_hypervariability <- read.csv("./population_specific_genome_graphs/AMR_genome_graph/datasets/collated_results/genome_graph_hypervariability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
AFR_hypervariability <- read.csv("./population_specific_genome_graphs/AFR_genome_graph/datasets/collated_results/genome_graph_hypervariability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
SAS_hypervariability <- read.csv("./population_specific_genome_graphs/SAS_genome_graph/datasets/collated_results/genome_graph_hypervariability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))
EAS_hypervariability <- read.csv("./population_specific_genome_graphs/EAS_genome_graph/datasets/collated_results/genome_graph_hypervariability.tsv", sep = "\t", header = FALSE, col.names = c("Chr", "Start", "End", "Variability"))

# merge the hypervariability values based on chr, start, end
total_hypervariability <- list(EUR_hypervariability, AMR_hypervariability, AFR_hypervariability, SAS_hypervariability, EAS_hypervariability) %>% reduce(inner_join, by=c("Chr", "Start", "End"))
colnames(total_hypervariability) <- c("Chr", "Start", "End", "EUR", "AMR", "AFR", "SAS", "EAS")



# Assing a 0 matrix for capturing the hypervariability distance between all combinations of population specific genome graphs
hypervariability_distance_matrix <- matrix(rep(0,240), 
                                      nrow = length(chromosomes), 
                                      ncol = length(population_combinations), 
                                      dimnames = list(chromosomes, population_combinations))


# Assing a 0 matrix for capturing the hypervariability correlation between all combinations of population specific genome graphs
hypervariability_correlation_matrix <- matrix(rep(0,240), 
                                      nrow = length(chromosomes), 
                                      ncol = length(population_combinations), 
                                      dimnames = list(chromosomes, population_combinations))


for(i in chromosomes){
    for(j in population_combinations){
        pop1 <- strsplit(j, "-")[[1]][1] # get EAS seperated from "EAS-SAS"
        pop2 <- strsplit(j, "-")[[1]][2] # get SAS seperated from "EAS-SAS"
        
        if(i == "All"){ # the metrics should be calculated for all chr
            hypervariability_distance_matrix[i,j] <- distance_calculator(total_hypervariability, pop1, pop2)
            hypervariability_correlation_matrix[i,j] <- correlation_calculator(total_hypervariability, pop1, pop2)
        }
        else{ # the metrics should be calculated for that particular chr
            hypervariability_distance_matrix[i,j] <- distance_calculator(total_hypervariability[total_hypervariability$Chr == i, ], pop1, pop2)
            hypervariability_correlation_matrix[i,j] <- correlation_calculator(total_hypervariability[total_hypervariability$Chr == i, ], pop1, pop2)
        }

        
    }
}

# remove chrY as it has Nan for both distance and correlation. Hypervariability is 0 for all subgraphs in it
hypervariability_distance_matrix <- hypervariability_distance_matrix[-length(chromosomes), ]
hypervariability_correlation_matrix <- hypervariability_correlation_matrix[-length(chromosomes), ]


# Melt the matrix into long format
hypervariability_distance_matrix_melted <- melt(hypervariability_distance_matrix)
colnames(hypervariability_distance_matrix_melted) <- c("Chromosome", "Population", "Distance")

# Create the heatmap
hypervariability_distance_plot <- ggplot(hypervariability_distance_matrix_melted, aes(Population, Chromosome, fill = Distance)) +
  geom_tile(color = "darkgrey") +  # Add white borders to the tiles
  geom_text(aes(label = round(Distance, 4)), vjust = 1) + # Add the values as text
  ggtitle("Distance in Hypervariability", subtitle = "between poputation-specific human genome graphs") +
  xlab("Combination of Populations") +
  ylab("Chromosome(s)") +
  scale_fill_gradient(low = "white", high = "red") + # Define color scale
  theme_minimal() + # Use a minimal theme
  coord_fixed() + # Ensure tiles are square
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 18, vjust = -1),
        axis.title.y = element_text(size = 18, vjust = 3),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())

ggsave("./population_specific_genome_graphs/comparatory_analysis/graph_distance_analysis/hypervariability_distance_plot.png", 
        plot = hypervariability_distance_plot, dpi = 600, width = 8, height = 20, units = "in")


# Melt the matrix into long format
hypervariability_correlation_matrix_melted <- melt(hypervariability_correlation_matrix)
colnames(hypervariability_correlation_matrix_melted) <- c("Chromosome", "Population", "Correlation")

# Create the heatmap
hypervariability_correlation_plot <- ggplot(hypervariability_correlation_matrix_melted, aes(Population, Chromosome, fill = Correlation)) +
  geom_tile(color = "white") +  # Add white borders to the tiles
  ggtitle("Correlation in Hypervariability", subtitle = "between poputation-specific human genome graphs") +
  xlab("Combination of Populations") +
  ylab("Chromosome(s)") +
  geom_text(aes(label = round(Correlation, 4)), vjust = 1, size = 2.5) + # Add the values as text
  scale_fill_gradient(high = "white", low = "red") + # Define color scale
  theme_minimal() + # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  coord_fixed() + # Ensure tiles are square 
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 18, vjust = -1),
        axis.title.y = element_text(size = 18, vjust = 3),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_blank())


ggsave("./population_specific_genome_graphs/comparatory_analysis/graph_distance_analysis/hypervariability_correlation_plot.png", 
        plot = hypervariability_correlation_plot, dpi = 600, width = 8, height = 20, units = "in")




