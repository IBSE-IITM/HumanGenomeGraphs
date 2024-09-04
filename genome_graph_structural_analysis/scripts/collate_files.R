
# R script to prepare files for visualizing the results of genome graphs structural analysis: 
#    1) Combine contig level summary created by genome_graph_structural_analysis.py
#    2) Create configuration files required for Circos plot generation
# 
# Author: Venkatesh Kamaraj
# Contact: ic11570@imail.iitm.ac.in



library(argparse)

# Function to read TSV files
read_tsv_files <- function(file) {
    data <- read.csv(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    return(data)
}


# Function to combine TSV files
combine_tsv_files <- function(files) {
  combined_data <- do.call(rbind, lapply(files, read_tsv_files))
  return(combined_data)
}


# Argument parser setup
parser <- ArgumentParser(description = "Prepare the files for genome graph visualisation")

# Adding argument for list of CSV files
parser$add_argument("--contig_variablity_files", nargs = "+", required = TRUE, help = "list of files describing contig level variability")
parser$add_argument("--contig_hypvar_files", nargs = "+", required = TRUE, help = "list of files describing contig level hypervariability")
parser$add_argument("--variability_output", required = TRUE, help = "Output filepath for the combined variability file")
parser$add_argument("--hyp_var_output", required = TRUE, help = "Output filepath for the combined hypervariability file")
parser$add_argument("--karyotype", required = TRUE, help = "Circos Files: karyotype filepath")
parser$add_argument("--ideogram_conf", required = TRUE, help = "Circos Files: Ideogram Configuration filepath")
parser$add_argument("--ideogram_label_conf", required = TRUE, help = "Circos Files: Ideogram Labels Configuration filepath")
parser$add_argument("--ticks_conf", required = TRUE, help = "Circos Files: Ticks Configuration filepath")
parser$add_argument("--main_conf", required = TRUE, help = "Circos Files: Main Configuration filepath")

parser$add_argument("--ticks_multiplier", required = TRUE, help = "Circos Parameter: Multiplier for Ticks")
parser$add_argument("--variability_color", required = TRUE, help = "Circos Parameter: Colour for the variability track")
parser$add_argument("--hypvar_color", required = TRUE, help = "Circos Parameter: Colour for the hypervariability track")
# Parse arguments
args <- parser$parse_args()

# # First combine the contig level summary files to one combined file.

# Combine the CSV files
combined_variability_data <- combine_tsv_files(args$contig_variablity_files)
combined_hypvar_data <- combine_tsv_files(args$contig_hypvar_files)


# Save the combined files to the specified output file
write.table(combined_variability_data, file =  args$variability_output, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(combined_hypvar_data, file =  args$hyp_var_output, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

cat( Sys.time(), " DONE: Variability and Hypervariability of genome graph collated", "\n")

# # Make configuration files for Circos

# Specify the lines to be written to the ideogra.conf file
ideogram_conf <- c(
  "<ideogram>",
  "",
  "    <spacing>",
  "        # spacing between ideograms",
  "        default = 0.005r",
  "    </spacing>",
  "",
  "    # ideogram position, thickness and fill",
  "    radius           = 0.90r",
  "    thickness        = 50p",
  "    fill             = yes",
  "",
  "    #stroke_thickness = 1",
  "    #stroke_color     = black",
  "",
  "    # configuration of ideogram labels",
  "    <<include ideogram.label.conf>>",
  "",
  "",
  "</ideogram>"
)

# Write the lines to the ideogram.conf file
writeLines(ideogram_conf, args$ideogram_conf)



# Specify the lines to be written to the ideogram.labels.conf file
ideogram_label_conf <- c(
  "# optional",
  "",
  "show_label       = yes",
  "label_radius     = dims(ideogram,radius) + 0.075r",
  "label_size       = 36",
  "label_parallel   = yes"
)

# Write the lines to ideogram.labels.conf file
writeLines(ideogram_label_conf, args$ideogram_label_conf)





# Specify the lines to be written to the ticks_config file
ticks_spacing <- combined_variability_data[1, 3] - combined_variability_data[1, 2]

ticks_conf <- c(
  "show_ticks          = yes",
  "show_tick_labels    = yes",
  "",
  "<ticks>",
  "    radius           = 1r",
  "    color            = black",
  "    thickness        = 2p",
  "",
  "    # the tick label is derived by multiplying the tick position by 'multiplier' and casting it in 'format':",
  "    # sprintf(format,position*multiplier)",
  "    #",
  "",
  paste("    multiplier       = ", args$ticks_multiplier),
  "",
  "    # %d   - integer",
  "    # %f   - float",
  "    # %.1f - float with one decimal",
  "    # %.2f - float with two decimals",
  "    # for other formats, see http://perldoc.perl.org/functions/sprintf.html",
  "",
  "    format           = %d",
  "",
  "",
  "",
  "    <tick>",
  paste("        spacing        = ", ticks_spacing) ,
  "        size           = 15p",
  "        show_label     = yes",
  "        label_size     = 20p",
  "        label_offset   = 10p",
  "        format         = %d",
  "    </tick>",
  "",
  "</ticks>"
)

# Write the lines to a text file
writeLines(ticks_conf, args$ticks_conf)


# Specify the lines to be written to the file
max_variability = max(combined_variability_data[,4])
max_hypvar = max(combined_hypvar_data[,4])

main_conf <- c(
  "# karyotype file contains the list of chromosomes with their sizes for an organism in a given reference build",
  paste("karyotype = ",args$karyotype),
  "",
  "<image>",
  "    <<include etc/image.conf>> # included from Circos distribution",
  "</image>",
  "",
  "<plots>",
  "    # define the plot parameters for histogram containing nodes with deg > 2",
  "    <plot>",
  "        type    = histogram",
  paste("        file    = ", args$variability_output),
  "",
  "        # edit these min max values accordingly",
  "        min=0",
  paste("        max=",max_variability),
  "",
  "        # for one histogram plot in the circos",
  "        r0  = 0.25r",
  "        r1  = 0.80r",
  "",
  "        # Defining the background colour",
  "        <backgrounds>",
  "            <background>",
  "                color = vvlgrey",
  "            </background>",
  "        </backgrounds>",
  "",
  "        # colour of outline",
  "        color         = black_a4",
  "",
  "        # color of hisogram bars",
  paste("        fill_color = ",args$variability_color),
  "        thickness     = 4",
  "    </plot>",
  "",
  "     # define the plot parameters for histogram containing nodes with deg >= 5",
  "    <plot>",
  "        type    = histogram",
  paste("        file    = ", args$hyp_var_output),
  "",
  "        # edit these min max values accordingly",
  "        min=0",
  paste("        max=",max_hypvar),
  "",
  "        # for histogram plot in the circos",
  "        r0  = 0.82r",
  "        r1  = 0.98r",
  "        orientation = in",
  "",
  "        # Defining the background colour",
  "        <backgrounds>",
  "            <background>",
  "                color = vvlgrey",
  "            </background>",
  "        </backgrounds>",
  "",
  "        # colour of outline",
  "        color         = black_a4",
  "",
  "        # color of hisogram bars",
  paste("        fill_color = ", args$hypvar_color),
  "        thickness     = 4",
  "    </plot>",
  "",
  "</plots>",
  "",
  "#######################################################################################################",
  "#  Include other config files that support the plots",
  "",
  "# ideogram configuration",
  "<<include ideogram.conf>>",
  "",
  "# configuration of ticks in ideogram",
  "<<include ticks.conf>>",
  "",
  "# RGB/HSV color definitions, color lists, location of fonts, fill patterns",
  "<<include etc/colors_fonts_patterns.conf>> # included from Circos distribution",
  "",
  "# debugging, I/O and other system parameters",
  "<<include etc/housekeeping.conf>> # included from Circos distribution"
)

# Write the lines to a text file
writeLines(main_conf, args$main_conf)

cat( Sys.time(), " DONE: Configuration files for Circos plot collated", "\n")