import os
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Get absolute paths for folder containing GFA files
gfa_files_path = "/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/10MB_subgraphs/gfa_files"

# Get absolute paths for files to write the outputs to
bin_summary_output_file = "/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/gg_metrics/10mb_summary_of_complete_1kgp_gg.csv"
bin_level_complexity_file = "/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/gg_metrics/bin_level_complete_1kgp_graph_complexity.txt"
bin_level_high_complexity_file = "/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/gg_metrics/bin_level_complete_1kgp_graph_highly_complexity.txt"
bin_level_max_degree_file = "/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/gg_metrics/bin_level_complete_1kgp_max_degree.txt"

# change working directory and get a list of all gfa files
os.chdir(gfa_files_path)
gfa_files = os.listdir()


# read gfa file
#f = open("x.gfa", "r")



def create_DAG_from_GFA(file_path, ref_path_name):
    '''
    Function that reads a gfa file and returns a list(DiGraph_networkx_object, nodes_in_reference_path)
    Parameters:
        file_path: relative/absolute path of the GFA file obtained from vg view
        ref_path_name: name of the reference path in the GFA file
    '''
    # define empty node and edges lists.
    node_list = []
    edge_list = []
    reference_path =[]
    f = open(file_path, "r")
    for line in f:
        # when encountering "S" lines in GFA extract the node ID and the sequence information. Push every node to node_list.
        if line.split("\t")[0] == 'S':
            node_id = line.split("\t")[1]
            node_sequence = line.split("\t")[2].split("\n")[0]
            node_object = (node_id, {"sequence" : node_sequence})
            node_list.append(node_object)
            del node_id, node_sequence, node_object
        # when encountering "L" lines in GFA extract the connecting nodes for each edge. Push every edge to edge_list.
        if line.split("\t")[0] == 'L':
            edge_left = line.split("\t")[1]
            edge_right = line.split("\t")[3]
            edge_object = (edge_left, edge_right)
            edge_list.append(edge_object)
            del edge_left, edge_right, edge_object
        # Get the reference path from GFA
        if line.split("\t")[0] == 'P':
            if line.split("\t")[1] == ref_path_name:
                reference_path = line.split("\t")[2].split("+,") #node ids are attached with + in ref path to indicate direction. So remove them.
                reference_path[len(reference_path) - 1] = reference_path[len(reference_path) -1].split("+")[0] #remove + in the last node.
 
    # define a graph object
    genome_graph = nx.DiGraph()

    # add all nodes and edges
    genome_graph.add_nodes_from(node_list)
    genome_graph.add_edges_from(edge_list)
    
    del node_list, edge_list
    return genome_graph, reference_path
    del genome_graph, reference_path



# create DAG from one GFA file directly
# genome_graph_Y = create_DAG_from_GFA('gi_200_graph_chrY.gfa')

chromosome_graphs = dict()
# create a dictionary with all DAGs from created from GFAs
for i in gfa_files:
    chr_path_name = i.split("_")[0]
    chromosome_graphs[i] = create_DAG_from_GFA(i, chr_path_name)
    print(i)

# analysing all the genome graphs
max(sorted((d for n, d in chromosome_graphs[gfa_files[1]][0].out_degree))) # Element 0 is graph and element 1 is reference path



# Get data to plot for individual bins

# initialize empty lists for each column. 
chromosome = []
start_position = []
end_position = []
max_deg = []
deg_1_count = []
deg_2_count = []
deg_3_count = []
deg_4_count = []
deg_5_count = []
deg_6_count = []
deg_7_count = []
deg_8_count = []
deg_9_count = []
deg_10_count = []
deg_11_count = []
deg_12_count = []
deg_13_count = []
deg_14_count = []
deg_15_count = []
deg_16_count = []
deg_17_count = []
deg_18_count = []
deg_19_count = []
deg_20_count = []

for i in gfa_files:
    out_degrees = sorted((d for n, d in chromosome_graphs[str(i)][0].out_degree), reverse=True) # 
    chromosome.append(i.split("_")[0])
    start_position.append(i.split("_")[1])
    end_position.append(i.split("_")[2].split(".")[0])
    #chromosome.append(i.split("_")[3].split(".")[0])
    deg_1_count.append(out_degrees.count(1))
    deg_2_count.append(out_degrees.count(2))
    deg_3_count.append(out_degrees.count(3))
    deg_4_count.append(out_degrees.count(4))
    deg_5_count.append(out_degrees.count(5))
    deg_6_count.append(out_degrees.count(6))
    deg_7_count.append(out_degrees.count(7))
    deg_8_count.append(out_degrees.count(8))
    deg_9_count.append(out_degrees.count(9))
    deg_10_count.append(out_degrees.count(10))
    deg_11_count.append(out_degrees.count(11))
    deg_12_count.append(out_degrees.count(12))
    deg_13_count.append(out_degrees.count(13))
    deg_14_count.append(out_degrees.count(14))
    deg_15_count.append(out_degrees.count(15))
    deg_16_count.append(out_degrees.count(16))
    deg_17_count.append(out_degrees.count(17))
    deg_18_count.append(out_degrees.count(18))
    deg_19_count.append(out_degrees.count(19))
    deg_20_count.append(out_degrees.count(20))
    max_deg.append(max(out_degrees))
    print(i)

bin_summary = pd.DataFrame(list(zip(chromosome, start_position, end_position, max_deg, 
deg_1_count, deg_2_count, deg_3_count, deg_4_count, deg_5_count, deg_6_count, deg_7_count, deg_8_count, 
deg_9_count, deg_10_count, deg_11_count, deg_12_count, deg_13_count, deg_14_count, deg_15_count, 
deg_16_count, deg_17_count, deg_18_count, deg_19_count, deg_20_count)),
columns= ["chromosome", "start_position", "end_position", "max_deg", 
"deg_1_count", "deg_2_count", "deg_3_count", "deg_4_count", "deg_5_count", "deg_6_count", "deg_7_count", "deg_8_count", 
"deg_9_count", "deg_10_count", "deg_11_count", "deg_12_count", "deg_13_count", "deg_14_count", "deg_15_count", 
"deg_16_count", "deg_17_count", "deg_18_count", "deg_19_count", "deg_20_count"])

bin_summary.to_csv(bin_summary_output_file, index=False, header=True)

# chr_summary = pd.DataFrame(list(zip(chromosome, start_position, end_position, max_deg, 
# deg_1_count, deg_2_count, deg_3_count, deg_4_count, deg_5_count, deg_6_count, deg_7_count, deg_8_count, 
# deg_9_count, deg_10_count, deg_11_count, deg_12_count, deg_13_count, deg_14_count, deg_15_count)),
# columns= ["chromosome", "start_position", "end_position", "max_deg", 
# "deg_1_count", "deg_2_count", "deg_3_count", "deg_4_count", "deg_5_count", "deg_6_count", "deg_7_count", "deg_8_count", 
# "deg_9_count", "deg_10_count", "deg_11_count", "deg_12_count", "deg_13_count", "deg_14_count", "deg_15_count"])

# chr_summary.to_csv("/cn2/data2/venkatesh_data2/visualizing_gg/gg_metrics/chr_summary_of_gg.csv", index=False, header=True)

max(bin_summary["max_deg"])



# Analysis of bin_summary file

binned_complexity = pd.DataFrame()
binned_complexity["chromosome"] = bin_summary["chromosome"]
binned_complexity["start_position"] = bin_summary["start_position"].astype(int) - 1
binned_complexity["end_position"] = bin_summary["end_position"].astype(int) - 1
binned_complexity["total_complex_nodes"] = bin_summary["deg_2_count"].astype(int) + bin_summary["deg_3_count"].astype(int) + bin_summary["deg_4_count"].astype(int) + bin_summary["deg_5_count"].astype(int) + bin_summary["deg_6_count"].astype(int) + bin_summary["deg_7_count"].astype(int) + bin_summary["deg_8_count"].astype(int) + bin_summary["deg_9_count"].astype(int) + bin_summary["deg_10_count"].astype(int) + bin_summary["deg_11_count"].astype(int) + bin_summary["deg_12_count"].astype(int) + bin_summary["deg_13_count"].astype(int) + bin_summary["deg_14_count"].astype(int) + bin_summary["deg_15_count"].astype(int) + bin_summary["deg_16_count"].astype(int) + bin_summary["deg_17_count"].astype(int) + bin_summary["deg_18_count"].astype(int) + bin_summary["deg_19_count"].astype(int) + bin_summary["deg_20_count"].astype(int)
binned_complexity.to_csv(bin_level_complexity_file, index=None, header=None, sep='\t')





bin_max_degree = pd.DataFrame()
bin_max_degree["chromosome"] = bin_summary["chromosome"]
bin_max_degree["start_position"] = bin_summary["start_position"].astype(int) - 1
bin_max_degree["end_position"] = bin_summary["end_position"].astype(int) - 1
bin_max_degree["max_degree"] = bin_summary["max_deg"]
bin_max_degree.to_csv(bin_level_max_degree_file, index=None, header=None, sep='\t')



# read bin summary output file 
# bin_summary = pd.read_csv(bin_summary_output_file)


# create binned summary for nodes with degree >= 5
# bin_summary = pd.read_csv("/cn4/data4/venkatesh_data4/gi_cbr_202_genome_graphs/vg_manual_construct/visualizing_genome_graphs/gg_metrics/10mb_summary_of_gg.csv")
binned_high_complexity = pd.DataFrame()
binned_high_complexity["chromosome"] = bin_summary["chromosome"]
binned_high_complexity["start_position"] = bin_summary["start_position"].astype(int) - 1
binned_high_complexity["end_position"] = bin_summary["end_position"].astype(int) - 1
binned_high_complexity["highly_complex_nodes"] = bin_summary["deg_5_count"].astype(int) + bin_summary["deg_6_count"].astype(int) + bin_summary["deg_7_count"].astype(int) + bin_summary["deg_8_count"].astype(int) + bin_summary["deg_9_count"].astype(int) + bin_summary["deg_10_count"].astype(int) + bin_summary["deg_11_count"].astype(int) + bin_summary["deg_12_count"].astype(int) + bin_summary["deg_13_count"].astype(int) + bin_summary["deg_14_count"].astype(int) + bin_summary["deg_15_count"].astype(int) + bin_summary["deg_16_count"].astype(int) + bin_summary["deg_17_count"].astype(int) + bin_summary["deg_18_count"].astype(int) + bin_summary["deg_19_count"].astype(int) + bin_summary["deg_20_count"].astype(int)
binned_high_complexity.to_csv(bin_level_high_complexity_file, index=None, header=None, sep='\t')


# create binned summary for nodes with degree == 5
# # bin_summary = pd.read_csv("/cn4/data4/venkatesh_data4/gi_cbr_202_genome_graphs/vg_manual_construct/visualizing_genome_graphs/gg_metrics/10mb_summary_of_gg.csv")
# binned_degree = pd.DataFrame()
# binned_degree["chromosome"] = bin_summary["chromosome"]
# binned_degree["start_position"] = bin_summary["start_position"].astype(int) - 1
# binned_degree["end_position"] = bin_summary["end_position"].astype(int) - 1
# binned_degree["highly_complex_nodes"] = bin_summary["deg_5_count"].astype(int)
# binned_degree.to_csv("/cn4/data4/venkatesh_data4/gi_cbr_824_genome_graphs/vg_manualconstruct/visualizing_genome_graphs/gg_metrics/bin_level_graph_degree_5.txt", index=None, header=None, sep='\t')

