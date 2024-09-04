'''
Python module to perform structural analysis of genome graphs: 
   1) Creates networkx version of genome graphs
   2) Builds a Node_ID to Genomic_location lookup table
   3) Identify Invariant zones in the genome graph
   4) Breaks the genome graph into binned subgraphs
   5) Measure the variability and hypervariability (degree-based complexities) of each of the subgraph 
   6) Output the relevant files

Author: Venkatesh Kamaraj
Contact: ic11570@imail.iitm.ac.in
'''


import networkx as nx
import pandas as pd
import argparse
from datetime import datetime


# Parse arguments from the command line
parser = argparse.ArgumentParser(description='A python module to perform structural analysis of genome graphs')
parser.add_argument("--gfa_file", help="Path to the genome graphs in GFA format")
parser.add_argument("--ref_path", help="Name of the path tracing the reference genome in the GFA file")
parser.add_argument("--bin_width", help="Width of a bin in bp - the reference path broken into bins and subgraphs made for each bin", type=int)
parser.add_argument("--invariance_cutoff", help="Minimum number of continious nodes with no variability to be termed as Invariant Zone", type=int)
parser.add_argument("--save_lookup", choices=['yes', 'no'], help="Choice to save the [Genomic_location, Node_id] lookup table")
parser.add_argument("--lookup_output", help="Path to the save the lookup table in .csv format")
parser.add_argument("--degree_output", help="Path to the save the detailed count of binned nodes with each out_degre (max 30) in .csv format")
parser.add_argument("--var_output", help="Path to the save the binned variability counts in .csv format")
parser.add_argument("--hypvar_output", help="Path to the save the binned hypervariability in .csv format")
parser.add_argument("--invar_output", help="Path to the save the invariant zones in .csv format")
args = parser.parse_args()



# Create NetworkX object of the GFA file 
def create_DAG_from_GFA(file_path, reference_path_name):
    '''
    Function that reads a gfa file and returns a DiGraph networkx object
    Parameters:
        file_path: relative/absolute path of the GFA file obtained from vg view
        reference_path_name: name of the reference path in the GFA file
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
            if line.split("\t")[1] == reference_path_name:
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



# Get lookup table for all the gfa file that matches Node ID with Genomic location
def create_node_reference_lookup(genome_graph, ref_path_nodes):
    '''
    Function that takes a networkx graph and returns the non-variant regions of the genome graph
    Parameters:
        genome_graph: the networkx object
        ref_path_nodes: list of nodes in the reference path in the GFA file 
    '''
    # Step 1: get the amount of nucleotides captured by each node in reference path
    cumulative_node_length = [0]     # starting cumulative distance from 0
    for i in ref_path_nodes:
        # add length of current node to the cumulative length at previous node
        cumulative_node_length.append(cumulative_node_length[-1] + len(genome_graph.nodes(data="sequence")[i]))
    del cumulative_node_length[0] # remove the first entry initialized with 0
    
    # Step 2: When reference location is less than a cumulative distance, it'll be contained in the latest node 
    ref_location = []
    ref_node = []
    flag = 0
    for i in range(1,cumulative_node_length[-1]+1,1):
        if (i > cumulative_node_length[flag]) and (flag < len(ref_path_nodes) - 1):
            flag = flag + 1
        if i <= cumulative_node_length[flag]:
            ref_location.append(i)
            ref_node.append(ref_path_nodes[flag])

    # Step 3: Make a lookup table for reference location and its corresponding node
    ref_lookup_table = pd.DataFrame()
    ref_lookup_table["Reference_location"] = ref_location
    ref_lookup_table["Genome_graph_node"] = ref_node
    
    # delete used lists to free space
    ## (cumulative_node_length)
    del cumulative_node_length
    del ref_location, ref_node

    return  ref_lookup_table
    del ref_lookup_table



def get_invariant_regions(genome_graph, ref_path_nodes, node_count_cutoff):
    '''
    Function that takes a networkx graph and returns the non-variant regions of the genome graph
    Parameters:
        genome_graph: the networkx object
        ref_path_nodes: list of nodes in the reference path in the GFA file
        node_count_cutoff: Min. no. of nodes to be present without variants to be considered an invariant region  
    '''
    i = 0
    invariant_zones = []
    while i < len(ref_path_nodes):
        if genome_graph.out_degree[ref_path_nodes[i]] == 1:
            current_invariant_zone = []
            # for the first node in ref_path_nodes
            if genome_graph.in_degree[ref_path_nodes[i]] == 0:
                #print("found 1st node at", ref_path_nodes[i])
                while True:
                    current_invariant_zone.append(ref_path_nodes[i])
                    if (genome_graph.out_degree[ref_path_nodes[i]] > 1) or (i == len(ref_path_nodes) - 1):
                        break
                    i = i + 1
                if len(current_invariant_zone) > node_count_cutoff:
                    invariant_zones.append(current_invariant_zone)
            #
            # for SNPs/deletion
            elif (genome_graph.out_degree[ref_path_nodes[i - 1]] > 1) and (genome_graph.in_degree[ref_path_nodes[i + 1]] > 1):
                #print("found snp/deletion at", ref_path_nodes[i])
                i = i + 1 #skip the SNP node
                if (genome_graph.out_degree[ref_path_nodes[i - 1]] == 1) and genome_graph.out_degree[ref_path_nodes[i]] == 1:
                    current_invariant_zone = []
                    while True:
                        current_invariant_zone.append(ref_path_nodes[i])
                        if (genome_graph.out_degree[ref_path_nodes[i]] > 1) or (i == len(ref_path_nodes) - 1):
                            break
                        i = i + 1
                else:
                    i = i + 1 # if adjacent variants are present, skip them too
                if len(current_invariant_zone) > node_count_cutoff:
                    invariant_zones.append(current_invariant_zone)
            # for nodes following deletions
            elif genome_graph.in_degree[ref_path_nodes[i]] > 1:
                #print("found post deletion node at", ref_path_nodes[i])
                while True:
                    current_invariant_zone.append(ref_path_nodes[i])
                    if (genome_graph.out_degree[ref_path_nodes[i]] > 1) or (i == len(ref_path_nodes) -1):
                        break
                    i = i + 1
                if len(current_invariant_zone) > node_count_cutoff:
                    invariant_zones.append(current_invariant_zone)
                del current_invariant_zone
        i = i + 1
    return invariant_zones
    del invariant_zones



# Create genome graph from the gfa file
genome_graph, nodes_in_ref_path = create_DAG_from_GFA(args.gfa_file, args.ref_path)

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "DONE: NetworkX object created for the GFA file: ", args.gfa_file)


# get the node-ref_location information
gg_ref_lookup = create_node_reference_lookup(genome_graph, nodes_in_ref_path)

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "DONE: Lookup table created for the GFA file: ", args.gfa_file)

# Save the lookup table if the user asks for it
if args.save_lookup == 'yes':
    gg_ref_lookup.to_csv(args.lookup_output, index = False)

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "DONE: Lookup table saved to ", args.lookup_output, "for the GFA file: ", args.gfa_file)


# Identify Invariant zones from the GFA file
invariant_zones = get_invariant_regions(genome_graph, nodes_in_ref_path, args.invariance_cutoff)

# Write information about invariant zones to a list
invariant_zone_summary = []
for j in invariant_zones:
    invariant_zone_summary.append([args.ref_path, len(j), j[0], j[-1]])
del invariant_zones

# Make dataframe of the invariant zones list
invariant_zones_in_genome_graph = pd.DataFrame(invariant_zone_summary, columns=["Chr", "Zone_length_in_nodes", "Zone_start_node", "Zone_end_node"])

invariant_zones_in_genome_graph.to_csv(args.invar_output, index=False)

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "DONE: Invariant zones extracted for the GFA file: ", args.gfa_file)


# Get the reference node-loc lookup table. Headers are predefined as [Reference_location, Genome_graph_node]
# gg_ref_lookup = pd.read_csv("ref_lookup.csv") 


# Get the start node and end node for subgraphs of each 20kbp genomic window
subgraph_start_nodes = []
subgraph_end_nodes = []
subgraph_start_loc = []
subgraph_end_loc = []
for ref_loc in range(1,len(gg_ref_lookup["Reference_location"]),args.bin_width):
    subgraph_start_nodes.append(gg_ref_lookup["Genome_graph_node"][ref_loc - 1]) # pandas indexes from 0 and ref_loc starts from 1. So subtract to match
    subgraph_start_loc.append(gg_ref_lookup["Reference_location"][ref_loc - 1])
    #untill the last bin is encountered, update the end node
    if ref_loc < len(gg_ref_lookup["Reference_location"]) - args.bin_width:
        subgraph_end_nodes.append(gg_ref_lookup["Genome_graph_node"][ref_loc + args.bin_width - 2]) # start from 1 and idexed 1. So subtract 2
        subgraph_end_loc.append(gg_ref_lookup["Reference_location"][ref_loc + args.bin_width - 2])
    # for the last bin, end node is the last node
    if ref_loc > len(gg_ref_lookup["Reference_location"]) - args.bin_width:
        subgraph_end_nodes.append(gg_ref_lookup["Genome_graph_node"][len(gg_ref_lookup["Reference_location"]) - 1])
        subgraph_end_loc.append(gg_ref_lookup["Reference_location"][len(gg_ref_lookup["Reference_location"]) - 1])


# get subgraphs for each 20kbp genomic window (subgraph command requires all the node ids in the subgraph)
# create one subgraph
# subgraph = genome_graph.subgraph(map(str, range(int(subgraph_start_nodes[1]), int(subgraph_end_nodes[1])+1))) # range omits the last number. So add 1.
# subgraph_2 = genome_graph.subgraph(map(str, range(1, 1056))) # range omits the last number. So add 1.


# Create a dictionary of subgraphs
binned_subgraphs = dict()
# create a dictionary with all DAGs from created from GFAs
for i in range(len(subgraph_start_nodes)):
    binned_subgraphs[str(i)] = genome_graph.subgraph(map(str, range(int(subgraph_start_nodes[i]), int(subgraph_end_nodes[i])+1))) # range omits the last number. So add 1.

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "DONE: Subgraphs created for the contig: ", args.ref_path)


# Get summary data for individual bins

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
deg_21_count = []
deg_22_count = []
deg_23_count = []
deg_24_count = []
deg_25_count = []
deg_26_count = []
deg_27_count = []
deg_28_count = []
deg_29_count = []
deg_30_count = []

for i in range(len(subgraph_start_nodes)):
    out_degrees = sorted((d for n, d in binned_subgraphs[str(i)].out_degree), reverse=True) # 
    start_position.append(subgraph_start_loc[i])
    end_position.append(subgraph_end_loc[i])
    chromosome.append(args.ref_path)
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
    deg_21_count.append(out_degrees.count(21))
    deg_22_count.append(out_degrees.count(22))
    deg_23_count.append(out_degrees.count(23))
    deg_24_count.append(out_degrees.count(24))
    deg_25_count.append(out_degrees.count(25))
    deg_26_count.append(out_degrees.count(26))
    deg_27_count.append(out_degrees.count(27))
    deg_28_count.append(out_degrees.count(28))
    deg_29_count.append(out_degrees.count(29))
    deg_30_count.append(out_degrees.count(30))
    max_deg.append(max(out_degrees))
    # print(i)

bin_summary = pd.DataFrame(list(zip(chromosome, start_position, end_position, max_deg, 
deg_1_count, deg_2_count, deg_3_count, deg_4_count, deg_5_count, deg_6_count, deg_7_count, deg_8_count, 
deg_9_count, deg_10_count, deg_11_count, deg_12_count, deg_13_count, deg_14_count, deg_15_count, 
deg_16_count, deg_17_count, deg_18_count, deg_19_count, deg_20_count, deg_21_count, 
deg_22_count, deg_23_count, deg_24_count, deg_25_count, deg_26_count, deg_27_count, deg_28_count, deg_29_count, deg_30_count)),
columns= ["chromosome", "start_position", "end_position", "max_deg", 
"deg_1_count", "deg_2_count", "deg_3_count", "deg_4_count", "deg_5_count", "deg_6_count", "deg_7_count", "deg_8_count", 
"deg_9_count", "deg_10_count", "deg_11_count", "deg_12_count", "deg_13_count", "deg_14_count", "deg_15_count", 
"deg_16_count", "deg_17_count", "deg_18_count", "deg_19_count", "deg_20_count", "deg_21_count", "deg_22_count", 
"deg_23_count", "deg_24_count", "deg_25_count", "deg_26_count", "deg_27_count", "deg_28_count", "deg_29_count", "deg_30_count"])

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "REPORT: Maximum degree observed in ", args.ref_path, "is ", max(bin_summary["max_deg"]))

bin_summary.to_csv(args.degree_output, index=False, header=True)

# chr_summary = pd.DataFrame(list(zip(chromosome, start_position, end_position, max_deg, 
# deg_1_count, deg_2_count, deg_3_count, deg_4_count, deg_5_count, deg_6_count, deg_7_count, deg_8_count, 
# deg_9_count, deg_10_count, deg_11_count, deg_12_count, deg_13_count, deg_14_count, deg_15_count)),
# columns= ["chromosome", "start_position", "end_position", "max_deg", 
# "deg_1_count", "deg_2_count", "deg_3_count", "deg_4_count", "deg_5_count", "deg_6_count", "deg_7_count", "deg_8_count", 
# "deg_9_count", "deg_10_count", "deg_11_count", "deg_12_count", "deg_13_count", "deg_14_count", "deg_15_count"])

# chr_summary.to_csv("chr_summary_of_gg.csv", index=False, header=True)





# Analysis of bin_summary file
# create binned summary for variability: nodes with degree >= 2
binned_variability = pd.DataFrame()
binned_variability["chromosome"] = bin_summary["chromosome"]
binned_variability["start_position"] = bin_summary["start_position"].astype(int) - 1
binned_variability["end_position"] = bin_summary["end_position"].astype(int) # Make the bins overlapping for compatibility with Circos
binned_variability["total_variable_nodes"] = bin_summary["deg_2_count"].astype(int) + bin_summary["deg_3_count"].astype(int) + bin_summary["deg_4_count"].astype(int) + bin_summary["deg_5_count"].astype(int) + bin_summary["deg_6_count"].astype(int) + bin_summary["deg_7_count"].astype(int) + bin_summary["deg_8_count"].astype(int) + bin_summary["deg_9_count"].astype(int) + bin_summary["deg_10_count"].astype(int) + bin_summary["deg_11_count"].astype(int) + bin_summary["deg_12_count"].astype(int) + bin_summary["deg_13_count"].astype(int) + bin_summary["deg_14_count"].astype(int) + bin_summary["deg_15_count"].astype(int) + bin_summary["deg_16_count"].astype(int) + bin_summary["deg_17_count"].astype(int) + bin_summary["deg_18_count"].astype(int) + bin_summary["deg_19_count"].astype(int) + bin_summary["deg_20_count"].astype(int) + bin_summary["deg_21_count"].astype(int) + bin_summary["deg_22_count"].astype(int) + bin_summary["deg_23_count"].astype(int) + bin_summary["deg_24_count"].astype(int) + bin_summary["deg_25_count"].astype(int) + bin_summary["deg_26_count"].astype(int) + bin_summary["deg_27_count"].astype(int) + bin_summary["deg_28_count"].astype(int) + bin_summary["deg_29_count"].astype(int) + bin_summary["deg_30_count"].astype(int)

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "REPORT: Maximum variability observed in ", args.ref_path, "is ", max(binned_variability["total_variable_nodes"]))

binned_variability.to_csv(args.var_output, index=None, header=None, sep='\t')



# create binned summary for hypervariability: nodes with degree >= 5
# bin_summary = pd.read_csv("summary_of_gg.csv")
binned_hypervariability = pd.DataFrame()
binned_hypervariability["chromosome"] = bin_summary["chromosome"]
binned_hypervariability["start_position"] = bin_summary["start_position"].astype(int) - 1
binned_hypervariability["end_position"] = bin_summary["end_position"].astype(int) # Make the bins overlapping for compatibility with Circos
binned_hypervariability["total_hypervariable_nodes"] = bin_summary["deg_5_count"].astype(int) + bin_summary["deg_6_count"].astype(int) + bin_summary["deg_7_count"].astype(int) + bin_summary["deg_8_count"].astype(int) + bin_summary["deg_9_count"].astype(int) + bin_summary["deg_10_count"].astype(int) + bin_summary["deg_11_count"].astype(int) + bin_summary["deg_12_count"].astype(int) + bin_summary["deg_13_count"].astype(int) + bin_summary["deg_14_count"].astype(int) + bin_summary["deg_15_count"].astype(int) + bin_summary["deg_16_count"].astype(int) + bin_summary["deg_17_count"].astype(int) + bin_summary["deg_18_count"].astype(int) + bin_summary["deg_19_count"].astype(int) + bin_summary["deg_20_count"].astype(int) + bin_summary["deg_21_count"].astype(int) + bin_summary["deg_22_count"].astype(int) + bin_summary["deg_23_count"].astype(int) + bin_summary["deg_24_count"].astype(int) + bin_summary["deg_25_count"].astype(int) + bin_summary["deg_26_count"].astype(int) + bin_summary["deg_27_count"].astype(int) + bin_summary["deg_28_count"].astype(int) + bin_summary["deg_29_count"].astype(int) + bin_summary["deg_30_count"].astype(int)

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "REPORT: Maximum hypervariability observed in ", args.ref_path, "is ", max(binned_hypervariability["total_hypervariable_nodes"]))

binned_hypervariability.to_csv(args.hypvar_output, index=None, header=None, sep='\t')


# bin_max_degree = pd.DataFrame()
# bin_max_degree["chromosome"] = bin_summary["chromosome"]
# bin_max_degree["start_position"] = bin_summary["start_position"].astype(int) - 1
# bin_max_degree["end_position"] = bin_summary["end_position"].astype(int) - 1
# bin_max_degree["max_degree"] = bin_summary["max_deg"]
# bin_level_max_degree_file = output_directory + "/" + "graph_" + "10mb_bin_max_degree.csv"
# bin_max_degree.to_csv(bin_level_max_degree_file, index=None, header=None, sep='\t')


print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "DONE: Structural analysis completed for the Genome Graph contained in: ", args.gfa_file)