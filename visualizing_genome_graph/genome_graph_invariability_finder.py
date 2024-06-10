import os
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


gfa_files_path = "/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files"
output_file_path = "/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/visualizing_genome_graph/gg_metrics/complete_1KGP_graph_invariant_zones.csv"

def create_DAG_from_GFA(file_path, ref_path_name):
    '''
    Function that reads a gfa file and returns a DiGraph networkx object
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
    # 
    # define a graph object
    genome_graph = nx.DiGraph()
    #
    # add all nodes and edges
    genome_graph.add_nodes_from(node_list)
    genome_graph.add_edges_from(edge_list)
    #
    del node_list, edge_list
    return genome_graph, reference_path
    del genome_graph, reference_path





def get_invariant_regions(genome_graph, ref_path, node_count_cutoff):
    '''
    Function that takes a networkx graph and returns the non-variant regions of the genome graph
    Parameters:
        genome_graph: the networkx object
        ref_path: name of the reference path in the GFA file
        node_count_cutoff: Min. no. of nodes to be present without variants to be considered an invariant region  
    '''
    i = 0
    invariant_zones = []
    while i < len(ref_path):
        if genome_graph.out_degree[ref_path[i]] == 1:
            current_invariant_zone = []
            # for the first node in ref_path
            if genome_graph.in_degree[ref_path[i]] == 0:
                #print("found 1st node at", ref_path[i])
                while True:
                    current_invariant_zone.append(ref_path[i])
                    if (genome_graph.out_degree[ref_path[i]] > 1) or (i == len(ref_path) - 1):
                        break
                    i = i + 1
                if len(current_invariant_zone) > node_count_cutoff:
                    invariant_zones.append(current_invariant_zone)
            #
            # for SNPs/deletion
            elif (genome_graph.out_degree[ref_path[i - 1]] > 1) and (genome_graph.in_degree[ref_path[i + 1]] > 1):
                #print("found snp/deletion at", ref_path[i])
                i = i + 1 #skip the SNP node
                if (genome_graph.out_degree[ref_path[i - 1]] == 1) and genome_graph.out_degree[ref_path[i]] == 1:
                    current_invariant_zone = []
                    while True:
                        current_invariant_zone.append(ref_path[i])
                        if (genome_graph.out_degree[ref_path[i]] > 1) or (i == len(ref_path) - 1):
                            break
                        i = i + 1
                else:
                    i = i + 1 # if adjacent variants are present, skip them too
                if len(current_invariant_zone) > node_count_cutoff:
                    invariant_zones.append(current_invariant_zone)
            # for nodes following deletions
            elif genome_graph.in_degree[ref_path[i]] > 1:
                #print("found post deletion node at", ref_path[i])
                while True:
                    current_invariant_zone.append(ref_path[i])
                    if (genome_graph.out_degree[ref_path[i]] > 1) or (i == len(ref_path) -1):
                        break
                    i = i + 1
                if len(current_invariant_zone) > node_count_cutoff:
                    invariant_zones.append(current_invariant_zone)
                del current_invariant_zone
        i = i + 1
    return invariant_zones
    del invariant_zones



# change working directory and get a list of all gfa files
os.chdir(gfa_files_path)
gfa_files = os.listdir()




# get a list to hold summary of invariant zones in [chr, length_of_zone, start_node, end_node] format 
invariant_zone_summary = []
for i in gfa_files:
    #get chromosome name from gfa file name: (if chr_name is the ref_path in gfa)
    #chr = "chrY"
    chr = i.split("_")[2].split(".")[0]
    # create genome graph from GFA file in lists
    genome_graph, ref_path = create_DAG_from_GFA(i, chr)
    # get invariant zones (no forks in 1000 nodes) in the genome graph
    invariant_zones = get_invariant_regions(genome_graph, ref_path, 1000)
    # write information about invariant zones to a list
    for j in invariant_zones:
        invariant_zone_summary.append([chr, len(j), j[0], j[-1]])
    print("Invariant zones captured for ", chr)
    del genome_graph, ref_path, invariant_zones

invariant_zones_in_genome_graph = pd.DataFrame(invariant_zone_summary, columns=["Chr", "Zone_length_in_nodes", "Zone_start_node", "Zone_end_node"])

invariant_zones_in_genome_graph.to_csv(output_file_path, index=False)


# # get invariant zones in toy data
# toy_invariant_zones = get_invariant_regions(sparse_gg, sparse_ref, 5)
# for i in toy_invariant_zones:
#      print(i)


# invariant_zones = get_invariant_regions(chr1_graph, chr1_ref_path, 1000)
# len(invariant_zones) 
# len_invariant_zones = []
# for i in invariant_zones:
#     len_invariant_zones.append(len(i))

# len_invariant_zones


