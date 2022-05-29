import cProfile
from partial_graph_match_s1 import Graph_Match
import json


def find_adjs(node, graph_match):
    set_i = set()
    for assitance_adj in graph_match.adj[node]:
        for full_adj in graph_match.adj[assitance_adj]:

            set_i.add(full_adj)

    return set_i


def find_bigest_set(match_graph):
    max_node = ''
    max_set = set()
    max_count = 0
    for node in match_graph.nodes():
        node_description =  match_graph.nodes[node]
        if node_description["alleles"]=='full_don':
            node_kir = node_description["kir"]
            node_count = 0
            set_i_tmp = find_adjs(node, match_graph)
            set_i = set()
            for adj in set_i_tmp:
                adj_description = match_graph.nodes[adj]
                if adj_description["kir"] != node_kir:
                    set_i.add(adj)
                    node_count += adj_description["num_of_occurrence"]

            if node_count > max_count:
                max_count = node_count
                max_node = node
                max_set = set_i
    return max_node, max_set
#
def move_from_graph_to_set(match_graph, max_set):
    weight = 0
    for node in max_set:
        weight += match_graph.nodes[node]["num_of_occurrence"]
        match_graph.remove_node(node)
    return weight


def call_greedy():

    with open('configuration_s1.json') as f:
        conf = json.load(f)

    pop_wanted_prec = conf["percentage_from_the_population"] /100 # the precentage of the population to cobvare
    don_in_file = conf["donors_input_file"]#'../../data/israel/don_israel_5_geno_kir_' + sim + '.csv'
    rec_in_file = conf["recipients_input_file"]


    g_m = Graph_Match(6, 0)
    match_graph, _, _,_, pop_size = g_m.create_graph(don_in_file, rec_in_file)#israel/israel_5_geno_kir.csv

    coverage_nodes = []

    number_of_coveraged = 0
    while number_of_coveraged/pop_size < pop_wanted_prec:
        node, max_set = find_bigest_set(match_graph)
        coverage_nodes.append(node)
        number_of_coveraged += move_from_graph_to_set(match_graph, max_set)
        print(f"##{len(coverage_nodes)},{number_of_coveraged/pop_size * 100}")


    print('Number of covareted genotypes:', number_of_coveraged )
    print("Percentage from the population: ", number_of_coveraged/pop_size * 100 )

    f_out = open(conf["output_file"], "w")
    f_out.write(f"Number of nodes that cover: {len(coverage_nodes)}\n")
    f_out.write(f"Percentage from the population: {number_of_coveraged / pop_size * 100}\n")
    f_out.write("selected genotypes:\n")
    print('selected genotypes:')
    for node in coverage_nodes:
        print(node)
        f_out.write(node + '\n')



if __name__ == '__main__':
    # Profiler start
    pr = cProfile.Profile()
    pr.enable()

    call_greedy()

    pr.disable()
    pr.print_stats(sort="time")
