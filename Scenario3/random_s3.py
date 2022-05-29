from partial_graph_match_s3 import Graph_Match
import timeit
import sys
import json

from random import randint, choice

# find the node with most adjs
# node_kind : full or knockout
"""def find_bigest_set(match_graph, node_kind='full_don', nodes_to_ignore=[]):
    dict_id_size = {}
    max_node = ''
    max_set = set()
    max_count = 0
    for node in match_graph.nodes():
        if match_graph.nodes[node]["alleles"] == node_kind:
            if node in nodes_to_ignore:
                continue
            node_count = 0
            set_i = set()
            for adj in match_graph.adj[node]:
                node_count += match_graph.nodes[adj]["num_of_occurrence"]
                set_i.add(adj)

            if node_count > max_count:
                max_count = node_count
                max_node = node
                max_set = set_i
    return max_node, max_set, max_count"""

def find_bigest_set(match_graph, dict_node_idx):
    node_count = 0
    set_i = set()
    node = choice(list(dict_node_idx.keys()))

    for adj in match_graph.adj[node]:
        node_count += match_graph.nodes[adj]["num_of_occurrence"]
        set_i.add(adj)
        to_continue = False

    return node, set_i, node_count


#
def move_from_graph_to_set(match_graph, set_coverage, max_set):
    for node in max_set:
        set_coverage.add(node)
        match_graph.remove_node(node)


def call_greedy():
    with open('configuration_s3.json') as f:
        conf = json.load(f)

    create_cost = 1
    #knocout_cost = conf["knockout_cost"] + create_cost
    pop_wanted_prec = 0.4# conf["percentage_from_the_population"] /100 # the precentage of the population to cobvare
    start1 = timeit.default_timer()
    # knocout_cost = int(sys.argv[1])
    count_knockout = 0
    g_m = Graph_Match(0, 0)
    #sim = sys.argv[1]
    don_in_file = conf["donors_input_file"]#'../../data/israel/don_israel_5_geno_kir_' + sim + '.csv'
    rec_in_file = conf["recipients_input_file"]
    match_graph, _, _,dict_node_idx, pop_size = g_m.create_graph(don_in_file, rec_in_file)

    set_coverage = set()
    coverage_nodes = []
    number_of_coveraged = 0
    count_full =0


    while number_of_coveraged/pop_size < pop_wanted_prec:

                node, max_set, count = find_bigest_set(match_graph,dict_node_idx)
                number_of_coveraged += count
                if len(node.split('~')) == 9:
                    count_knockout += 1
                else:
                    count_full += 1
                
                del dict_node_idx[node]



                print("#number of _nodes_that_cover:", len(coverage_nodes))
                print("#Number of covareted genotypes:", number_of_coveraged)
                print("#Percentage from the population: ", number_of_coveraged / pop_size * 100)
                print("#number of choose knockout:", count_knockout)
                print("cost:", count_knockout * 10 + count_full)

    f_out = open(conf["output_file"], "w")
    f_out.write(f"Number of _nodes_that_cover: {len(coverage_nodes)}\n")
    f_out.write(f"Percentage from the population: {number_of_coveraged / pop_size * 100}\n")
    f_out.write(f"Number of choose knockout: {count_knockout}\n")
    f_out.write("selected genotypes:\n")
    print('selected genotypes:')
    for node in coverage_nodes:
        print(node)
        f_out.write( node + '\n')
    print("number of _nodes_that_cover:", len(coverage_nodes))
    print('Number of covareted genotypes:', number_of_coveraged)
    print("Percentage from the population: ", number_of_coveraged / pop_size * 100)
    print("number of choose knockout:", count_knockout)
    print("cost:", count_knockout * 10 + count_full)
    print("time", timeit.default_timer() - start1)




if __name__ == '__main__':
    call_greedy()



