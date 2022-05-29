from partial_graph_match_s3 import Graph_Match
import timeit
import sys
import json



# find the node with most adjs
# node_kind : full or knockout
def find_bigest_set(match_graph, node_kind='full_don', nodes_to_ignore=[]):
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
    return max_node, max_set, max_count


#
def move_from_graph_to_set(match_graph, set_coverage, max_set):
    for node in max_set:
        set_coverage.add(node)
        match_graph.remove_node(node)


def call_greedy():
    with open('configuration_s3.json') as f:
        conf = json.load(f)

    create_cost = 1
    knocout_cost = conf["knockout_cost"] + create_cost
    pop_wanted_prec = conf["percentage_from_the_population"] /100 # the precentage of the population to cobvare
    start1 = timeit.default_timer()
    # knocout_cost = int(sys.argv[1])
    count_knockout = 0
    g_m = Graph_Match(0, 0)
    #sim = sys.argv[1]
    don_in_file = conf["donors_input_file"]#'../../data/israel/don_israel_5_geno_kir_' + sim + '.csv'
    rec_in_file = conf["recipients_input_file"]
    match_graph, _, _,_, pop_size = g_m.create_graph(don_in_file, rec_in_file)

    set_coverage = set()
    coverage_nodes = []
    number_of_coveraged = 0

    while number_of_coveraged / pop_size < pop_wanted_prec:
        node_knockout, max_set_knockout, covered_knockout_size = find_bigest_set(match_graph, 'knockout_don')

        set_size_without_knockout = 0
        list_of_nodes_that_cover = []
        max_set_full = set()
        while set_size_without_knockout < covered_knockout_size:
            node, max_set, covered_without_knockout_size = find_bigest_set(match_graph,
                                                                           nodes_to_ignore=list_of_nodes_that_cover)
            set_size_without_knockout += covered_without_knockout_size
            max_set_full = max_set_full.union(max_set)
            if len(max_set) == 0:
                break
            list_of_nodes_that_cover.append(node)
            # set_size_without_knockout += move_from_graph_to_set(match_graph, set_coverage, max_set)
        if knocout_cost / covered_knockout_size < (
                create_cost * len(list_of_nodes_that_cover)) / set_size_without_knockout:
            count_knockout += 1
            coverage_nodes.append(node_knockout)
            number_of_coveraged += covered_knockout_size
            move_from_graph_to_set(match_graph, set_coverage, max_set_knockout)
        else:
            coverage_nodes += list_of_nodes_that_cover
            move_from_graph_to_set(match_graph, set_coverage, max_set_full)
            number_of_coveraged += set_size_without_knockout


        print("#number of _nodes_that_cover:", len(coverage_nodes))
        print("#Number of covareted genotypes:", number_of_coveraged)
        print("#Percentage from the population: ", number_of_coveraged / pop_size * 100)
        print("#number of choose knockout:", count_knockout)

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
    print("time", timeit.default_timer() - start1)




if __name__ == '__main__':
    call_greedy()



