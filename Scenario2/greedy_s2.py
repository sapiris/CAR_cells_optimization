import cProfile
from partial_graph_match_s2 import Graph_Match
import json


def find_adjs(node, graph_match):
    list_homozygot = []
    locus_node = node.split('~')
    for i in range(0,len(locus_node), 2):
        if locus_node[i] == locus_node[i+1]:
            list_homozygot.append(locus_node[i])
    dict_i = {}
    for assitance_adj in graph_match.adj[node]:
        for full_adj in graph_match.adj[assitance_adj]:
            pat_with_homozygot_allele_of_donor = False
            #check if the pateint have allele thar it homozygot in the donor, if so dont includ him as adj
            for locus_homozygot in list_homozygot:
                if locus_homozygot in full_adj:
                    pat_with_homozygot_allele_of_donor = True
                    continue
            if not pat_with_homozygot_allele_of_donor:
                if full_adj in dict_i:
                    dict_i[full_adj] = max(dict_i[full_adj], len(assitance_adj.split('~')))
                else:
                    dict_i[full_adj] = len(assitance_adj.split('~'))

    return dict_i

#find the node with most adjs
def find_bigest_set(match_graph, prior_4):
    dict_id_size = {}
    max_node = ''
    max_set = set()
    max_count = 0
    for node in match_graph.nodes():
        node_description =  match_graph.nodes[node]
        if node_description["alleles"]=='full_don':
            node_count = 0
            dict_i = find_adjs(node, match_graph) #dict of adj and the number of matching alleles
            for adj, num_of_match in dict_i.items():
                adj_description = match_graph.nodes[adj]

                node_count += (adj_description["num_of_occurrence"] * prior_4)

            if node_count > max_count:
                max_count = node_count
                max_node = node
                max_set = dict_i
    return max_node, max_set
#
def move_from_graph_to_set(match_graph, set_coverage, max_set):
    weight = 0
    for node in max_set:
        set_coverage.add(node)
        weight += match_graph.nodes[node]["num_of_occurrence"]
        match_graph.remove_node(node)
    return weight

def update_match_level(max_set, dict_level_of_match, match_graph):
    for node, level_match in max_set.items():
        dict_level_of_match[level_match] +=  (1 *  match_graph.nodes[node]["num_of_occurrence"])

def call_greedy():
    with open('configuration_s2.json') as f:
        conf = json.load(f)

    pop_wanted_prec = conf["percentage_from_the_population"] /100 # the precentage of the population to cobvare
    don_in_file = conf["donors_input_file"]#'../../data/israel/don_israel_5_geno_kir_' + sim + '.csv'
    rec_in_file = conf["recipients_input_file"]
    prior_4 = int(conf.get("prior_for_4_match", 1))


    g_m = Graph_Match(int(conf["mismatch_a_b"]), 0, 'A*02:01', 'A*02:01')
    match_graph, _, _,_, pop_size = g_m.create_graph(don_in_file, rec_in_file)

    dict_level_of_match = {1:0, 2:0, 3:0, 4:0}
    set_coverage = set()
    coverage_nodes = []

    number_of_coveraged = 0

    while number_of_coveraged/pop_size < pop_wanted_prec:
        node, max_set = find_bigest_set(match_graph, prior_4)
        update_match_level(max_set, dict_level_of_match, match_graph)
        coverage_nodes.append(node)
        number_of_coveraged += move_from_graph_to_set(match_graph, set_coverage,max_set)
        print(f"##{len(coverage_nodes)},{number_of_coveraged/pop_size * 100}")

    f_out = open(conf["output_file"], "w")
    f_out.write(f"Number of nodes that cover: {len(coverage_nodes)}\n")
    f_out.write(f"Percentage from the population: {number_of_coveraged / pop_size * 100}\n")
    f_out.write(f"Percentage with match of 3: {dict_level_of_match[3] / pop_size * 100}\n"
                f"Percentage with match of 4: {dict_level_of_match[4] / pop_size * 100}\n")
    f_out.write("selected genotypes:\n")
    print('selected genotypes:')
    for node in coverage_nodes:
        print(node)
        f_out.write(node + '\n')

    print('selected genotypes:', len(coverage_nodes))
    print('Number of covareted genotypes:', number_of_coveraged )
    print("Percentage from the population: ", number_of_coveraged/pop_size * 100 )
    for match_level, count in dict_level_of_match.items():
        print(f"{match_level}, {count}")

if __name__ == '__main__':
    # Profiler start
    pr = cProfile.Profile()
    pr.enable()

    call_greedy()

    pr.disable()
    pr.print_stats(sort="time")
