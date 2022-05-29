from partial_graph_match_s2 import Graph_Match
import pulp
import cProfile
import json

def find_adjs(node, graph_match):

    dict_i = {}
    for assitance_adj in graph_match.adj[node]:
        for full_adj in graph_match.adj[assitance_adj]:
            if not "pat" in full_adj:
                pat_with_homozygot_allele_of_donor = False
                #check if the pateint have allele thar it homozygot in the donor, if so dont includ him as adj
                locus_adj = full_adj.split('~')
                list_homozygot = []
                for i in range(0, len(locus_adj), 2):
                    if locus_adj[i] == locus_adj[i + 1]:
                        list_homozygot.append(locus_adj[i])
                for locus_homozygot in list_homozygot:
                    if locus_homozygot in node:
                        pat_with_homozygot_allele_of_donor = True
                        continue
                if not pat_with_homozygot_allele_of_donor:
                    if full_adj in dict_i:
                        dict_i[full_adj] = max(dict_i[full_adj], len(assitance_adj.split('~')))
                    else:
                        dict_i[full_adj] = len(assitance_adj.split('~'))

    return dict_i

def init_lp_prolem(match_graph, pop_size, dict_idx_node_don, dict_idx_node_pat , dict_node_idx, lim1):

    lp_problem = pulp.LpProblem("LP_Problem", pulp.LpMinimize)

    ##variables boundaries
    var = []
    for i in range(len(dict_idx_node_don)):
        var.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Integer'))
    for i in range(len(dict_idx_node_pat)):
        var.append(pulp.LpVariable('y' + str(i), lowBound=0, upBound=1, cat='Integer'))

    ##equation to minimize
    z = 0
    for idx, node in dict_idx_node_don.items():
        z += var[idx]
    lp_problem += z

    ##Constraints
    #(sum(x_i) <= k
    k = pop_size * lim1

    sum_x_i = -1 * var[len(dict_idx_node_don)] * match_graph.nodes[dict_idx_node_pat[0]]["num_of_occurrence"]
    for i in range(1 + len(dict_idx_node_don), len(dict_idx_node_don) + len(dict_idx_node_pat)):
        sum_x_i -= (var[i] * match_graph.nodes[dict_idx_node_pat[i - len(dict_idx_node_don)]]["num_of_occurrence"])
    lp_problem += sum_x_i <= (-1 * k)

    sum_x_i = var[len(dict_idx_node_don)] * match_graph.nodes[dict_idx_node_pat[0]]["num_of_occurrence"]
    for i in range(1 + len(dict_idx_node_don), len(dict_idx_node_don) + len(dict_idx_node_pat)):
        sum_x_i += (var[i] * match_graph.nodes[dict_idx_node_pat[i - len(dict_idx_node_don)]]["num_of_occurrence"])
    lp_problem += sum_x_i <= (k)



    for i in range(len(dict_idx_node_pat)):
        cons_vec = 0
        node = dict_idx_node_pat[i]
        #kir_node = match_graph.nodes[node]["kir"]
        set_i = find_adjs(node, match_graph)
        for adj in set_i:
            #if (adj != node and match_graph.nodes[adj]["kir"] != kir_node) :
                    cons_vec += (-1 * var[dict_node_idx[adj]])
        cons_vec += var[len(dict_idx_node_don) + i]

        lp_problem += cons_vec <= 0



    print(lp_problem)

    return lp_problem, var







def check_simlarity(geno1, geno2):

    count_similar = 0

    for i in range(0, 4, 2):
        if geno1[i] == geno2[i]:
            if geno1[i + 1] == geno2[i + 1]:
                count_similar += 2
            else:
                count_similar += 1
        elif geno1[i + 1] == geno2[i + 1]:
            count_similar += 1
        elif geno1[i] == geno2[i + 1]:
            count_similar += 1
        elif geno1[i + 1] == geno2[i]:
            count_similar += 1

    return count_similar


if __name__ == '__main__':
    pr = cProfile.Profile()
    pr.enable()

    with open('configuration_s2.json') as f:
        conf = json.load(f)

    pop_wanted_prec = conf["percentage_from_the_population"] /100 # the precentage of the population to cobvare
    don_in_file = conf["donors_input_file"]#'../../data/israel/don_israel_5_geno_kir_' + sim + '.csv'
    rec_in_file = conf["recipients_input_file"]


    f_out = open(conf["output_file"], "w")


    g_m = Graph_Match(int(conf["mismatch_a_b"]), 0, 'A*02:01','A*02:01', False)
    match_graph, dict_idx_node_don, dict_idx_node_pat, dict_node_idx, pop_size =  g_m.create_graph(don_in_file, rec_in_file)

    lp_problem, var = init_lp_prolem(match_graph, pop_size, dict_idx_node_don, dict_idx_node_pat, dict_node_idx, pop_wanted_prec)

    lp_problem.solve()

    ll = lp_problem.variables()


    num_of_geno_that_covarge = 0
    res_value = pulp.value(lp_problem.objective)

    dict_y_res = {}
    dict_num_matchs = {3:0.0, 4:0.0}
    dict_geno_covaregs = {}
    f_out.write("selected genotypes:\n")
    for i, variable in enumerate(lp_problem.variables()):
        if "x" in variable.name:
            index = int(variable.name.replace('x', ''))
            if variable.varValue == 1:
                f_out.write(dict_idx_node_don[index] + '\n')
                dict_geno_covaregs[dict_idx_node_don[index]] = variable.varValue
        if 'y' in variable.name:
            index = int(variable.name.replace('y', ''))
            if variable.varValue > 0:
                geno_y = dict_idx_node_pat[index]
                max_match = 0
                geno_y = geno_y.replace('pat-','').split('~')[:4]
                max_x = 0
                for geno_x in dict_geno_covaregs:
                    genox = geno_x.split('~')[:4]
                    match = check_simlarity(geno_y, genox)
                    if match > max_match:
                        max_match = match
                        max_x = dict_geno_covaregs[geno_x]
                    max_match = max(max_match, match)

                dict_num_matchs[max_match] += (match_graph.nodes[dict_idx_node_pat[index]]["num_of_occurrence"] * max_x)


    print("status", lp_problem.status)
    print("num of covared:" , res_value)

    print("num that covarge: ", num_of_geno_that_covarge)

    for k, v in dict_num_matchs.items():
        print(f"{k},{v}")

    f_out.write(f"Cost/Number of needed genotypes: {res_value}\n")

    f_out.write(f"Percentage from the population: {pop_wanted_prec * 100}\n")
    f_out.write(f"Percentage with match of 3: {dict_num_matchs[3]/pop_size * 100}\n"
                f"Percentage with match of 4: {dict_num_matchs[4]/pop_size * 100}\n")

    pr.disable()
    pr.print_stats(sort="time")
