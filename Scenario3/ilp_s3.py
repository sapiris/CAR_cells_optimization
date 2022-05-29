from partial_graph_match_s3 import Graph_Match
import pulp
import timeit
import json


def init_lp_prolem(match_graph, pop_size, dict_idx_node_don, dict_idx_node_pat, dict_node_idx, lim1, knockout_cost):
    lp_problem = pulp.LpProblem("LP_Problem", pulp.LpMinimize)

    ##variables boundaries
    var = []
    for i in range(len(dict_idx_node_don)):
        var.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Integer'))
    for i in range(len(dict_idx_node_pat)):
        var.append(pulp.LpVariable('y' + str(i), lowBound=0, upBound=1, cat='Integer'))
    ##equation to maximize
    z = 0
    for idx, node in dict_idx_node_don.items():
        cost = 1
        if len(node.split('~')) == 9:
            cost = knockout_cost

        z += (var[idx] * cost)
    lp_problem += z

    ##Constraints
    # (sum(x_i) <= k
    k = pop_size * lim1
    sum_x_i = -1 * var[len(dict_idx_node_don)] * match_graph.nodes[dict_idx_node_pat[0]]["num_of_occurrence"]
    for i in range(1 + len(dict_idx_node_don), len(dict_idx_node_don) + len(dict_idx_node_pat)):
        sum_x_i -= (var[i] * match_graph.nodes[dict_idx_node_pat[i - len(dict_idx_node_don)]]["num_of_occurrence"])
    lp_problem += sum_x_i <= (-1 * k)

    k = pop_size * lim1
    sum_x_i = var[len(dict_idx_node_don)] * match_graph.nodes[dict_idx_node_pat[0]]["num_of_occurrence"]
    for i in range(1 + len(dict_idx_node_don), len(dict_idx_node_don) + len(dict_idx_node_pat)):
        sum_x_i += (var[i] * match_graph.nodes[dict_idx_node_pat[i - len(dict_idx_node_don)]]["num_of_occurrence"])
    lp_problem += sum_x_i <= (k)

    for i in range(len(dict_idx_node_pat)):
        cons_vec = 0
        node = dict_idx_node_pat[i]
        # kir_node = match_graph.nodes[node]["kir"]
        # set_i = find_adjs(node, match_graph)
        set_i = match_graph.adj[node]
        for adj in set_i:
            if not "pat" in adj:
                cons_vec += (-1 * var[dict_node_idx[adj]])
        cons_vec += var[len(dict_idx_node_don) + i]

        lp_problem += cons_vec <= 0

    # A*02:01~a2~b2~c1~c2~q1~q2~r1~r1

    print(lp_problem)

    return lp_problem, var




if __name__ == '__main__':

    with open('configuration_s3.json') as f:
        conf = json.load(f)

    create_cost = 1
    knockout_cost = conf["knockout_cost"] + create_cost
    pop_wanted_prec = conf["percentage_from_the_population"] /100 # the precentage of the population to cobvare
    don_in_file = conf["donors_input_file"]#'../../data/israel/don_israel_5_geno_kir_' + sim + '.csv'
    rec_in_file = conf["recipients_input_file"]


    f_out = open(conf["output_file"], "w")


    g_m = Graph_Match(0, 0, directed= False)
    match_graph, dict_idx_node_don, dict_idx_node_pat, dict_node_idx, pop_size = g_m.create_graph(don_in_file, rec_in_file)
    lp_problem, var = init_lp_prolem(match_graph, pop_size, dict_idx_node_don, dict_idx_node_pat,
                                                            dict_node_idx, pop_wanted_prec, knockout_cost)

    lp_problem.solve()

    ll = lp_problem.variables()

    num_of_geno_that_covarge = 0
    res_value = pulp.value(lp_problem.objective)
    sum_9 = 0
    sum_full = 0
    dict_y_res = {}
    f_out.write("selected genotypes:\n")
    for i, variable in enumerate(lp_problem.variables()):
        print("{} = {}".format(variable.name, variable.varValue))
        # if variable.varValue == 1:
        if "x" in variable.name:
            index = int(variable.name.replace('x', ''))
            if variable.varValue == 1:
                if len(dict_idx_node_don[index].split('~')) == 9:
                    sum_9 += variable.varValue
                else:
                    sum_full += variable.varValue

                f_out.write(dict_idx_node_don[index] + '\n')

    print("sum_9", sum_9)
    print("sum_full", sum_full)
    print("status", pulp.LpStatus[lp_problem.status])
    print("num of covared:", res_value)

    f_out.write(f"Cost: {res_value}\n")
    f_out.write(f"Percentage from the population: {pop_wanted_prec*100}\n")
    f_out.write(f"Number of choose knockout: {sum_9}\n")


