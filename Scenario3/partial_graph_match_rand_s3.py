import time
import networkx as nx

class Graph_Match(object):

    def __init__(self,  missmatch_c1,missmatch_c2, directed = True):
        self.miss_c1 = missmatch_c1
        self.miss_c2 = missmatch_c2
        self.graph = nx.DiGraph() if directed else nx.Graph()
        self.list_assistance_nodes = []


    def create_sub_geno_node(self,full_node,list_alleles,class_idx, miss_c, don_or_pat, max_class_1 = 6, max_class_2=4):
        class_node = ('~').join(list_alleles)
        if not class_node in self.graph.nodes():
            self.graph.add_node(class_node, alleles=class_idx)
            self.list_assistance_nodes.append(class_node) #save in list of assistance to remove at the end

        num_of_allels = len(list_alleles)
        if don_or_pat == 'don':
            self.graph.add_edge(full_node, class_node, weight=class_idx)
        else:
            self.graph.add_edge(class_node, full_node, weight=class_idx)
        # create temporal nodes of all combination of class_node
        if num_of_allels > ((max_class_1-miss_c) if class_idx == 'class1' else (max_class_2-miss_c)):
            for i in range(len(list_alleles)):
                sub_class_node = list_alleles[0:i] + list_alleles[i + 1:num_of_allels]
                self.create_sub_geno_node(full_node, sub_class_node, class_idx, miss_c, don_or_pat,max_class_1,max_class_2)


    def find_adjs(self , node):
        dict_adjs = {}
        for assitance_adj in self.graph.adj[node]:
            edge = self.graph.get_edge_data(node, assitance_adj)
            # edage info- for know from which class
            for full_adj in self.graph.adj[assitance_adj]:
                if node != full_adj:

                    if full_adj in dict_adjs:
                        if not edge["weight"] in dict_adjs[full_adj]:
                            dict_adjs[full_adj] = ('_').join(sorted([dict_adjs[full_adj], edge["weight"]]))
                    else:
                        dict_adjs[full_adj] = edge["weight"]

        for adj, classes in dict_adjs.items():
            if classes == 'class1_class2':
                self.graph.add_edge(node, adj, weight=classes)

    #create graph: nodes of genotypes, edges- between 2 genotypes with 9 identical alleles
    def create_graph(self, f_data_don, f_data_pat):
        dict_idx_node_don = {}
        dict_idx_node_pat = {}
        dict_node_idx = {}
        t=time.time()


        idx = 0
        with open(f_data_don) as f_in:
            for i,line in enumerate(f_in):
                id, geno = line.strip().split(',')[0:2]
                geno = geno.replace('+', '~').replace('^', '~')
                geno = ('~').join(sorted(geno.split('~')))
                if geno in self.graph:
                    self.graph.nodes[geno]["num_of_occurrence"] +=1
                else:
                    dict_idx_node_don[idx] = geno
                    dict_node_idx[geno] = idx
                    idx+=1
                    self.graph.add_node(geno, MUUG= geno, alleles='full_don', num_of_occurrence=1 )
                    geno_class1 = geno.split('~')[0:6]
                    geno_class2 = geno.split('~')[6:10]
                    self.create_sub_geno_node(geno, geno_class1, 'class1', self.miss_c1, 'don')
                    self.create_sub_geno_node(geno, geno_class2, 'class2', self.miss_c2, 'don')
                    list_geno = geno.split('~')
                    for i in range(len(list_geno)):
                        geno_9 = list_geno[0:i] + list_geno[i + 1:10]
                        geno_9_joind = ('~').join(list_geno[0:i] + list_geno[i + 1:10])
                        if geno_9_joind in self.graph:
                            self.graph.nodes[geno_9_joind]["num_of_occurrence"] += 1
                        else:
                            dict_idx_node_don[idx] = geno_9_joind
                            dict_node_idx[geno_9_joind] = idx
                            idx += 1
                            self.graph.add_node(geno_9_joind, MUUG=geno_9_joind, alleles='knockout_don', num_of_occurrence=1)
                            geno_class1 = geno_9[0:5 + int(i/5)]
                            geno_class2 = geno_9[5 + int(i/5):9]
                            self.create_sub_geno_node(geno_9_joind, geno_class1, 'class1', self.miss_c1, 'don',len(geno_class1),len(geno_class2))
                            self.create_sub_geno_node(geno_9_joind, geno_class2, 'class2', self.miss_c2, 'don',len(geno_class1),len(geno_class2))


        idx  = 0
        with open(f_data_pat) as f_in:
            for i,line in enumerate(f_in):
                id, geno = line.strip().split(',')[0:2]
                geno = geno.replace('+', '~').replace('^', '~')
                geno = ('~').join(sorted(geno.split('~')))
                geno_name = 'pat-' + geno
                if geno_name in self.graph:
                    self.graph.nodes[geno_name]["num_of_occurrence"] +=1
                else:
                    dict_idx_node_pat[idx] = geno_name
                    idx+=1
                    self.graph.add_node(geno_name, MUUG= geno_name, alleles='full_pat', num_of_occurrence=1 )
                    geno_class1 =  sorted(geno.split('~'))[0:6]
                    geno_class2 =  sorted(geno.split('~'))[6:10]
                    self.create_sub_geno_node(geno_name, geno_class1, 'class1', self.miss_c1 + 1, 'pat')
                    self.create_sub_geno_node(geno_name, geno_class2, 'class2', self.miss_c2 + 1, 'pat')


        print(len(self.graph.nodes()))


        for node in dict_node_idx:
            self.find_adjs(node)


        #remove temporal nodes
        for node in self.list_assistance_nodes:
            self.graph.remove_node(node)


        print("graph creation time: ", time.time()-t)

        return self.graph, dict_idx_node_don, dict_idx_node_pat, dict_node_idx, i+1


