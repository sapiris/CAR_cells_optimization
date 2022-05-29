import time
import pickle
import networkx as nx

class Graph_Match(object):

    def __init__(self,  missmatch_c1,missmatch_c2, alleles_in_donor, alleles_in_patient, directed = True):
        self.miss_c1 = missmatch_c1
        self.miss_c2 = missmatch_c2
        self.graph = nx.DiGraph() if directed else nx.Graph()
        self.list_assistance_nodes = []
        self.allele_in_don = alleles_in_donor
        self.allele_in_pat = alleles_in_patient

    def create_sub_geno_node(self,full_node,list_alleles,class_idx, miss_c, don_or_pat):
        class_node = ('~').join(list_alleles)
        if not class_node in self.graph.nodes():
            self.graph.add_node(class_node, alleles=class_idx)

        num_of_allels = len(list_alleles)
        if don_or_pat == 'don':
            self.graph.add_edge(full_node, class_node, weight=class_idx)
        else:
            self.graph.add_edge(class_node, full_node, weight=class_idx)
        # create temporal nodes of all combination of class_node
        if num_of_allels > ((4-miss_c) if class_idx == 'class1' else (4-miss_c)):
            for i in range(len(list_alleles)):
                sub_class_node = list_alleles[0:i] + list_alleles[i + 1:num_of_allels]
                self.create_sub_geno_node(full_node, sub_class_node, class_idx, miss_c, don_or_pat)



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
                if self.allele_in_don in geno:
                    geno = geno.replace('+', '~').replace('^', '~')
                    geno = ('~').join(sorted(geno.split('~')))
                    if geno in self.graph:
                        self.graph.nodes[geno]["num_of_occurrence"] +=1
                    else:
                        dict_idx_node_don[idx] = geno
                        dict_node_idx[geno] = idx
                        idx+=1
                        self.graph.add_node(geno, MUUG= geno,  alleles='full_don', num_of_occurrence=1 )
                        geno_class1 = sorted(geno.split('~'))[0:4]
                        #geno_class2 = geno.split('~')[6:10]
                        self.create_sub_geno_node(geno, geno_class1, 'class1', self.miss_c1, 'don')
                        #self.create_sub_geno_node(geno, geno_class2, 'class2', self.miss_c2, 'don')

        idx = 0
        with open(f_data_pat) as f_in:
            for i,line in enumerate(f_in):
                id, geno = line.strip().split(',')[0:2]
                if self.allele_in_pat in geno:
                    geno = geno.replace('+', '~').replace('^', '~')
                    geno = ('~').join(sorted(geno.split('~')))
                    geno_name = 'pat-' + geno
                    if geno_name in self.graph:
                        self.graph.nodes[geno_name]["num_of_occurrence"] +=1
                    else:
                        dict_idx_node_pat[idx] = geno_name
                        idx += 1
                        self.graph.add_node(geno_name, MUUG= geno_name,alleles='full_pat', num_of_occurrence=1 )
                        geno_class1 = sorted(geno.split('~'))[0:4]
                        #geno_class2 = geno.split('~')[6:10]
                        self.create_sub_geno_node(geno_name, geno_class1, 'class1', self.miss_c1, 'pat')
                        #self.create_sub_geno_node(geno, geno_class2, 'class2', self.miss_c2, 'pat')




        print("graph creation time: ", time.time()-t)

        return self.graph, dict_idx_node_don, dict_idx_node_pat, dict_node_idx, i+1


