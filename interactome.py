"""Handling of Interactome"""

__author__ = 'Ahmed Abdou'

import networkx as nx
import markov_clustering as mc
import statistics
import os
import utils

class interactome:
    
    def __init__(self,data,weight=None):
        '''
        Initalizes an interactome object from provided data.
        :param data: Data to be fed to networkx.
        '''
        if type(data) is nx.Graph:
            self.G = data

        self.G = nx.Graph(data)

        # Check to see if graph is weighted
        self.weight = weight
        self.graphStatistics()

        self.clustered = False

    def write_graph(self):
        '''
        writes the graph to a data.txt file to be used for processing.

        :return: True if successful.
        '''
        fh = open('data.txt','wb')
        nx.write_edgelist(self.G,fh,data=False)
        fh.close()

        return True

    def cluster(self,cluster_type='markov',resolution=1.1):
        '''
        :param cluster_type: Clustering method, the options are Louvain, Markov, and ClusterOne.
        '''

        print(f'Clustering using the {cluster_type} method')

        if cluster_type == 'louvain':
            self.clusters = nx.algorithms.community.louvain_communities(self.G,weight=self.weight,resolution= 1 + 5.0 * resolution)
            print(f'Detected {len(self.clusters)} complexes')

        elif cluster_type == 'markov':
            matrix = nx.to_scipy_sparse_matrix(self.G) 
            result = mc.run_mcl(matrix,inflation=(1.5 + resolution*5.0))
            nodes = list(self.G.nodes)
            self.clusters = [[ nodes[g] for g in group] for group in mc.get_clusters(result)]
            print(f'Detected {len(self.clusters)} complexes')
            
        elif cluster_type == 'clusterONE':
            result = self.write_graph()
            if result:
                command = 'java -jar clusterONE/cluster_one-1.0.jar data.txt > clusters.txt'
                os.system(command)

                self.clusters = []
                fh = open('clusters.txt','r')
                for line in fh:
                    self.clusters.append(line.split())
                fh.close()

                os.system('rm clusters.txt data.txt')

        # Update Graph Statistics
        cluster_lengths = [len(x) for x in self.clusters]

        self.num_clusters = len(self.clusters)
        self.cluster_min = min(cluster_lengths)
        self.cluster_max = max(cluster_lengths)
        self.cluster_median = statistics.median(cluster_lengths)
        self.clustered = True

    def graphStatistics(self):
        '''
        This function finds the graph statistics.
        '''
        self.num_edges = nx.number_of_edges(self.G)
        self.num_nodes = nx.number_of_nodes(self.G)
        self.num_clusters = None
        self.cluster_min = None
        self.cluster_max = None
        self.cluster_median = None

    def subGraphNeighbors(self,nbrs,hops=1):
        '''
        Returns the subgraph of all nodes connected to the specified list.
        :param nbrs: the list of nodes to return the subgraph.

        :return: A networkx subgraph of nodes in and connected to nbrs.
        '''
        nodes = set([x for x in nbrs if len(x) != 1])
        updated_nodes = set([x for x in nbrs if len(x) != 1])

        for i in range(hops):
            nodes = updated_nodes.copy()

            for nbr in nodes:
                try:
                    for edge in self.G.edges(nbr):
                        color= self.G.get_edge_data(*edge)['color']

                        if color != 'grey':
                            updated_nodes.add(edge[1])
                except:
                    continue

        return self.G.subgraph(list(updated_nodes))

    def setAttribute(self,labels,name,attr='node'):
        '''
        Assigns an attribute to the node
        '''
        if attr == 'node':
            nx.set_node_attributes(self.G, labels, name)
        elif attr == 'edge':
            nx.set_edge_attributes(self.G, labels, name)
        else:
            raise TypeError('Attribute type must be \'node\' or \'edge\'')


    def __add__(self,G2):
        '''
        This overides the add operator to add two graphs together.
        '''
        H = nx.compose(self.G,G2.G)

        s1 = set(self.G.nodes)
        s2 = set(G2.G.nodes)

        color1 = nx.get_node_attributes(self.G,'color')
        color2 = nx.get_node_attributes(G2.G,'color')
        intersect_nodes = s1.intersection(s2)

        color_labels = dict(color1, **color2)

        for node in intersect_nodes:
            color_labels[node] = 'green'

        H = interactome(H)
        H.setAttribute(color_labels,'color')

        return H

    def subtract(self,G2):
        '''
        '''
        e1 = set([tuple(sorted(x)) for x in self.G.edges])
        e2 = set([tuple(sorted(x)) for x in G2.edges])

        return list(e1.difference(e2))

    def intersection(self,G2):
        '''
        Find the intersection of two graphs.
        :param G2: networkx graph to analyze in conjuction with.
        :return: nodes that are interesecting, edges that are intersecting, 
        networkxgraph from interesecting edges
        '''

        s1 = set(self.G.nodes)
        s2 = set(G2.nodes)

        e1 = set([tuple(sorted(x)) for x in self.G.edges])
        e2 = set([tuple(sorted(x)) for x in G2.edges])

        intersect_nodes = s1.intersection(s2)
        intersect_edges = e1.intersection(e2)

        return intersect_nodes, intersect_edges, nx.from_edgelist(intersect_edges)

    def induced(self):
        '''
        Finds an induced graph from the clustered graph for plotting purposes.
        :return: true or false based on success.
        '''
        if self.clustered == False:
            print('You must cluster the interactome first!')
            return False
        else:
            self.partition = {}

            for i,cluster in enumerate(self.clusters):
                for protein in cluster:
                    self.partition[protein] = i

            self.induced_g = nx.Graph()
            self.induced_g.add_nodes_from(self.partition.values())

            for node1, node2, datas in self.G.edges(data=True):
                edge_weight = datas.get('weight', 1)
                edge_color = datas.get('color','black')
                try:
                    com1 = self.partition[node1]
                    com2 = self.partition[node2]
                    w_prec = self.induced_g.get_edge_data(com1, com2, {'weight': 0}).get('weight', 1)
                    connecting = self.induced_g.get_edge_data(com1,com2,{'Connecting_Protein':''}).get('Connecting_Protein')
                    self.induced_g.add_edge(com1, com2, **{'weight': w_prec + edge_weight,'Connecting_Protein': connecting + str(node1) + ' to ' + str(node2) + ' , ','color':edge_color})
                except:
                    pass

            # calculate comm color
            node_colors = nx.get_node_attributes(self.G,'color')
            node_color_m = {}
            for idx,cluster in enumerate(self.clusters):
                colors = [node_colors[protein] for protein in cluster]
                node_color = max(colors,key=colors.count)
                node_color_m[idx] = node_color
            
            self.clusters = [[str(x) for x in y] for y in self.clusters]
            nx.set_node_attributes(self.induced_g,node_color_m,'color')
            nx.set_node_attributes(self.induced_g,{k:str(sorted(v)) for k,v in enumerate(self.clusters)},'Proteins')

            utils.name_complexes(self)
            map_unknown = {}
            for node in self.induced_g.nodes():
                prots = self.induced_g.nodes()[node]['Proteins'].strip('][').split(', ')

                if type(node) is int and len(prots) > 1:
                    map_unknown[node] = 'Unmapped Compelex ' + str(node)
                elif type(node) is int and len(prots) == 1:
                    map_unknown[node] = prots[0].replace("'",'')

            nx.relabel_nodes(self.induced_g,map_unknown,copy=False)

            return True

    def resolution_combination(self,G2,baits,corum_support=False):
        '''
        Combines IP and SEC data in a way to enrich resolution
        :param:
        :param:
        :param:
        '''
        H = nx.Graph()

        for bait in baits:
            bait_subgraph = self.subGraphNeighbors([bait])
            ip_nodes = set(bait_subgraph.nodes())
            ip_edges = set([tuple(sorted(x)) for x in bait_subgraph.edges()])

            # sec subgraph on bait IP.
            sec_subgraph = G2.G.subgraph(ip_nodes)
            sec_edges = set([tuple(sorted(x)) for x in sec_subgraph.edges()])

            ip_unique_edges = ip_edges - sec_edges
            sec_unique_edges = sec_edges - ip_edges
            intersecting_edges = ip_edges.intersection(sec_edges)
            all_edges = ip_edges.union(sec_edges)

            edge_label = {}
            for edge in all_edges:
                if edge in ip_unique_edges:
                    edge_label[edge] = 'solid'#'dotted'
                elif edge in sec_unique_edges:
                    edge_label[edge] = 'dashed'
                if edge in intersecting_edges:
                    edge_label[edge] = 'dotted'

            # Add edges and style
            H.add_edges_from(list(all_edges))
            nx.set_edge_attributes(H,edge_label,'line_style')

        nx.set_node_attributes(H,{n:'lightsteelblue' for n in H.nodes()},'color')
        return H
        

    def complexScore(self,reference='data/references/corum.txt'):
        '''
        Calculates the clustering efficacy using clusterONE metrics.
        :return: true or false based on success.
        '''
        if self.clustered == False:
            print('You must cluster the interactome first!')
            return False

        self.write_graph()

        # write cluster results
        fh = open('clusters.txt','w')
        for cluster in self.clusters:
            fh.write(' '.join(list(cluster)) + '\n')
        fh.close()

        # execute mwmatching statistics
        command = "python clusterONE/match_standalone.py -q -n data.txt -m frac -m cws -m ppv -m acc -m mmr " + reference + " clusters.txt"
        os.system(command)
        os.system('rm clusters.txt data.txt')

        return True

    def return_information():
        '''
        '''
        info = {}
    
