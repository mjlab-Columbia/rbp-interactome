"""Utils for Project"""

__author__ = 'Ahmed Abdou'
__author__ = 'Wenhao Jin'

from typing import DefaultDict
import pandas as pd
from collections import Counter
import networkx as nx
from bokeh.models import CustomJS

def to_adjacency_list(origin,destination):
    '''
    Returns an adjacency list from lists of origin and destination.
    '''
    return list(zip(origin,destination))

def to_adj_dict(df,weighted=True,feature='Bait_main_location'):
    '''
    Takes dataframe of Bait-Prey pairs and returns an adjcaency dictionary.
    '''
    # Find unique baits
    adj_dict = {}
    features = {}
    # For each bait, find the list of prey it catches.
    for bait in df['Bait'].unique().tolist(): 
        targets = {}
        # For each prey find it's weight
        for prey in df[df['Bait'] == bait]['Prey'].tolist():
            if bait != prey:

                log2FC = df[(df['Bait'] == bait) & (df['Prey'] == prey)]['log2FC'].tolist()[0]
                colors = df[(df['Bait'] == bait) & (df['Prey'] == prey)]['Interaction_type'].tolist()[0]

                if feature=='location':
                    try:
                        features[prey] =  df[(df['Bait'] == bait) & (df['Prey'] == prey)]['Prey_main_location'].tolist()[0].split(';')[0]
                    except:
                        pass
                elif feature =='lifecycle stage':
                    try:
                        features[prey] =  df[(df['Bait'] == bait) & (df['Prey'] == prey)]['Life_cycle_step'].tolist()[0].split(';')[0]
                    except:
                        pass
                elif feature == 'disease':
                    try:
                        features[prey] =  df[(df['Bait'] == bait) & (df['Prey'] == prey)]['Prey_Disgenet_disease_#1'].tolist()[0].split(';')[0]
                    except:
                        pass

                # Set colors for connection Type
                if colors == 'Direct':
                    colors = 'red'
                elif colors == 'RNA mediated':
                    colors = 'green'
                else:
                    colors = 'blue'
                
                # Add log2FC scores
                if weighted:
                    targets[prey] = {'weight': log2FC,'color':colors}
                else:
                    targets = df[(df['Bait'] == bait) &(df['Prey'] == prey)]['Bait'].tolist()
                
        adj_dict[bait] = targets

        if feature =='location':
            try:
                features[bait] =  df[(df['Bait'] == bait)]['Bait_main_location'].tolist()[0].split(';')[0]
            except:
                pass
        
        elif feature == 'lifecycle stage':
            try:
                features[bait] =  df[(df['Bait'] == bait)]['Life_cycle_step'].tolist()[0].split(';')[0]
            except:
                pass

        elif feature == 'disease':
            try:
                features[bait] =  df[(df['Bait'] == bait)]['Bait_Disgenet_disease_#1'].tolist()[0].split(';')[0]
            except:
                pass
        
    return adj_dict,features

def name_complexes(graph):
    '''
    '''
    CORUM = pd.read_csv('data/references/humanComplexes.txt',delimiter='\t')

    assigns= DefaultDict(list)
    for genes in CORUM['subunits(Gene name)']:
        lst = genes.split(';')
        lst = [graph.partition.get(x,None) for x in lst]
        prots = Counter(lst)

        if prots[None]/len(lst) <= .9:
            if prots.most_common(1)[0][0] == None:
                comm = int(prots.most_common(2)[1][0])
                
                try:
                    size = assigns[comm][1]
                    if int(prots.most_common(2)[1][1]) > size:
                        assigns[comm] = (CORUM[CORUM['subunits(Gene name)'] == genes]['ComplexName'].tolist()[0], int(prots.most_common(2)[1][1]))
                except:
                    assigns[comm] = (CORUM[CORUM['subunits(Gene name)'] == genes]['ComplexName'].tolist()[0], int(prots.most_common(2)[1][1]))
            
            else:
                
                #pdb.set_trace()
                comm = int(prots.most_common(1)[0][0])
                
                try:
                    size = assigns[comm][1]
                    if int(prots.most_common(1)[0][1]) > size:
                        assigns[comm] = (CORUM[CORUM['subunits(Gene name)'] == genes]['ComplexName'].tolist()[0], int(prots.most_common(1)[0][1]))
                except:
                    assigns[comm] = (CORUM[CORUM['subunits(Gene name)'] == genes]['ComplexName'].tolist()[0], int(prots.most_common(1)[0][1]))

                
        
        assigns = {k:str(v[0]) for k,v in assigns.items()}

        nx.relabel_nodes(graph.induced_g,assigns,copy=False)

JAVSCRIPT_CODE="""
var data = source.data;
var filetext = data
var filetext = 'protein1,protein2,interactionType\\n'

for (let i=0; i < data['x'].length; i++) {

    var currRow = [data['x'][i].toString(), data['y'][i].toString(),data['edge_type'][i].toString().concat('\\n')];
	var joined = currRow.join();
	filetext = filetext.concat(joined);
}	

var filename = 'interactome_data.csv';
var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename);
}

else {
    var link = document.createElement("a");
    link = document.createElement('a')
    link.href = URL.createObjectURL(blob);
    link.download = filename
    link.target = "_blank";
    link.style.visibility = 'hidden';
    link.dispatchEvent(new MouseEvent('click'))
}
"""