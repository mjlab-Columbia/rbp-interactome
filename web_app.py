import utils
from interactome import interactome as ict
import pandas as pd
from os.path import dirname, join
import math
import numpy as np

import networkx as nx
import bokeh
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource,CheckboxButtonGroup, Div,\
    Slider,ColumnDataSource,LabelSet,TextInput,Circle,Button,Select,\
    PreText,CheckboxGroup,LegendItem,Legend,CustomJS,HoverTool
from bokeh.layouts import column, row
from bokeh.plotting import figure

from operator import itemgetter
from mycolorpy import colorlist as mcp

##############
###DATASETS###
##############

NW = pd.read_csv('data/input/full_network.txt',delimiter='\t')
NW['IP_interaction_type'] = NW['IP_interaction_type'].fillna('SEC')
LC = pd.read_csv('data/input/life_cycle_annotations.csv')

LOC = dict(zip(NW['Prey'],NW['Prey_main_location']))
LOC.update(dict(zip(NW['Bait'],NW['Bait_main_location'])))
LOC = {str(k): str(LOC[k]) for k in LOC if ((LOC[k] != np.nan) and (k != np.nan))}

DIS = dict(zip(NW['Prey'],NW['Prey_Disgenet_disease_.1']))
DIS.update(dict(zip(NW['Bait'],NW['Bait_Disgenet_disease_.1'])))
DIS = {str(k): str(DIS[k]) for k in DIS if ((DIS[k] != np.nan) and (k != np.nan))}

baits = pd.read_csv('data/input/baits.txt',delimiter='\t',header=None)[0]
CORUM = pd.read_csv('data/references/corum_interactors.csv')

FEATURE_MAP = {'lifecycle stage': dict(zip(LC['protein'],LC['life_cycle_stage'])),'location': LOC,'disease': DIS} 

# Create Column Data Source that will be used by the plot
source_edges = ColumnDataSource(dict(xs=[], ys=[],style=[],color=[],label=[]))
source_nodes = ColumnDataSource(dict(xs=[], ys=[],color=[],name=[],node_size=[],label=[],Proteins=[]))
source_graph = ColumnDataSource(pd.DataFrame([]))
source_legend = ColumnDataSource(dict(xs=[],ys=[],color=[]))
source_line_legend = ColumnDataSource(dict(xs=[],ys=[],color=[],style=[],label=[]))

# Hover Capabilities
EDGE_HOVER_TOOLTIPS=[]
NODE_HOVER_TOOLTIPS=[]
node_hover_tool = HoverTool(tooltips=NODE_HOVER_TOOLTIPS)

#############
###BUTTONS###
#############

# Dataset
LABELS_DS= ['IP', 'SEC','CORUM']
dataset = CheckboxButtonGroup(labels=LABELS_DS,active=[0,1,2])

# Interaction Type
LABELS_IT= ["direct","mediated","shielded"]
interactionType = CheckboxButtonGroup(labels=LABELS_IT,active=[0,1,2])

# layout
layout = Select(title="Graph Layout", options=['Spring', 'Circle'], value="Spring")

# Subgraph search
proteins = TextInput(title="Protein")

# Clustering options
clustering = Select(title="Clustering method", options=['no clustering','louvain','markov','clusterONE'], value="no clustering")

# Clustering resolution
resolution = Slider(title="Clustering resolution", value=0.0, start=0.0, end=1.0, step=0.1)

# Overlay options
colorby = Select(title="Color node by", options=['none','lifecycle stage','location','disease'], value="none")

# Number neighbors
neighbors = Select(title="Number of neighbors", options=['1','2','3'], value="1")

# checkbox 
checkbox = CheckboxGroup(labels=['LABELS'], active=[0])

# Run button
button = Button(label="Create Interactome", button_type="success")

# Download button
button_d = Button(label="Download", button_type="success")

# Statistics
stats = PreText(text='Graph Statistics:', width=300)

################
###Formatting###
################

desc = Div(text=open(join(dirname(__file__), "Interactome.html")).read(), sizing_mode="stretch_width")

plt = figure(height=800,width=1200,tools=["pan","wheel_zoom","save","reset"], active_scroll='wheel_zoom',tooltips=[("Proteins", "@Proteins")])
plt.xgrid.grid_line_color = None
plt.ygrid.grid_line_color = None
plt.axis.visible = False

# Renders edges and nodes
plt.multi_line('xs', 'ys', line_color='color', line_dash='style',source=source_edges, alpha=0.5,width=2)
plt.circle('xs', 'ys', fill_color='color', line_color='black', size='node_size',source = source_nodes)

dum_fig = figure(width=None,height=None,toolbar_location=None)
crcls = dum_fig.circle('xs', 'ys', fill_color='color', line_color='black',size=0, legend_field='label',source = source_legend)

# Renders labels
plt.renderers.append(LabelSet(x='xs', y='ys', text='name', source=source_nodes, background_fill_color='white', text_font_size='8px', background_fill_alpha=.7))

# Dummy Figures
dum_fig.legend[0].location = "bottom"
dum_fig.legend[0].orientation = "horizontal"
dum_fig.legend[0].padding = 0

dum_fig2 = figure(width=None,height=None,toolbar_location=None)
lns = dum_fig2.multi_line('xs','ys',line_color='color',line_dash='style',legend_field='label',source=source_line_legend)
dum_fig2.legend[0].location = 'top_right'


plt.add_layout(dum_fig.legend[0],'below')
plt.add_layout(dum_fig2.legend[0])
plt.renderers.extend([crcls,lns])

################
## Functions ###
################

def retrieve_attr(attr,func,graph,node_or_edge='Node'):
    '''
    Retrieves attribute and sets default attribute value if attribute isn't present.
    :param: attribute to check.
    :param: default value.
    '''
    if node_or_edge == 'Node':
        attr = list(nx.get_node_attributes(graph,attr).values())

        if len(attr) == 0:
            attr = [func(node) for node in graph.nodes()]
    else:
        attr = list(nx.get_edge_attributes(graph,attr).values())

        if len(attr) == 0:
            attr = [func(edge) for edge in graph.edges()]
    
    return attr

def get_edge_information(graph, layout):
    '''
    Returns NetworkX graph edge information for graphing.
    :param: graph: The networkX graph object to retrieve information from.
    :param: layout: The networkX layout type to retrieve coordinate information.
    '''

    attributes = {}

    # Axis information
    attributes['xs'] = [[layout[edge[0]][0], layout[edge[1]][0]] for edge in graph.edges()]
    attributes['ys'] = [[layout[edge[0]][1], layout[edge[1]][1]] for edge in graph.edges()]

    # Line Style and Color
    attributes['style'] = retrieve_attr('style',lambda x: 'solid',graph,node_or_edge='Edge')
    attributes['color'] = retrieve_attr('color',lambda x: 'black',graph,node_or_edge='Edge')

    return attributes

def get_node_information(graph,layout,feature):
    '''
    Returns NetworkX graph Node information for graphing.
    :param: graph: The networkX graph object to retrieve information from.
    :param: layout: The networkX layout type to retrieve coordinate information.
    '''
    attributes = {}

    name, coords = zip(*layout.items())
    attributes['xs'],attributes['ys'] = zip(*coords)

    # Attributes you want to retrieve
    attributes['names'] = [list(graph.nodes())[i] for i in range(len(name))]

    attributes['Proteins'] = retrieve_attr('Proteins',lambda x: x, graph)
    attributes['size'] = [ math.log(len(x.strip('][').split(', ')),10) * 12 if '[' in str(x) else 12 for x in attributes['Proteins'] ]
    attributes['size'] =  [i if i >= 12 else 12 for i in attributes['size']]

    # Handles coloring of Nodes based on Attribute
    node_legend = retrieve_attr('color',lambda x: 'unspecified',graph)
    attributes['label'] = [x if x != 'lightsteelblue' else 'unspecified' for x in node_legend]
    
    # unique functions
    if feature != 'none':
        node_attr = set(FEATURE_MAP[feature].values())
        unique_c = dict(zip(node_attr,mcp.gen_color(cmap="tab20c",n=len(node_attr))))
        colors = [x if x == 'lightsteelblue' else unique_c[x] for x in node_legend]
    else:
        colors = [x for x in node_legend]
    
    if len(colors) == 0:
        colors = ['lightsteelblue' for node in graph.nodes()]

    attributes['color'] = colors

    return attributes

def select_criteria():
    '''
    Selects Criteria from dataFrame to display on interactive tool.
    :return: NetworkX Graph object, Layout, Bool to show if graph is clustered.
    '''

    CLUSTERED = False

    # Retrieve values from buttons
    interactionType_val = [LABELS_IT[x] for x in interactionType.active]
    dataset_val = [LABELS_DS[x] for x in dataset.active]
    layout_val = layout.value
    neighbors_val = int(neighbors.value)
    proteins_val = proteins.value.strip()
    resolution_val = resolution.value
    clustering_val = clustering.value
    feature = colorby.value

    display_s = {'direct':'solid','mediated':'dashed','shielded':'dotted','undetermined':'dashdot','SEC':'solid'}
    display_c = {'IP':'blue','SEC':'grey','both':'orange',True:'green'}

    ProteinA = NW[ NW['IP_interaction_type'].isin(interactionType_val + ['SEC'])]['Bait']
    ProteinB = NW[ NW['IP_interaction_type'].isin(interactionType_val + ['SEC'])]['Prey']

    CORUM = NW[ NW['IP_interaction_type'].isin(interactionType_val + ['SEC'])]['in_CORUM_2022']
    SUPPORT = NW[ NW['IP_interaction_type'].isin(interactionType_val + ['SEC'])]['Interaction_support']
    interaction = NW[ NW['IP_interaction_type'].isin(interactionType_val + ['SEC'])]['IP_interaction_type']

    styles = dict( zip( zip(ProteinA,ProteinB) , interaction) )
    colors = dict( zip( zip(ProteinA,ProteinB) , SUPPORT) )

    styles = {k:display_s[v] for k,v in styles.items()}
    colors = {k:display_c[v] for k,v in colors.items()}

    graph = ict(utils.to_adjacency_list(ProteinA,ProteinB))

    graph.setAttribute({node:'lightsteelblue' for node in graph.G.nodes},'color')
    graph.setAttribute(colors,'color',attr='edge')
    graph.setAttribute(styles,'style',attr='edge')

    if feature != 'none':
         graph.setAttribute(FEATURE_MAP[feature],'color')

    # Protein search
    if proteins_val != '':
        val = proteins_val.split(',')
        graph = ict(graph.subGraphNeighbors(val,hops=neighbors_val))

    # Clustering control
    if clustering_val != 'no clustering':
        graph.cluster(clustering_val,resolution=resolution_val)
        graph.induced()
        CLUSTERED = True

    remove = None

    if 'CORUM' in dataset_val:
        corum = {k:display_c[v] for k,v in dict( zip( zip(ProteinA,ProteinB) , CORUM) ).items() if v == True}
        colors = nx.get_edge_attributes(graph.G,'color')
        colors.update(corum)
        graph.setAttribute(colors,'color',attr='edge')

    if 'IP' in dataset_val and 'SEC' not in dataset_val:
        remove = [k for k,v in colors.items() if v == display_c['SEC']]

    if 'SEC' in dataset_val and 'IP' not in dataset_val:
        remove = [k for k,v in  colors.items() if v != display_c['SEC'] and v != display_c['both']]

    if remove:
        graph.G.remove_edges_from(remove)

    return graph,CLUSTERED,layout_val,feature

def graph_statistics(graph_full,edges,CLUSTERED):
    '''
    Retrieves specific Graph Statistics

    '''
    if CLUSTERED:
        graph = ict(graph_full.induced_g).G
    else:
        graph = graph_full.G

    # Populate Statistics string
    stat_string = 'Graph Statistics:\n\nNumber of Nodes: ' + str(len(graph.nodes())) + '\n' \
                  'Number of Interactions: ' + str(len(graph.edges())) + '\n\n' 
    labels = []
    counts = {'Direct: ' :0,'RNA mediated: ':0, 'RNA shielded: ':0,'Undetermined: ':0,'Spacer':None, \
        'Shared interactions: ':0, 'IP only interactions: ':0,'SEC only interactions: ':0,'CORUM interactions: ':0}

    for idx,edge in enumerate(edges['color']):

        int_type = ''
        if edge == 'green':
            int_type ='CORUM ' 
            counts['CORUM interactions: '] += 1

        elif edge == 'blue':
            int_type = 'IP '
            counts['IP only interactions: '] += 1

        elif edge == 'grey':
            int_type = 'SEC '
            counts['SEC only interactions: '] += 1

        elif edge == 'orange':
            int_type = 'Shared '
            counts['Shared interactions: '] += 1

        label = ''
        if edges['style'][idx] == 'solid' and edge != 'grey':
            label = 'Direct'
            counts['Direct: '] += 1

        elif edges['style'][idx] == 'dashed':
            label = 'mediated'
            counts['RNA mediated: '] += 1

        elif edges['style'][idx] == 'dotted':
            label = 'shielded'
            counts['RNA shielded: '] += 1
        elif edges['style'][idx] == 'dashdot':
            label = 'undetermined'
            counts['Undetermined: '] += 1

        labels.append(int_type + label)

    for type,c in counts.items():
        if c == None:
            stat_string += '\n'
        elif c > 0:
            stat_string += (type + str(c) + '\n')

    if CLUSTERED:
          # find betweeness nodes
        bw_dict = nx.betweenness_centrality(graph_full.G) 
        stat_string += '\nTop Five Betweeness Proteins:\n' + '\n'.join([str(x[0]) for x in sorted(bw_dict.items(), key=itemgetter(1), reverse=True)[:5]])

        stat_string += '\n\nNumber of Complexes: ' + str(len(graph.nodes())) + '\n'
        bw_dict = nx.betweenness_centrality(graph)
        stat_string += '\nTop Five Betweeness Complexes:\n' + '\n'.join([str(x[0]) for x in sorted(bw_dict.items(), key=itemgetter(1), reverse=True)[:5]])

    else:
        # find betweeness nodes
        bw_dict = nx.betweenness_centrality(graph) 
        stat_string += '\nTop Five Betweeness Proteins:\n' + '\n'.join([str(x[0]) for x in sorted(bw_dict.items(), key=itemgetter(1), reverse=True)[:5]])


    return stat_string,labels

def update():
    '''
    This function updates the Data Sources to display the new graph based on information retrieved from select criteria.
    '''
    # Retrieve filtered information
    graph,CLUSTERED,layout_val,feature = select_criteria()

    # Network X graph retrieval
    if CLUSTERED:
        graph_g = ict(graph.induced_g).G
    else:
        graph_g = graph.G

    # Layout control
    if layout_val == 'Spring':
        plot_layout = nx.spring_layout(graph_g) 
    elif layout_val == 'Circle':
        plot_layout = nx.circular_layout(graph_g)

    # Node and Edge information
    _nodes = get_node_information(graph_g,plot_layout,feature)
    _edges = get_edge_information(graph_g, plot_layout)

    # Retrieve Graph Statistics
    stats.text,labels = graph_statistics(graph,_edges,CLUSTERED)

    # Download Button
    fh = open('IP_and_SEC.txt','wb')
    nx.write_edgelist(graph_g,fh,data=True)

    # check if labels need to be shown
    if 0 not in checkbox.active:
        names = ['' for node in _nodes['names']] 
    else:
        names = _nodes['names']

    # Legend handling
    color_map = set(zip(_nodes['label'],_nodes['color']))
    source_legend.data = dict(xs=[0 for x in color_map], ys=[0 for x in color_map], color=[c[1] for c in color_map], \
        label = [c[0] for c in color_map])

    fig_colors = ['green','blue','grey','orange','black','black','black','black']
    fig_styles = ['solid','solid','solid','solid','solid','dotted','dashed','dashdot']
    fig_labels = ['corum','IP','SEC','Shared','Direct','RNA Mediated','RNA Shielded','Undetermined']

    source_line_legend.data = dict(xs=[0 for x in fig_colors], ys=[0 for x in fig_colors], color=fig_colors, \
        label = fig_labels,style=fig_styles)

    # Update Column Data Source
    source_edges.data = dict(xs=_edges['xs'], ys=_edges['ys'],style=_edges['style'],color=_edges['color'],label=labels)
    source_nodes.data = dict(xs=_nodes['xs'], ys=_nodes['ys'],name=names,color=_nodes['color'],node_size=_nodes['size'],label=_nodes['label'],Proteins=_nodes['Proteins'])

    # Download button
    source_graph.data = pd.DataFrame({'x':[x[0] for x in graph_g.edges()],'y':[x[1] for x in graph_g.edges()],'edge_type':labels})

# Populate Buttons
controls= [dataset,interactionType,proteins,neighbors,clustering,resolution,colorby,layout,checkbox,button,button_d]
button.on_click(update)

# Download Button
button_d.js_on_click(CustomJS(args=dict(source=source_graph), code=utils.JAVSCRIPT_CODE))

# Place components on website.
inputs = column(*controls, width=300)
stats_input = column(*[stats],width=300)
l = column(desc, row(inputs, plt,stats_input), sizing_mode="scale_both")

# Bokeh handling
curdoc().add_root(l)
curdoc().title = "Interactome"