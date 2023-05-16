import networkx as nx
import numpy as np
import pandas as pd
import sys

def main(sixteens_path, eighteens_path, graph_path, focus_group, top_n):
    sixteenS = pd.read_csv(sixteens_path, sep='\t')
    eighteenS = pd.read_csv(eighteens_path, sep='\t')
    sixteenS.index = pd.Series(sixteenS.index).apply(lambda x: 'potu_' + str(x))
    eighteenS.index = pd.Series(eighteenS.index).apply(lambda x: 'eotu_' + str(x))
    gph = nx.read_gml(graph_path)
    # Add taxonomy info to the graph
    taxonDict = {}
    for node, taxon in sixteenS['silva_Taxon'].iteritems():
        taxonDict[node] = {'Taxon': taxon}
    for node, taxon in eighteenS['pr2_Taxon'].iteritems():
        taxonDict[node] = {'Taxon': taxon}
    nx.set_node_attributes(gph, taxonDict)
    # Split taxonomy information into domain, kingdom, phylum, class, etc.
    new_attrs = {}
    tax_attrs = nx.get_node_attributes(gph, 'Taxon')
    for node in tax_attrs:
        new_attrs[node] = {}
        for i, t in enumerate(tax_attrs[node].split(';')):
            new_attrs[node]['t' + str(i)] = t.strip()
    nx.set_node_attributes(gph, new_attrs)
    # Get edgelist
    edgelist = get_edgelist(sixteens_path, eighteens_path, graph_path)
    # Get diatom edges
    top_edge_counts = get_relevant_edges(new_attrs, edgelist, gph, focus_group, top_n)
    if top_edge_counts == -1:
        return
    focus_edges_df = pd.DataFrame({'Edge': top_edge_counts.keys(), 'Counts': top_edge_counts.values()})
    focus_edges_df['Node1'] = focus_edges_df['Edge'].apply(lambda x: x[0])
    focus_edges_df['Node2'] = focus_edges_df['Edge'].apply(lambda x: x[1])
    focus_edges_df['Positive'] = focus_edges_df['Counts'].apply(lambda x: x['positive'])
    focus_edges_df['Negative'] = focus_edges_df['Counts'].apply(lambda x: x['negative'])
    focus_edges_df = focus_edges_df.drop(['Counts', 'Edge'], axis=1)
    focus_edges_df['Node1'] = focus_edges_df['Node1'].apply(lambda x: x.strip().strip('c__'))
    focus_edges_df['Node2'] = focus_edges_df['Node2'].apply(lambda x: x.strip().strip('c__'))
    focus_edges_df['Node1'] = focus_edges_df['Node1'].apply(lambda x: x.replace('(', '_').replace(')', '_'))
    focus_edges_df['Node2'] = focus_edges_df['Node2'].apply(lambda x: x.replace('(', '_').replace(')', '_'))
    focus_edges_df = focus_edges_df.sort_values(['Node1', 'Node2']).reset_index(drop=True)
    return focus_edges_df
    
def get_edgelist(sixteens_path, eighteens_path, graph_path):
    """
    Generate edgelist node1, node2, node1Taxonomy, node2Taxonomy, weight from sixteens, eighteens, graph files
    """
    # Load files
    sixteenS = pd.read_csv(sixteens_path, sep='\t')
    eighteenS = pd.read_csv(eighteens_path, sep='\t')
    sixteenS.index = pd.Series(sixteenS.index).apply(lambda x: 'potu_' + str(x))
    eighteenS.index = pd.Series(eighteenS.index).apply(lambda x: 'eotu_' + str(x))
    gph = nx.read_gml(graph_path)
    # Add taxonomy info to the graph
    taxonDict = {}
    for node, taxon in sixteenS['silva_Taxon'].iteritems():
        taxonDict[node] = {'Taxon': taxon}
    for node, taxon in eighteenS['pr2_Taxon'].iteritems():
        taxonDict[node] = {'Taxon': taxon}
    nx.set_node_attributes(gph, taxonDict)
    # build edgelist
    nodetaxons = nx.get_node_attributes(gph, 'Taxon')
    edgelist = nx.to_pandas_edgelist(gph, source='node1', target='node2')
    edgelist['node1Taxonomy'] = edgelist['node1'].apply(lambda x: nodetaxons[x] if x in nodetaxons else np.nan)
    edgelist['node2Taxonomy'] = edgelist['node2'].apply(lambda x: nodetaxons[x] if x in nodetaxons else np.nan)
    edgelist = edgelist[['node1', 'node2', 'node1Taxonomy', 'node2Taxonomy', 'weight']]
    return edgelist

def get_relevant_edges(new_attrs, edgelist, gph, focus_group, top_n):
    focus_index = -1
    for edge in nx.edges(gph):
        if edge[0] in new_attrs and focus_group in new_attrs[edge[0]].values():
            focus_index = list(new_attrs[edge[0]].values()).index(focus_group)
            break
        elif edge[1] in new_attrs and focus_group in new_attrs[edge[1]].values():
            focus_index = list(new_attrs[edge[1]].values()).index(focus_group)
            break
    if focus_index == -1:
        print('Focus group not found in graph')
        return -1
    edgelist_new = edgelist.copy()
    edgelist_new['node1TaxonomyInterest'] = edgelist_new['node1Taxonomy'].apply(get_group_of_interest,
                                                                                args=[focus_group, focus_index])
    edgelist_new['node2TaxonomyInterest'] = edgelist_new['node2Taxonomy'].apply(get_group_of_interest,
                                                                                args=[focus_group, focus_index])
    edgelist_new = edgelist_new.dropna(subset=['node1TaxonomyInterest','node2TaxonomyInterest'])
    # Get set of all nodes that are the focus group or are connected to the focus group, eventually will remove all 
    # entries from edgelist that are not the focus group or connected to the focus group
    # Also get count of how many edges each group has to the focus group
    nodes_connected_to_focus = set()
    group_counts = {}
    for _, row in edgelist_new.iterrows():
        # Add node 2 if node 1 is focus group
        if row['node1TaxonomyInterest'] == focus_group:
            nodes_connected_to_focus.add(row['node2'])
            nodes_connected_to_focus.add(row['node1'])
            if row['node2TaxonomyInterest'] not in group_counts:
                group_counts[row['node2TaxonomyInterest']] = 0
            group_counts[row['node2TaxonomyInterest']] += 1
        if row['node2TaxonomyInterest'] == focus_group:
            nodes_connected_to_focus.add(row['node1'])
            nodes_connected_to_focus.add(row['node2'])
            if row['node1TaxonomyInterest'] != focus_group:
                if row['node1TaxonomyInterest'] not in group_counts:
                    group_counts[row['node1TaxonomyInterest']] = 0
                group_counts[row['node1TaxonomyInterest']] += 1
    # Get top n groups connected to focus group
    group_counts_df = pd.DataFrame({'Taxonomic Group': group_counts.keys(), 'Edges to Focus': group_counts.values()})
    group_counts_df.sort_values(by='Edges to Focus', ascending=False, inplace=True)
    top_groups = set(group_counts_df.reset_index(drop=True)[:top_n]['Taxonomic Group'])
    top_groups.add(focus_group)
    # Iterate through edgelist and drop edges where either node is not in top_n
    droprows = []
    for i, row in edgelist_new.iterrows():
        if row['node1'] not in nodes_connected_to_focus or row['node2'] not in nodes_connected_to_focus:
            droprows.append(i)
        if row['node1TaxonomyInterest'] not in top_groups or row['node2TaxonomyInterest'] not in top_groups:
            droprows.append(i)
    edgelist_new = edgelist_new.drop(index=droprows)
    # Get total counts of every group-group edge for nodes in the top_n groups
    top_edge_counts = {}
    for i, row in edgelist_new.iterrows():
        node1 = row['node1TaxonomyInterest']
        node2 = row['node2TaxonomyInterest']
        weight = 'positive' if row['weight'] >= 0 else 'negative'
        # store each group, group pair in dict with 
        if node1 > node2:
            node1, node2 = node2, node1
        if (node1, node2) not in top_edge_counts:
            top_edge_counts[(node1, node2)] = {'positive': 0, 'negative': 0}
        top_edge_counts[(node1, node2)][weight] += 1
    return top_edge_counts

def get_group_of_interest(taxonomy, focus_group, focus_index):
    if taxonomy != taxonomy:
        return taxonomy
    tax_split = taxonomy.split(';')
    if len(tax_split) >= focus_index + 1:
        if tax_split[focus_index].strip() == focus_group:
            return focus_group
    # Return third level for Eukaryotes, except for Bacillariophyta and Syndiniales
    if tax_split[0] == 'Eukaryota':
        if len(tax_split) >= 4:
            if tax_split[3].strip() == 'Bacillariophyta':
                return 'Bacillariophyta'
            elif tax_split[3].strip() == 'Syndiniales':
                return 'Syndiniales'
        if len(tax_split) >= 3:
            return tax_split[2].strip()
        else:
            return tax_split[-1].strip()
    # Return third level for Bacteria/Archaea
    else:
        if len(tax_split) >= 3:
            return tax_split[2].strip()
        else:
            return tax_split[-1].strip()