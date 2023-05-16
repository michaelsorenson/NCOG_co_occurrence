import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import sys

def main(sixteens_path, eighteens_path, graph_path, focus_group, top_n):
    """
    Note that top_n might get 1 less than the number specified if the focus_group itself is in the top_n groups
    """
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
    # Plot distribution of number of edges
    fig = plt.figure(figsize=(15, 4))
    ax = fig.add_subplot(111)
    y = nx.degree_histogram(gph)
    sns.barplot(x=list(range(len(y))), y=y, ax=ax)
    fig.savefig('Graph_edgecounts.png')
    # Get edge weights
    edgeweights = nx.get_edge_attributes(gph, 'weight')
    # Get number of positive and negative edges
    pos_edges = 0
    neg_edges = 0
    for edge, weight in edgeweights.items():
        if weight >= 0:
            pos_edges += 1
        else:
            neg_edges += 1
    # Build graph information table
    table1 = {
        'Nodes': [nx.number_of_nodes(gph)],
        'Edges': [nx.number_of_edges(gph)],
        'Positive Edges (%)': [str(pos_edges) + ' ' + str(pos_edges/(pos_edges + neg_edges))],
        'Negative Edges': [neg_edges],
        'Average clustering coefficient': [nx.average_clustering(gph)],
        'Average path length': [nx.average_shortest_path_length(gph)],
        'Diameter': [nx.diameter(gph)],
        'Average betweenness': [np.nan],
        'Modularity of positive network': [np.nan],
        'Number of modules in positive network': [np.nan]
    }
    pd.DataFrame(table1).to_csv('Graph_information.csv')
    # Split taxonomy information into domain, kingdom, phylum, class, etc.
    new_attrs = {}
    tax_attrs = nx.get_node_attributes(gph, 'Taxon')
    for node in tax_attrs:
        new_attrs[node] = {}
        for i, t in enumerate(tax_attrs[node].split(';')):
            new_attrs[node]['t' + str(i)] = t.strip()
    nx.set_node_attributes(gph, new_attrs)
    # Build counts of each taxonomic group
    tax_counts = {}
    t1_tax_counts = {}
    for edge in nx.edges(gph):
        # Get weight
        weight_num = edgeweights[edge]
        weight = ''
        if weight_num >= 0:
            weight = 'positive'
        else:
            weight = 'negative'
        # Get the highest level taxon edge pair, ignore edges involving metadata
        if edge[0] in new_attrs and edge[1] in new_attrs:
            node1_tax1 = new_attrs[edge[0]]['t0']
            node2_tax1 = new_attrs[edge[1]]['t0']
            if node1_tax1 <= node2_tax1:
                if (node1_tax1, node2_tax1) not in t1_tax_counts:
                    t1_tax_counts[(node1_tax1, node2_tax1)] = {'positive': 0, 'negative': 0}
                t1_tax_counts[(node1_tax1, node2_tax1)][weight] += 1
            else:
                if (node2_tax1, node1_tax1) not in t1_tax_counts:
                    t1_tax_counts[(node2_tax1, node1_tax1)] = {'positive': 0, 'negative': 0}
                t1_tax_counts[(node2_tax1, node1_tax1)][weight] += 1
            # Add 1 to each taxonomic subgroup
            for t in new_attrs[edge[0]]:
                t_val = new_attrs[edge[0]][t]
                if t_val != '':
                    if t_val not in tax_counts:
                        tax_counts[t_val] = {'positive': 0, 'negative': 0}
                    tax_counts[t_val][weight] += 1
            for t in new_attrs[edge[1]]:
                t_val = new_attrs[edge[1]][t]
                if t_val != '':
                    if t_val not in tax_counts:
                        tax_counts[t_val] = {'positive': 0, 'negative': 0}
                    tax_counts[t_val][weight] += 1
    tax_counts_df = pd.DataFrame({'Taxonomic Group': tax_counts.keys(), 'Counts': tax_counts.values()})
    tax_counts_df['Positive'] = tax_counts_df['Counts'].apply(lambda x: x['positive'])
    tax_counts_df['Negative'] = tax_counts_df['Counts'].apply(lambda x: x['negative'])
    tax_counts_df = tax_counts_df.drop('Counts', axis=1)
    tax_counts_df = tax_counts_df.sort_values('Positive', ascending=False).reset_index(drop=True)
    tax_counts_df.to_csv(f'Graph_taxonomic_group_counts.csv')
    # Get edgelist
    edgelist = get_edgelist(sixteens_path, eighteens_path, graph_path)
    edgelist.to_csv('NCOG_network_edgelist.csv', index=False)
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
    # Write circos data
    write_circos_data(top_edge_counts, focus_edges_df, focus_group)

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
    edgelist_new = edgelist.copy()
    edgelist_new['node1TaxonomyInterest'] = edgelist_new['node1Taxonomy'].apply(get_group_of_interest,
                                                                                args=[focus_group, focus_index])
    edgelist_new['node2TaxonomyInterest'] = edgelist_new['node2Taxonomy'].apply(get_group_of_interest,
                                                                                args=[focus_group, focus_index])
    # Get set of all nodes connected to focus group, remove all entries from edgelist not including those nodes
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
    # Drop rows that have nodes not connected to the focus group or don't have both nodes in the top_n groups
    droprows = []
    for i, row in edgelist_new.iterrows():
        if row['node1'] not in nodes_connected_to_focus or row['node2'] not in nodes_connected_to_focus:
            droprows.append(i)
    edgelist_new = edgelist_new.drop(index=droprows).dropna(subset=['node1TaxonomyInterest','node2TaxonomyInterest'])
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
        elif row['node1TaxonomyInterest'] not in top_groups or row['node2TaxonomyInterest'] not in top_groups:
            droprows.append(i)
    edgelist_new = edgelist_new.drop(index=droprows)
    # Get total counts of every group-group edge for nodes in the top_n groups
    top_edge_counts = {}
    for i, row in edgelist_new.iterrows():
        node1 = row['node1TaxonomyInterest']
        node2 = row['node2TaxonomyInterest']
        weight = 'positive' if row['weight'] >= 0 else 'negative'
        if node1 > node2:
            node1, node2 = node2, node1
        if (node1, node2) not in top_edge_counts:
            top_edge_counts[(node1, node2)] = {'positive': 0, 'negative': 0}
        top_edge_counts[(node1, node2)][weight] += 1
    return top_edge_counts
    
def get_group_counts(focus_df):
    for i, row in focus_df.iterrows():
        group1 = row['Edge'][0]
        group2 = row['Edge'][1]
        if not group1 in tax_group_counts:
            tax_group_counts[group1] = {'positive': 0, 'negative': 0}
        if not group2 in tax_group_counts:
            tax_group_counts[group2] = {'positive': 0, 'negative': 0}
        tax_group_counts[group1]['positive'] += row['Positive']
        tax_group_counts[group1]['negative'] += row['Negative']
        tax_group_counts[group2]['positive'] += row['Positive']
        tax_group_counts[group2]['negative'] += row['Negative']
    return tax_group_counts

def write_circos_data(top_edge_counts, focus_edges_df, focus_group):
    # Get counts for each group
    tax_group_counts = {}
    for i, row in focus_edges_df.iterrows():
        group1 = row['Node1']
        group2 = row['Node2']
        if not group1 in tax_group_counts:
            tax_group_counts[group1] = {'positive': 0, 'negative': 0}
        if not group2 in tax_group_counts:
            tax_group_counts[group2] = {'positive': 0, 'negative': 0}
        tax_group_counts[group1]['positive'] += row['Positive']
        tax_group_counts[group1]['negative'] += row['Negative']
        tax_group_counts[group2]['positive'] += row['Positive']
        tax_group_counts[group2]['negative'] += row['Negative']
    tax_group_counts_df = pd.DataFrame({'Group': tax_group_counts.keys(),
                                        'Positive': [x['positive'] for x in tax_group_counts.values()],
                                        'Negative': [x['negative'] for x in tax_group_counts.values()]
                                       })
    tax_group_counts_df['Total'] = tax_group_counts_df['Positive'] + tax_group_counts_df['Negative']
    # Build list of which color goes with which group, so when a link between two groups is defined, we know to make
    # it the color of the group with the higher total
    colors = [f'chr{x}_a2' for x in range(4, 22)]
    color_weights = tax_group_counts_df.sort_values('Total', ascending=False)
    color_col = [colors[i % len(colors)] for i in range(len(color_weights))]
    color_weights['Color'] = color_col
    color_weights.set_index('Group', drop=True, inplace=True)
    # Get dataframe of self-connected edges
    self_edges_df = (focus_edges_df[focus_edges_df['Node1'] == focus_edges_df['Node2']]
                     .set_index('Node1', drop=True)
                     .drop('Node2', axis=1)
                    )
    print(focus_edges_df['Node1'].unique())
    # Build circos circle
    with open('circos/simple_data.txt', 'w') as f:
        for i, group in enumerate(color_weights.index):
            # Write group length
            # Add 3 pixels between the self-edges for positive and negative
            if not group in self_edges_df:
                offset = 0
            else:
                offset = (int(self_edges_df.loc[group, 'Positive'] > 0) + int(color_weights.loc[group, 'Negative'] > 0)
                          * 3)
            group_count = tax_group_counts_df[tax_group_counts_df['Group'] == group]['Total'].iloc[0] + offset
            group_color = color_weights.loc[group, 'Color']
            f.write('chr - ' + group + ' ' + group + ' 0 ' + str(group_count) + 
                    ' ' + group_color + '\n')
    with open('circos/positive_links.txt', 'w') as pos, open('circos/negative_links.txt', 'w') as neg:
        links = {
            group: {'pos': 0, 'neg': tax_group_counts_df[tax_group_counts_df['Group'] == group]['Total'].iloc[0]}
            for group in tax_group_counts_df['Group']
        }
        for _, row in focus_edges_df.iterrows():
            if row['Positive'] > 0:
                node1_color = color_weights.loc[row['Node1'], 'Color']
                node2_color = color_weights.loc[row['Node2'], 'Color']
                node1_weight = color_weights.loc[row['Node1'], 'Total']
                node2_weight = color_weights.loc[row['Node2'], 'Total']
                # Assign color to color of higher weight group
                color = node1_color if node1_weight >= node2_weight else node2_color
                node1_start = links[row['Node1']]['pos']
                node2_start = links[row['Node2']]['pos']
                if row['Node1'] == row['Node2']:
                    node2_start = node1_start + row['Positive'] + 3
                pos.write(f'{row["Node1"]} {node1_start} {node1_start + row["Positive"]}' +
                          f' {row["Node2"]} {node2_start} {node2_start + row["Positive"]}' +
                          f' color={color}\n')
                links[row['Node1']]['pos'] += row['Positive']
                links[row['Node2']]['pos'] += row['Positive']
                if row['Node1'] == row['Node2']:
                    links[row['Node1']]['pos'] += 3
            if row['Negative'] > 0:
                node1_start = links[row['Node1']]['neg']
                node2_start = links[row['Node2']]['neg']
                if row['Node1'] == row['Node2']:
                    node2_start = node1_start - row['Negative'] - 3
                neg.write(f'{row["Node1"]} {node1_start - row["Negative"]} {node1_start}' +
                          f' {row["Node2"]} {node2_start - row["Negative"]} {node2_start}' +
                          ' color=vvdred\n')
                links[row['Node1']]['neg'] -= row['Negative']
                links[row['Node2']]['neg'] -= row['Negative']
                if row['Node1'] == row['Node2']:
                    links[row['Node1']]['neg'] -= 3
    
    
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

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Usage: python3 analysis.py sixteens_path eighteens_path graph_path focus_group')
        sys.exit(0)
    sixteens_path, eighteens_path, graph_path, focus_group = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    main(sixteens_path, eighteens_path, graph_path, focus_group)