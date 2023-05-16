# NCOG Co-occurence network creation and visualization
The NOAA CalCOFI Ocean Genomics (NCOG) co-occurence network pipeline builds a co-occurence network of marine organisms in the California coastal region using [FlashWeave](https://github.com/meringlab/FlashWeave.jl), and generates circos plots of interactions between organisms using the [circlize](https://jokergoo.github.io/circlize_book/book/) package in R. 

## Setup
To start, you must have a folder with 16S, 18S, and metadata. The 16S and 18S data should have columns as samples, with rows as ASV's. In the data we used, the 16S dataset has a *silva_Taxon* and a *silva_Confidence* column for SILVA taxonomy prediction. The 18S datasets have a *pr2_Taxon* and a *pr2_Confidence* column for PR2 taxonomy prediction. These taxonomy columns are important as we use them to group ASVs by taxonomy, which is used in the circos plots to show the number of relations between species/taxonomic groups.

## Usage

### Step 1: FlashWeave ETL
To prepare the data for FlashWeave, we must have a 16S, 18S, and metadata that all represent the same samples and the same ASV's. We also must have the 16S and 18S tables transposed, so that columns represent ASV's and rows represent samples. Finally, we will filter out ASV's following the method published in [Chaffron, et al.](https://www.science.org/doi/10.1126/sciadv.abg1921): we calculate the Q3 of the nonzero values for each sample, and each ASV is retained if its abundance level is above the Q3 in at least five samples.

To run the FlashWeave ETL, run:

```bash
python3 flashweave_etl.py [sixteenS_path] [eighteenS_path] [metadata_path] [sixteenS_outpath] [eighteenS_outpath] [metadata_outpath] [id_map_outpath]
```

sixteenS_path: The full or relative path to the 16S file
eighteenS_path: The full or relative path to the 18S file
metadata_path: The full or relative path to the metadata file
sixteenS_outpath: The full or relative output path (including filename) to save the 16S filtered file
eighteenS_outpath: The full or relative output path (including filename) to save the 18S filtered file
metadata_outpath: The full or relative output path (including filename) to save the metadata filtered file
id_map_outpath: The full or relative output path (including filename) to save the id map file (this maps Feature.ID values to their potu/eotu value for later reference)

Example:

```bash
python3 flashweave_etl.py final_data/NCOG_21_16S_redo2_asv_count_tax.tsv final_data/NCOG_18sV4_asv_count_tax.tsv final_data/NCOG_sample_log_DNA_stvx_meta_2014-2020_mod.tsv final_filtered_data/NCOG_21_16S_v4_redo2_filtered_asv_count_tax.tsv final_filtered_data/NCOG_18sV4_filtered_asv_count_tax.tsv final_filtered_data/NCOG_v4_sample_log_DNA_stvx_meta_filtered_2014-2020_mod.tsv final_filtered_data/NCOG_v4_id_map.tsv
```

### Step 2: Run FlashWeave
To run FlashWeave, you likely need a supercomputer center that has julia installed or loaded. Follow instructions on the [FlashWeave GitHub](https://github.com/meringlab/FlashWeave.jl) to install the FlashWeave package in julia. Then run:

```julia
using FlashWeave
DIR = {Input Directory}
OUTDIR = {Output Directory}
data_path = [string(DIR, {16S File}), string(DIR, {18S File})]
meta_data_path = string(DIR, {Metadata File})
netw_results = learn_network(data_path, meta_data_path, sensitive=true, heterogeneous=false)
G = graph(netw_results)
save_network(string(OUTDIR, "NCOG_network_output.gml"), netw_results)
```

Example:
```julia
using FlashWeave
DIR = "./projdata/masorens/final_filtered_data/"
OUTDIR = "./projdata/masorens/graph_output/"
data_path = [string(DIR, "NCOG_21_16S_v4_redo2_filtered_asv_count_tax.tsv"), string(DIR, "NCOG_18sV4_filtered_asv_count_tax.tsv")]
meta_data_path = string(DIR, "NCOG_v4_sample_log_DNA_stvx_meta_filtered_2014-2020_mod.tsv")
netw_results = learn_network(data_path, meta_data_path, sensitive=true, heterogeneous=false)
G = graph(netw_results)
save_network(string(OUTDIR, "NCOG_network_output_v4.gml"), netw_results)
```

Once this has finished, copy the .gml file to a local directory where you will do the analysis

### Step 3: Create Circos Plots
To create the circlize circos plots, create a jupyter notebook like *Circlize ETL.ipynb* notebook found in *jupyter_notebooks*. Here are the circlize_etl functions:

1. `circlize_etl.get_edgelist(16s_filtered, 18s_filtered, metadata_filtered)`: this returns an edgelist with the columns `node1, node2, node1_Taxonomy, node2_Taxonomy, weight`
2. `circlize_etl.main(16s_filtered, 18s_filtered, network_file.gml, focus_group, top_n)`: This builds the edgelist for the circos plot of the focus group and the top_n groups connected to it

Once you have run `circlize_etl.main` then run the R code in `circlize_plot.R` inside RStudio. Because it can only display one plot at a time, you can uncomment each chordDiagram() line one at a time, zoom into the plot and right click it and save it as an image, then run the next chordDiagram() line.

### Step 4: Network Analysis
The final step is to do some network analysis. In the network_analysis.R script you will find calculations for the following network properties:

1. [Average clustering coefficient](https://en.wikipedia.org/wiki/Clustering_coefficient) (aka transitivity): measure of the degree to which nodes tend to cluster together - a higher value means nodes have highly inter-connected neighborhoods, lower value means less connections within neighborhoods
2. [Average path length](https://en.wikipedia.org/wiki/Average_path_length): average number of steps along the shortest paths for all possible pairs of network nodes
3. [Diameter](https://mathworld.wolfram.com/GraphDiameter.html): the longest path of the shortest paths between any two network nodes
https://igraph.org/r/doc/distances.html
4. [Average betweenness](https://en.wikipedia.org/wiki/Betweenness_centrality): average number of shortest paths that pass through each node
5. [Modularity of positive network](https://en.wikipedia.org/wiki/Modularity_(networks)): measures the strength of division of a network into modules. A higher number represents a positive network (edges with weigth > 0) that is highly connected to itself and more separable from the overall network
6. Modularity of negative network: A higher number represents a negative network (edges with weight < 0) that is highly connected to itself and more separable from the overall network