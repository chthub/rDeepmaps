# Controller for iris3api
options(digits = 3)
require(iris3api, quietly = TRUE)

#* @apiTitle Plumber Example API

#* Log requests
#* filter log_request
#* log_request

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
echo_back

#* Plot a histogram
#* @serializer png
#* @get /plot
plot_histogram

#* Return the sum of two numbers
#* @param a The first number to add
#* @param b The second number to add
#* @post /sum
sum_numbers

####################################

#* Read data into  object
#* @param filename
#* @param type Upload expression file type, CellGene, 10X h5, 10X folder
#* @post /load
load_single_rna

#* Read multiple scRNA-seq
#* @param filename
#* @param type Upload expression file type, CellGene, 10X h5, 10X folder
#* @post /load-multi-rna
load_multi_rna

#* Read data into  object
#* @param filename
#* @param type Upload expression file type, CellGene, 10X h5, 10X folder
#* @post /load-multiome
load_multiome

#################################### Clustering part

#* Run RNA sample clustering
#* @param nPCs
#* @param resolution
#* @param neighbor
#* @post /cluster
cluster_single_rna

#* Run multiome clustering
#* @param nPCs
#* @param resolution
#* @param neighbor
#* @post /cluster-multiome
cluster_multiome

#* Plot umap
#* @get /umap-cluster
active_label

#* Merge idents
#* @param newClusterIds
#* @post /merge-idents
merge_idents

#* Rename idents
#* @param old_name
#* @param new_name
#* @post /rename-idents
rename_idents

#* Select-category
#* @param categoryName
#* @post /select-category
select_category

#* Select-cells
#* @param categoryName
#* @param filterPayload
#* @post /select-cells
select_cells

#* Subset cells
#* @param selectionPayload
#* @post /subset-cells
subset_cells

####################################

#* Get Variable genes list
#* @get /var-genes-list
rna_qc_list

#* Plot Variable genes statc plot
#* @get /var-genes-plot
#* @serializer png list(width = 800, height = 600)
rna_qc_plot

#* Get ATAC QC list
#* @get /atac-qc-list
atac_qc_list


#################################### Regulon

#* Calculate regulon
#* @param cluster string
#* @serializer unboxedJSON
#* @post /regulon
example_regulon_network

#* Return all regulons
#* @serializer unboxedJSON
#* @post /all-regulon
list_regulon_network

#* Calculate regulon clusters umap
#* @post /regulon-clusters
example_cluster_coords

#* Calculate RAS
#* @param gene
#* @param assay
#* @post /ras
example_ras

#* Calculate GAS
#* @param gene
#* @param assay
#* @post /gas
example_gas

#* Calculate DR
#* @param tf
#* @param ct1
#* @param ct2
#* @post /dr
example_dr

#* Calculate DR
#* @param genes
#* @post /ri-heatmap
example_ri_heatmap
#################################### DEG

#* Calculate DEG for two selections
#* @post /deg
calc_deg

#################################### Enrichment

#* Calculate GSEA table
#* @param genes
#* @param database
#* @post /gsea-table
calc_gsea_table

#* Calculate GSEA table
#* @param genes
#* @param database
#* @post /regulon-enrichment
calc_enrichr_table

#* Enrichr dot plot
#* @param df
#* @param isPvalLog
#* @post /enrichment-dotplot
#* @serializer png list(width = 4000, height = 2000, res = 300)
plot_enrichr_dot

#* Enrichr bar plot
#* @param df
#* @param isPvalLog
#* @post /enrichment-barplot
#* @serializer png list(width = 4000, height = 2000, res = 300)
plot_enrichr_bar

#* GSEA plot
#* @param df
#* @param isPvalLog
#* @post /enrichment-gsea-plot
#* @serializer png list(width = 3000, height = 2000, res = 300)
plot_gsea

#################################### Object manipulation

#* Set object to full or subset
#* @param type
#* @post /set-obj
set_obj

#* Get all gene names
#* @get /genes
get_all_genes

#* Get allIdents names
#* @get /idents
get_all_idents

#* Set idents
#* @param name
#* @post /idents
set_idents

#* Get all assay names
#* @get /assays
get_all_assays

#* Set assay names
#* @param name
#* @post /set-assay
set_assay

#* Get all embedding names
#* @get /embeddings
get_all_embeddings

#* Get all  Idents names
#* @param name
#* @post /set-embedding
set_embedding

#################################### Visualization

#* Plot gene-correlation plot
#* @post /gene-correlation-plot
#* @serializer png list(width = 1600, height = 1200, res = 300)
gene_cor_plot

#* Plot umap
#* @param categoryName
#* @post /umap-static
#* @serializer png list(width = 700, height = 600)
umap_plot

#* Plot umap gene plot
#* @param gene
#* @post /gene-umap-static
#* @serializer png list(width = 700, height = 600)
gene_umap_plot

#* Plot violin gene plot
#* @param gene
#* @param split
#* @param group
#* @post /violin-gene
#* @serializer png list(width = 700, height = 600)
violin_gene_plot

#* Heatmap
#* @param genes
#* @post /static-heatmap
#* @serializer png list(width = 1000, height = 600)
static_heatmap

#* Heatmap
#* @param genes
#* @post /interactive-heatmap
#* @serializer html
interactive_heatmap

#* Coverage plot
#* @param type
#* @param gene
#* @param flank
#* @param chr
#* @param start
#* @param end
#* @param is_annotation
#* @param is_peak
#* @post /coverage-plot
#* @serializer png list(width = 1200, height = 1000)
coverage_plot

#* Return embeddings coordinates from a given embedding name
#* @param embedding
#* @post /embedding-coords
#* @serializer json
cluster_coords

#* Return embeddings coordinates from a given embedding name
#* @param gene
#* @param assay
#* @post /feature-coords
#* @serializer json
feature_coords


#################################### CellChat
#* Run
#* @param job
#* @param active_idents
#* @post /run-cellchat
run_cellchat

#* Plot
#* @param mode
#* @param slot
#* @param pathway_show
#* @param lr_show
#* @param ident1
#* @param ident2
#* @post /plot-cellchat
#* @serializer png list(width = 700, height = 600)
plot_cellchat
