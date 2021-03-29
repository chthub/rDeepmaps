# Controller for iris3api

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

#* Read data into Seurat object
#* @param filename
#* @param type Upload expression file type, CellGene, 10X h5, 10X folder
#* @post /load
load_seurat

#* Read data into Seurat object
#* @param nPCs
#* @param resolution
#* @param neighbor
#* @post /cluster
cluster_seurat

#* Plot umap
#* @get /umap-cluster
active_label

#* Merge idents
#* @param newClusterIds
#* @post /merge-idents
merge_idents

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

#* Plot umap
#* @param categoryName
#* @post /umap-static
#* @serializer png list(width = 700, height = 600)
umap_plot

#* Get Variable genes list
#* @get /var-genes-list
rna_qc_list

#* Get all Seurat Idents names
#* @get /idents
get_all_idents

#* Get all gene names
#* @get /genes
get_all_genes

#* Get all Seurat Idents names
#* @post /idents
set_idents

#* Get all Seurat Idents names
#* @param type
#* @post /set-obj
set_obj

#* Calculate DEG for two selections
#* @post /deg
calc_deg
