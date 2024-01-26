#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:09:46 2024

@author: gouthamvasam
"""

import scanpy as sc

def scRNAseq_workflow(counts_matrix_dir, output_dir):
    # Load the count matrix into an AnnData object
    adata = sc.read_10x_mtx(
        counts_matrix_dir,  # the directory with the `.mtx` file
        var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
        cache=True)    # write a cache file for faster subsequent reading

    # Quality control
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Calculate the percentage of mitochondrial genes expressed
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)

    # Filter out cells with high mitochondrial gene expression
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Logarithmize the data
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]

    # Dimensionality reduction
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.tsne(adata, n_pcs=40)

    # Clustering
    sc.tl.leiden(adata, resolution=0.5)

    # Differential expression analysis
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    # Annotate cell clusters using known marker genes (this requires a list of marker genes)
    # marker_genes_dict = {'CellType1': ['GeneA', 'GeneB'], 'CellType2': ['GeneC', 'GeneD']}
    # sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden')

    # Save the results
    adata.write(output_dir + '/scRNAseq_analysis_results.h5ad')

    # Plotting the results
    sc.pl.umap(adata, color=['leiden'], save='_clusters.png')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_markers.png')

def main():
    # Define the directory containing the count matrix and the directory to save the output
    counts_matrix_dir = 'path/to/counts_matrix'
    output_dir = 'path/to/output'

    # Run the scRNA-seq workflow
    scRNAseq_workflow(counts_matrix_dir, output_dir)

if __name__ == "__main__":
    main()
