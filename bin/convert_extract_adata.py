#!/usr/bin/env python3

import argparse
import sys
import gzip

import anndata as ad
import fast_matrix_market as fmm
import pandas as pd
import scanpy as sc
import scipy
import decoupler as dc


def parse_args(args=None):
    Description = "Extract and save raw counts as matrix, barcodes, features, metadata, and umap coordinates from an adata object seperately"
    Epilog = "Example usage: python convert_extract_adata.py <adata>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--adata")
    return parser.parse_args(args)

def pseudobulk(adata,groups_col,sample_col):
    # Pseudobulk by sample type 
  
    layer="counts"
    pseudobulk = [
        (
            group,
            dc.get_pseudobulk(
                adata,
                sample_col=sample_col,
                groups_col=groups_col,
                layer=layer,
                mode="sum",
                min_prop=0.05,
                min_cells=10,
                min_counts=500,
                min_smpls=10,
            ),
        )
        for group in adata.obs[groups_col].unique()
    ]
    pdata = pseudobulk[0][1]
    #pdata.layers['counts']=pdata.X.copy()
    return pdata



def convert_to_CSC(data_mat: scipy.sparse.spmatrix) -> scipy.sparse.csc_matrix:
    """
    Converts a sparse matrix to Compressed Sparse Column (CSC) format.
    See https://github.com/MarioniLab/scran/issues/70

    Parameters:
    data_mat (scipy.sparse.spmatrix): The input sparse matrix to be converted.

    Returns:
    scipy.sparse.csc_matrix: The input matrix converted to CSC format.

    Notes:
    - If the input matrix is already in CSC format, it is returned unchanged.
    - If the input matrix is in another sparse format (e.g., COO), it is converted to CSC format.
    - If the input matrix is not sparse, an error will be raised.
    """

    if scipy.sparse.issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()

    return data_mat


def extract_adata(adata):


    #fmm.mmwrite(
    #    "counts_matrix.mtx",
    #    convert_to_CSC(scipy.sparse.csr_matrix(adata.layers["counts"].T)),
    #)

    adata.copy().T.to_df().to_csv('counts.tsv')
    
    pd.DataFrame(adata.var.index).to_csv(
        "features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )
   
    pd.DataFrame(adata.obs.index).to_csv(
        "barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )
 
    adata.obs.to_csv("metadata.tsv", sep="\t", index=True)
    adata.var.to_csv("metadata_var.tsv", sep="\t", index=True)

    #pd.DataFrame(adata.obsm["X_umap"]).to_csv(
    #    "umap.tsv.gz", sep="\t", index=False, compression="gzip"
    #)


def main(args=None):
    args = parse_args(args)
    sample_col ="sample_id"
    groups_col="sample_type"

    adata=sc.read_h5ad(args.adata)
    adata =  adata[adata.obs["enrichment_cell_types"]=="naive"]
   
    pdata = pseudobulk(adata,groups_col,sample_col)
    extract_adata(pdata)


if __name__ == "__main__":
    sys.exit(main())
