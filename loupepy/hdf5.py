"""
reimplement
https://github.com/10XGenomics/loupeR/blob/v1.1.0/R/hdf5.R
"""

import re
import os
import sys
import h5py
import numpy as np
import pandas as pd
import anndata
from scipy.sparse import csc_matrix


def create_hdf5(adata, outfile="out.h5"):
    with h5py.File(outfile, "w") as out:
        matrix_to_h5(out, adata)
        clusters_to_h5(out, adata)
        projections_to_h5(out, adata)
        metadata_to_h5(out, adata)


def matrix_to_h5(h5_obj, adata):
    features = adata.var_names
    barcodes_unmodified = adata.obs_names.to_numpy()
    barcodes_formatted = format_barcodes(adata)

    sparse_mat = csc_matrix(adata.X)

    matrix_group = h5_obj.create_group("matrix")
    features_group = matrix_group.create_group("features")

    matrix_group.create_dataset("barcodes", data=barcodes_formatted)
    matrix_group.create_dataset("barcodes_unmodified", data=barcodes_unmodified)
    matrix_group.create_dataset("data", data=sparse_mat.data.astype(int))
    matrix_group.create_dataset("indices", data=sparse_mat.indices.astype(int))
    matrix_group.create_dataset("indptr", data=sparse_mat.indptr.astype(int))
    matrix_group.create_dataset("shape", data=np.array(adata.shape).astype(int))

    features_group.create_dataset("name", data=features.values)
    features_group.create_dataset("id", data=adata.var.gene_ids.values)
    features_group.create_dataset("feature_type", data=adata.var.feature_type.values)
    features_group.create_dataset("_all_tag_keys", data=np.empty())


def validate_barcodes(barcodes):
    bc_pattern = re.compile("^([ACTG]{6,})(-.*?)?$")
    return all([bc_pattern.match(bc) for bc in barcodes])


def format_barcodes(adata):
    barcodes = adata.obs_names
    if validate_barcodes(barcodes):
        return barcodes.to_numpy()

    bc_pattern = re.compile("^(.*?)(_|-|:)?([ACTG]{6,})(-[0-9]+)?(_|-|:)?(.*?)$")
    expanded = barcodes.str.extractall(bc_pattern).fillna("")
    expanded.columns = ["prefix", "sep1", "barcode", "dashnum", "sep2", "suffix"]
    if expanded.prefix.any():
        expanded.prefix = "-" + expanded.prefix
    if expanded.suffix.any():
        expanded.suffix = "-" + expanded.suffix

    formated = expanded.barcode + expanded.dashnum + expanded.prefix + expanded.suffix
    if len(formated) != len(barcodes):
        return barcodes.to_numpy()

    return formated.to_numpy()


def clusters_to_h5(h5_obj, adata):
    # clusters (aka all categorical metadata
    clusters_group = h5_obj.create_group("clusters")
    for colname in adata.obs_keys():
        dtype = adata.obs[colname].dtype
        if not isinstance(dtype, np.dtypes.ObjectDType()):
            continue
        if not isinstance(dtype, pd.CategoricalDtype()):
            continue

        as_cat = adata.obs[colname].astype("category")

        group = clusters_group.create_group(colname)
        group.create_dataset(group, "name", data=colname)
        group.create_dataset(group, "group_names", data=as_cat.cat.categories)
        group.create_datasets(group, "assignments", data=as_cat.cat.codes)
        group.create_datasets(group, "score", data=0.0)
        group.create_datasets(group, "clustering_type", data="unknown")


def projections_to_h5(h5_obj, adata):
    projections_group = h5_obj.create_group("projections")
    for proj_name in adata.obsm_keys():
        proj = adata.obsm[proj_name]
        is_umap = "umap" in proj_name
        is_tsne = ("tsne" in proj_name) or ("t-sne" in proj_name)

        method = proj_name
        if is_umap:
            method = "UMAP"
        if is_tsne:
            method = "t-SNE"

        group = projections_group.create_group("proj_name")
        group.create_dataset_like("name", data=proj_name)
        group.create_dataset_like("method", data=method)
        group.create_dataset_like("data", data=proj)


def metadata_to_h5(h5_obj, adata):
    info = sys.version_info()
    version = "{}.{}.{}".format(info.main, info.minor, info.micro)

    metadata = {
        "tool": "loupepy",
        "tool_version": "0.0.0",
        "os": os.name,
        "system": sys.platform,
        "language": "Python",
        "language_version": version,
        "extra": {
            "loupepy_anndata_vesion": anndata.__version__,
            "loupepy_hdf5_version": h5py.version.hdf5_version,
            "loupepy_h5py_version": h5py.version.verse,
        },
    }

    dict_to_h5(h5_obj, metadata, "metadata")


def dict_to_h5(parent, data, groupname):
    group = parent.create_group(groupname)

    for key, val in data.items():
        if isinstance(val, dict):
            dict_to_h5(group, val, key)
        else:
            group.create_dataset(key, data=val)
