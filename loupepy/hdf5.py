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
from anndata import AnnData
from scipy.sparse import csc_matrix
from pandas.api.types import is_numeric_dtype
from pathlib import Path


def create_hdf5(adata: AnnData, outfile: Path | str = Path("out.h5")) -> None:
    with h5py.File(outfile, "w") as out:
        matrix_to_h5(out, adata)
        clusters_to_h5(out, adata)
        projections_to_h5(out, adata)
        metadata_to_h5(out)


def matrix_to_h5(h5_obj: h5py.File, adata: AnnData) -> None:
    n_bcs, n_vars = adata.shape
    barcodes_unmodified = adata.obs_names.to_numpy()
    barcodes_formatted = format_barcodes(adata)

    sparse_mat = csc_matrix(adata.X.T)

    matrix_group = h5_obj.create_group("matrix")
    features_group = matrix_group.create_group("features")

    create_string_dataset(matrix_group, "barcodes", barcodes_formatted)
    create_string_dataset(matrix_group, "barcodes_unmodified", barcodes_unmodified)
    matrix_group.create_dataset("data", data=sparse_mat.data.astype(np.int32))
    matrix_group.create_dataset("indices", data=sparse_mat.indices.astype(np.int32))
    matrix_group.create_dataset("indptr", data=sparse_mat.indptr.astype(np.int32))
    matrix_group.create_dataset(
        "shape", data=np.array([n_vars, n_bcs]).astype(np.int32)
    )

    if "gene_ids" not in adata.var_keys():
        gene_ids = np.array([f"feature_{k:05d}" for k in range(adata.shape[1])])
    else:
        gene_ids = adata.var.gene_ids.values

    if "feature_type" not in adata.var_keys():
        feat_types = np.repeat(["Gene Expression"], adata.shape[1])
    else:
        feat_types = adata.var.feature_type.values

    create_string_dataset(features_group, "name", adata.var_names.values)
    create_string_dataset(features_group, "id", gene_ids)
    create_string_dataset(features_group, "feature_type", feat_types)
    create_string_dataset(features_group, "_all_tag_keys", [])


def validate_barcodes(barcodes: np.ndarray | pd.Index | list[str]):
    bc_pattern = re.compile("^([ACTG]{6,})(-.*?)?$")
    return all([bc_pattern.match(bc) for bc in barcodes])


def format_barcodes(adata: AnnData) -> np.ndarray:
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


def clusters_to_h5(h5_obj: h5py.File, adata: AnnData) -> None:
    # clusters (aka all categorical metadata
    clusters_group = h5_obj.create_group("clusters")
    for colname in adata.obs_keys():
        dtype = adata.obs[colname].dtype
        if is_numeric_dtype(dtype):
            continue

        as_cat = adata.obs[colname].astype("category")

        group = clusters_group.create_group(colname)
        create_string_dataset(group, "name", [colname])
        create_string_dataset(group, "group_names", as_cat.cat.categories.values)
        group.create_dataset("assignments", data=as_cat.cat.codes.astype(np.int32))
        group.create_dataset("score", data=np.zeros(1, dtype=np.float64))
        create_string_dataset(group, "clustering_type", ["unknown"])


def projections_to_h5(h5_obj: h5py.File, adata: AnnData) -> None:
    projections_group = h5_obj.create_group("projections")
    for proj_name in adata.obsm_keys():
        proj = adata.obsm[proj_name].T
        is_umap = "umap" in proj_name
        is_tsne = ("tsne" in proj_name) or ("t-sne" in proj_name)

        method = proj_name
        if is_umap:
            method = "UMAP"
        if is_tsne:
            method = "t-SNE"

        group = projections_group.create_group(proj_name)
        create_string_dataset(group, "name", [proj_name])
        create_string_dataset(group, "method", [method])
        group.create_dataset("data", data=proj[:2, :])


def metadata_to_h5(h5_obj: h5py.File) -> None:
    info = sys.version_info
    version = "{}.{}.{}".format(info.major, info.minor, info.micro)

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
            "loupepy_h5py_version": h5py.version.version,
        },
    }

    dict_to_h5(h5_obj, metadata, "metadata")


def dict_to_h5(parent: h5py.Group, data: dict, groupname: str) -> None:
    group = parent.create_group(groupname)

    for key, val in data.items():
        if isinstance(val, dict):
            dict_to_h5(group, val, key)
        else:
            create_string_dataset(group, key, [val])


def create_string_dataset(
    group: h5py.Group, dset_name: str, strs: np.array, **kwargs
) -> None:
    if len(strs) < 1:
        longest = 1
    else:
        longest = max(map(len, strs))
    converted = np.array(strs).astype(f"S{longest}")
    group.create_dataset(dset_name, data=converted, **kwargs)


def read_h5(h5_path: Path | str):
    h5_path = Path(h5_path)
    assert h5_path.exists()

    with h5py.File(h5_path, "r") as h5:
        x = h5["/matrix/data"]
        i = h5["/matrix/indices"]
        p = h5["/matrix/indptr"]
        mat = csc_matrix((x, i, p)).T

        bcs = h5["/matrix/barcodes"]
        genes = h5["/matrix/features/name"]
        gids = h5["/matrix/features/id"]
        fts = h5["/matrix/features/feature_type"]
        var = pd.DataFrame(
            {"gene_ids": gids, "feature_type": fts}, index=pd.Index(genes)
        )

        adata = AnnData(mat, obs=pd.DataFrame(index=bcs), var=var)

        for key in h5["/projections"].keys():
            adata.obsm[f"{key}"] = np.array(h5[f"/projections/{key}/data"]).T

        for key in h5["/clusters"].keys():
            cats = np.array(list(map(bytes.decode, h5[f"/clusters/{key}/group_names"])))
            codes = np.array(h5[f"/clusters/{key}/assignments"])
            if (codes < 0).any():
                cats = np.append(cats, ["NULL"])
                nan_code = codes.max() + 1
                codes[codes < 0] = nan_code
            adata.obs[key] = cats[codes]
            adata.obs[key] = adata.obs[key].astype(pd.CategoricalDtype(categories=cats))

    return adata


def teardown_h5(h5path: Path | str):
    h5path = Path(h5path)
    with h5py.File(h5path, "r") as h5:
        print(h5["/"].keys())
        print(h5["/matrix"].keys())
        print(h5["/matrix/data"][:5])
        print(h5["/matrix/data"].dtype)
        print(h5["/matrix/indices"][:5])
        print(h5["/matrix/indices"].dtype)
        print(h5["/matrix/indptr"][:5])
        print(h5["/matrix/indptr"].dtype)
        print(h5["/matrix/features"].keys())
        print(h5["/matrix/features/id"][:5])
        print(h5["/matrix/features/name"][:5])
        print(h5["/matrix/features/feature_type"][:5])
        print(h5["/matrix/features/_all_tag_keys"])
        print(h5["/matrix/barcodes"])
        print(h5["/matrix/barcodes_unmodified"])
        print(h5["/projections"].keys())
        print(h5["/projections/umap"].keys())
        print(h5["/projections/umap/data"][:5])
        print(h5["/clusters"].keys())
        for key in h5["/clusters"].keys():
            print(h5[f"/clusters/{key}/assignments"][:5])
            print(h5[f"/clusters/{key}/name"])
            print(h5[f"/clusters/{key}/group_names"][:])
