#!/usr/bin/env python3

import json
import logging
from anndata import AnnData
import scanpy as sc
import pandas as pd
from vitessce import (
    VitessceConfig,
    Component as cm,
    AnnDataWrapper
)
from vitessce.data_utils import (
    optimize_adata,
    VAR_CHUNK_SIZE
)

sc.settings.verbosity = 3
sc.logging.print_header()

logger = logging.getLogger()


def vitessce():

    # Read the raw read counts
    input_fp = "${INPUT}"
    logger.info(f"Input file: {input_fp}")
    if input_fp.endswith(".h5ad"):
        logger.info("Parsing as h5ad")
        adata = sc.read_h5ad(input_fp)
    elif input_fp.endswith(".loom"):
        logger.info("Parsing as loom")
        adata = sc.read_loom(input_fp)
    elif input_fp.endswith(".h5"):
        logger.info("Parsing as 10x_h5")
        adata = sc.read_10x_h5(input_fp)
    else:
        logger.info("Did not recognize file extension, stopping")
        return

    logger.info("Making variable names unique")
    adata.var_names_make_unique()

    # Read the results from running Azimuth
    logger.info("Reading ${AZIMUTH} (as h5ad)")
    azimuth = sc.read_h5ad("${AZIMUTH}")

    # Apply the azimuth annotations and embeddings to the raw data
    logger.info("Adding azimuth annotations to the raw data")
    adata.obs = pd.merge(
        adata.obs,
        azimuth.obs,
        left_index=True,
        right_index=True
    )
    adata.obsm["X_umap"] = azimuth.obsm["X_ref.umap"]

    # annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # Basic QC metrics
    logger.info("Calculating QC metrics")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Normalize counts
    logger.info("Normalizing counts")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Annotate highly variable genes
    logger.info("Annotating highly variable genes")
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=300)

    assert "predicted.celltype.l1" in adata.obs.columns.values

    # Save the data to Zarr store
    logger.info("Saving to zarr")
    zarr_filepath = "azimuth.zarr"
    adata: AnnData = optimize_adata(
        adata,
        obs_cols=["predicted.celltype.l1"],
        obsm_keys=["X_umap"],
        optimize_X=True,
        var_cols=["highly_variable"],
    )
    adata.write_zarr(zarr_filepath, chunks=[adata.shape[0], VAR_CHUNK_SIZE])

    # Start to set up the vitessce viz
    logger.info("Setting up vitessce")
    vc = VitessceConfig(
        schema_version="1.0.15",
        name='${params.name}',
        description='${params.description}'
    )

    dataset = (
        vc
        .add_dataset(name='Azimuth')
        .add_object(
            AnnDataWrapper(
                adata_path=zarr_filepath,
                obs_embedding_paths=["obsm/X_umap"],
                obs_embedding_names=["UMAP"],
                obs_set_paths=["obs/predicted.celltype.l1"],
                obs_set_names=["Cell Type"],
                obs_feature_matrix_path="X",
                initial_feature_filter_path="var/highly_variable"
            )
        )
    )

    scatterplot = vc.add_view(cm.SCATTERPLOT, dataset=dataset, mapping="UMAP")
    cell_set_sizes = vc.add_view(cm.OBS_SET_SIZES, dataset=dataset)
    heatmap = vc.add_view(cm.HEATMAP, dataset=dataset)

    vc.layout((scatterplot | cell_set_sizes) / heatmap)

    vc_dict = vc.to_dict(base_url=".")
    for dataset in vc_dict["datasets"]:
        for file in dataset["files"]:
            file["url"] = zarr_filepath
    logger.info("Saving vitessce config")
    with open("vitessce.config.0.json", "w") as handle:
        json.dump(
            vc_dict,
            handle,
            indent=4
        )


if __name__ == "__main__":
    vitessce()
