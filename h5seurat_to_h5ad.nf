process h5seurat_to_h5ad {
    publishDir params.output, mode: 'copy', overwrite: true, enabled: params.publish
    container "${params.container__azimuth}"
    input:
    path INPUT

    output:
    path "*.h5ad"

    """#!/usr/bin/env Rscript
library(Seurat)
library(SeuratData)
library(SeuratDisk)

input = "${INPUT}"
print(paste("Input:", input))
print("Converting to h5ad")
Convert(input, dest = "h5ad")
"""
}
