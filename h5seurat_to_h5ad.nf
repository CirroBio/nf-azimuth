process h5seurat_to_h5ad {
    publishDir params.output, mode: 'copy', overwrite: true, enabled: params.publish
    container "${params.container__azimuth}"
    input:
    path INPUT

    output:
    path "*.h5ad"

    script:
    """#!/usr/bin/env Rscript
library(Seurat)
library(SeuratData)
library(SeuratDisk)

input = "${INPUT}"
print(paste("Input:", input))
dest = gsub("h5seurat", "h5ad", input)
print(paste("Converting to h5ad:", dest))
Convert(input, dest=dest)
"""
}
