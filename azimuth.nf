process azimuth {
    publishDir params.output, mode: 'copy', overwrite: true
    container "${params.container__azimuth}"
    input:
    path INPUT
    val REFERENCE

    output:
    path "*.h5seurat"

    script:
    """#!/usr/bin/env Rscript
library(Azimuth)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tools)

az_ref = "${REFERENCE}"
input = "${INPUT}"

# Run Azimuth
print(paste("Analyzing input file", input))
print(paste("Using reference Azimuth dataset", az_ref))
res <- RunAzimuth(
    input,
    reference = az_ref
)

output_filename <- paste(file_path_sans_ext(input), az_ref, "h5seurat", sep=".")
print(paste("Saving to", output_filename))
SaveH5Seurat(res, filename = output_filename)
"""
}