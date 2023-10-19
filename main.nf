#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process azimuth {
    publishDir params.output, mode: 'copy', overwrite: true
    container "${params.container__azimuth}"
    input:
    path INPUT

    output:
    path "*.h5ad", emit: h5ad
    path "*.h5seurat", emit: h5seurat

    """#!/usr/bin/env Rscript
library(Azimuth)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

az_ref = "${params.reference}"
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
print("Converting to h5ad")
Convert(output_filename, dest = "h5ad")
"""
}

process vitessce {
    container "${params.container__vitessce}"
    publishDir "${params.web_output}", mode: 'copy', overwrite: true

    input:
        path INPUT
        path AZIMUTH

    output:
        path "*"

    script:
    template "vitessce.py"

}

workflow {
    if(!params.input){error "Must provide --input"}
    if(!params.output){error "Must provide --output"}

    input = file(params.input, checkIfExists: true)
    azimuth(input)

    if (params.web_output){
        vitessce(
            input,
            azimuth.out.h5ad
        )
    }
}