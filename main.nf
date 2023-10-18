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

# Run Azimuth
print("Analyzing input file ${INPUT}")
print("Using reference Azimuth dataset ${params.reference}")
res <- RunAzimuth(
    "${INPUT}",
    reference = "${params.reference}"
)

print("Saving to ${params.reference}.h5seurat")
SaveH5Seurat(res, filename = "${params.reference}.h5seurat")
print("Converting to h5ad")
Convert("${params.reference}.h5seurat", dest = "h5ad")
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