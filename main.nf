#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process azimuth {
    publishDir params.output, mode: 'copy', overwrite: true
    input:
    path INPUT

    output:
    path "*.h5*"

    """#!/usr/bin/env Rscript

# Run Azimuth
print("Analyzing input file ${INPUT}")
print("Using reference Azimuth dataset ${params.reference}")
res <- RunAzimuth(
    "${INPUT}",
    reference = "${params.reference}"
)

print("Saving to ${params.reference}.h5seurat")
SaveH5Seurat(res, filename = "${params.reference}.h5Seurat")
print("Converting to h5ad")
Convert("${params.reference}.h5Seurat", dest = "h5ad")
"""
}

workflow {
    Channel
        .fromPath(
            params.input.split(',').toList(),
            checkIfExists: true
        )
        .flatten()
        | azimuth
}