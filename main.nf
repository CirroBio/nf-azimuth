#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

include { h5seurat_to_h5ad as input_to_h5ad } from "./h5seurat_to_h5ad" addParams(publish: false)
include { h5seurat_to_h5ad as output_to_h5ad } from "./h5seurat_to_h5ad" addParams(publish: true)
include { azimuth } from "./azimuth"
include { vitessce } from "./vitessce"

workflow {
    if(!params.input){error "Must provide --input"}
    if(!params.output){error "Must provide --output"}

    input = file(params.input, checkIfExists: true)
    azimuth(input)
    output_to_h5ad(azimuth.out)

    if (params.web_output){

        if (params.input =~ /h5seurat/){
            input_to_h5ad(input)
            vitessce(
                input_to_h5ad.out,
                output_to_h5ad.out
            )
        } else {
            vitessce(
                input,
                output_to_h5ad.out
            )
        }
    }
}