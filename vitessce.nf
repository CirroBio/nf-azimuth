process vitessce {
    container "${params.container__vitessce}"
    publishDir "${params.web_output}", mode: 'copy', overwrite: true

    input:
        path INPUT
        path AZIMUTH

    output:
        path "vitessce.config.*.json", optional: true
        path "chart.manifest.json", optional: true
        path "azimuth.zarr/**", optional: true, hidden: true

    script:
    template "vitessce.py"

}