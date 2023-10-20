process vitessce {
    container "${params.container__vitessce}"
    publishDir "${params.web_output}", mode: 'copy', overwrite: true

    input:
        path INPUT
        path AZIMUTH

    output:
        path "*", optional: true, hidden: true

    script:
    template "vitessce.py"

}