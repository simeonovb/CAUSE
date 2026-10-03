nextflow.enable.dsl=2

process sayHello {
    input:
    val x

    output:
    path "output.txt"

    script:
    """
    echo "${x} World!" > output.txt
    """
}

workflow {
    input_ch = Channel.of(params.greeting)
    sayHello(input_ch)
}
