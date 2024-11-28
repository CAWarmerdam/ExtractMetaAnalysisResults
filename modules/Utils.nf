#!/bin/bash nextflow


process CollectGzippedFiles {
    publishDir "${params.output}", mode: 'copy', overwrite: true, enabled: true

    input:
        path input
        val label

    output:
        path "${label}.out.txt.gz"

    shell:
        '''
        first=1
        for f in !{input.join(' ')}
        do
            if [ "$first" ]
            then
                gzcat "$f"
                first=
            else
                gzcat "$f"| tail -n +2
            fi
        done | gzip > "!{label}.collectgz.txt.gz"
        '''
}