#!/bin/bash nextflow


process GetUncorrelatedVariants {
    tag "${bcf_file.simpleName}"

    input:
        path bcf_file

    output:
        path "corrective_variants.prune.in"

    script:
        """
        echo '6\t28510120\t33480577\tHLA\n' > hla_range.bed

        plink \
            --bcf ${bcf_file} \
            --out "corrective_variants"\
            --exclude 'range' hla_range.bed \
            --maf 0.05 \
            --geno 0.01 \
            --hwe 0.05 \
            --bp-space 100000 \
            --indep-pairwise 100 5 0.1

        """
}

process GetUncorrelatedVariantsRestricted {
    input:
        path bcf_file

    output:
        path "corrective_variants.prune.in"

    script:
        """
        echo '6\t28510120\t33480577\tHLA\n' > hla_range.bed

        plink \
            --bcf ${bcf_file} \
            --out "corrective_variants"\
            --exclude 'range' hla_range.bed \
            --geno 0.01 \
            --hwe 0.01 \
            --bp-space 500000 \
            --indep-pairwise 100 5 0.1

        """
}
