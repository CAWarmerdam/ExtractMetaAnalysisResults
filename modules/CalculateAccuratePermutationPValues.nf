#!/bin/bash nextflow


process CalculateAccuratePermutationPValues {

    input:
      path empiricalResults
      path permutationResults
      path minorAlleleFrequencies
      path breakpoints
      path uncorrelatedVariants
      path andersonDarlingTable

    output:
      path "output"

    script:
    """
    accurate_pvalues.py \
    --ad-test-table ${andersonDarlingTable} \
    --empirical ${empiricalResults} \
    --permuted ${permutationResults} \
    --mafs ${minorAlleleFrequencies} \
    --breaks ${breakpoints} \
    --variants ${uncorrelatedVariants} \
    --output-prefix "output"
    """
}

process GetBreakpoints {

    input:
        path maf_table

    output:
        path "breaks.txt"

    script:
    """
    
    """

}
