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
    '''
    python2 accurate_pvalues.py \
    --ad-test-table ${andersonDarlingTable} \
    --empirical ${empiricalResults} \
    --permuted ${permutationResults} \
    --mafs ${minorAlleleFrequencies} \
    --breaks ${breakpoints} \
    --variants ${uncorrelatedVariants} \
    --output-prefix "output"
    '''
}

process GetUncorrelatedVariantsWithLdDataset {

    input:
      path ldDataset

    output:
      path "uncorrelated.tsv"

    script:
    '''

    '''

}

process GetBreakpoints
    input:
      path
