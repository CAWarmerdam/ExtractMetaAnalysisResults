
process GlobalPrincipalExpressionPca {

    input:
        path permuted
        val nThreshold
        val nPcs

    output:
        path "global_expression_components.rds"
        path "uncorrelated_genes.rds"

    script:
    """
    Rscript --vanilla $baseDir/bin/global_expression_pcs.R \
    --permuted ${permuted} \
    --n-threshold ${nThreshold} \
    --n-pcs ${nPcs}
    """
}

process EstimateComponentEffectsOnVariants {
    input:
        path empirical
        path global_expression_components
        path gene_reference
        path variant_reference
        val variants
        val uncorrelated_genes

    output:
        path "component_variant_estimates.rds"

    script:
    """
    Rscript --vanilla $baseDir/bin/remove_global_pcs.R \
    --empirical ${empirical} \
    --eigenvectors ${global_expression_components} \
    --gene-reference ${gene_reference} \
    --variant-reference ${variant_reference} \
    --variant-selection ${variants} \
    --uncorrelated-genes ${uncorrelated_genes}
    """
}

process CorrectQtlEffects {
    input:
        path empirical
        path global_expression_components
        path component_variant_estimates
        val genes
        val unbiased_variants

    output:
        path "corrected_eqtls"
        path "residual_variance_explained.rds"

    script:
    Rscript --vanilla $baseDir/bin/remove_global_pcs.R \
    --empirical ${empirical} \
    --eigenvectors ${global_expression_components} \
    --component-variant-estimates ${component_variant_estimates} \
    --genes ${genes} \
    --unbiased-variants ${unbiased_variants}
}