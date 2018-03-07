@file:JvmName("VariantContextUtils")

package com.hartwig.hmftools.extensions.samtools

import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder

@JvmName("splitMultiAlleleVariant")
fun VariantContext.split(): List<VariantContext> {
    val ad = this.getGenotype(0).ad
    val dp = this.getGenotype(0).dp
    return this.alternateAlleles.mapIndexed { index, alt ->
        val alleles = listOf(this.reference, alt)
        val altADIndex = if (this.isWronglyPostProcessedIndel()) {
            1
        } else {
            index + 1
        }
        val genotype = GenotypeBuilder(this.sampleNamesOrderedByName[0], alleles).DP(dp)
                .AD(intArrayOf(ad[0], ad[altADIndex]))
                .make()
        VariantContextBuilder(this).alleles(alleles).genotypes(genotype).make()
    }
}

private fun VariantContext.isWronglyPostProcessedIndel(): Boolean {
    return this.isIndel && this.getGenotype(0).ad.size < this.alleles.size
}
