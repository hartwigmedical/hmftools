package com.hartwig.hmftools.portal.converter.extensions

import com.google.common.collect.Lists
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder

fun VariantContext.split(): List<VariantContext> {
    val variants = Lists.newArrayList<VariantContext>()
    val ad = this.getGenotype(0).ad
    val dp = this.getGenotype(0).dp
    for (index in 0 until this.alternateAlleles.size) {
        val alt = this.getAlternateAllele(index)
        val alleles = reformatAlleles(this.reference, alt)
        val altADIndex = if (this.isWronglyPostProcessedIndel()) {
            1
        } else {
            index + 1
        }
        val genotype = GenotypeBuilder(this.sampleNamesOrderedByName[0], alleles).DP(dp)
                .AD(intArrayOf(ad[0], ad[altADIndex]))
                .make()
        variants.add(VariantContextBuilder(this)
                .loc(this.contig, this.start.toLong(), (this.start + alleles[0].length() - 1).toLong()).alleles(alleles)
                .genotypes(genotype).make())
    }
    return variants
}

private fun VariantContext.isWronglyPostProcessedIndel(): Boolean {
    return this.isIndel && this.getGenotype(0).ad.size < this.alleles.size
}

//MIVO: transforms variants like TAAAATAA -> TCGGAGGTCGCCGAAAATAA to regular format: T -> TCGGAGGTCGCCG
private fun reformatAlleles(ref: Allele, alt: Allele): List<Allele> {
    if (isInsertionWithTrailingBases(ref.baseString, alt.baseString)) {
        val newRef = Allele.create(ref.baseString[0].toString(), true)
        val lastIndexOfRefBasesInAlt = alt.baseString.lastIndexOf(ref.baseString.substring(1))
        val newAlt = Allele.create(alt.baseString.substring(0, lastIndexOfRefBasesInAlt))
        return listOf(newRef, newAlt)
    }
    return listOf(ref, alt)
}

fun VariantContext.mutationType(): String {
    return when {
        this.isSNP -> "1"
        this.isSimpleInsertion -> "2"
        this.isSimpleDeletion -> "3"
        this.isMNP -> "4"
        else -> "0"
    }
}

//MIVO: detects variants like TAAAATAA -> TCGGAGGTCGCCGAAAATAA which occur in cases when there is both an insertion and deletion at same
//      base location
private fun isInsertionWithTrailingBases(ref: String, alt: String): Boolean {
    val refSubstring = ref.substring(1)
    val altSubstring = alt.substring(1)
    return alt.length > ref.length && ref[0] == alt[0] && altSubstring.lastIndexOf(refSubstring) == altSubstring.length - refSubstring.length
}