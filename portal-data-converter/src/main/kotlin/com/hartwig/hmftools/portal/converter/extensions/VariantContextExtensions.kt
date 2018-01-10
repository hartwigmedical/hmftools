package com.hartwig.hmftools.portal.converter.extensions

import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import org.apache.logging.log4j.LogManager

private val logger = LogManager.getLogger("VariantContextExtensions")


fun VariantContext.reformatAlleles(): VariantContext {
    assert(this.isBiallelic)
    val alleles = reformatAlleles(this.reference, this.alternateAlleles[0])
    return VariantContextBuilder(this)
            .loc(this.contig, this.start.toLong(), (this.start + alleles[0].length() - 1).toLong())
            .alleles(alleles).make()
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
        else -> {
            logger.warn("Determined mutation type 0 for variant ${this.contig}:${this.start}: " +
                    "${this.reference.baseString} -> ${this.alternateAlleles[0].baseString}}")
            "0"
        }
    }
}

//MIVO: detects variants like TAAAATAA -> TCGGAGGTCGCCGAAAATAA which occur in cases when there is both an insertion and deletion at same
//      base location
private fun isInsertionWithTrailingBases(ref: String, alt: String): Boolean {
    val refSubstring = ref.substring(1)
    val altSubstring = alt.substring(1)
    return alt.length > ref.length && ref[0] == alt[0] && altSubstring.lastIndexOf(refSubstring) == altSubstring.length - refSubstring.length
}