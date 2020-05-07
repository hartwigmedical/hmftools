package com.hartwig.hmftools.extensions

import htsjdk.variant.variantcontext.Genotype

fun Genotype.allelicFrequency(isSingleBreakEnd: Boolean, isShortDelDup: Boolean): Double {
    val fragmentSupport = fragmentSupport(isSingleBreakEnd)
    val readPairSupport = if (isSingleBreakEnd || !isShortDelDup) this.refSupportReadPair() else 0
    val totalSupport = fragmentSupport + this.refSupportRead() + readPairSupport

    return fragmentSupport.toDouble() / totalSupport.toDouble();
}

fun Genotype.splitRead(): Int = requireAttributeAsInt("SR");

fun Genotype.refSupportRead(): Int = requireAttributeAsInt("REF");

fun Genotype.refSupportReadPair(): Int = requireAttributeAsInt("REFPAIR");

fun Genotype.fragmentSupport(isSingleBreakEnd: Boolean): Int {
    val attribute = when (isSingleBreakEnd) {
        true -> "BVF"
        false -> "VF"
    }
    return requireAttributeAsInt(attribute)
}

fun Genotype.requireAttributeAsInt(attribute: String): Int {
    this.checkAttribute(attribute)

    val value = this.extendedAttributes[attribute]
    when (value) {
        is Int -> return value
        is String -> return value.toInt()
    }

    throw IllegalStateException("Unable to interpret ${value} as Int")
}

private fun Genotype.checkAttribute(attribute: String) {
    if (!this.hasExtendedAttribute(attribute)) {
        throw IllegalStateException("Genotype must contain attribute $attribute")
    }
}
