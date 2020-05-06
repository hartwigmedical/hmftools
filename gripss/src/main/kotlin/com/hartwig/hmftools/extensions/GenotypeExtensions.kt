package com.hartwig.hmftools.extensions

import htsjdk.variant.variantcontext.Genotype

fun Genotype.refSupportRead(): Int = requireAttributeAsInt("REF");

fun Genotype.refSupportReadPair(): Int = requireAttributeAsInt("REFPAIR");

fun Genotype.fragmentSupport(isSingleBreakend: Boolean): Int {
    return if (isSingleBreakend) {
        requireAttributeAsInt("BVF")
    } else {
        requireAttributeAsInt("VF")
    }
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
