package com.hartwig.hmftools.extensions

import htsjdk.variant.variantcontext.Genotype

fun Genotype.allelicFrequency(isSingleBreakEnd: Boolean, isShort: Boolean): Double {
    val fragmentSupport = fragmentSupport(isSingleBreakEnd)
    val readPairSupport = if (isSingleBreakEnd || !isShort) this.refSupportReadPair() else 0
    val totalSupport = fragmentSupport + this.refSupportRead() + readPairSupport

    return fragmentSupport.toDouble() / totalSupport.toDouble();
}

fun Genotype.assemblyReadPairs(): Int = attributeAsInt("ASRP", 0);

fun Genotype.readPairs(): Int = attributeAsInt("RP", 0);

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
    return attributeAsInt(attribute, 0);
}

fun Genotype.attributeAsInt(attribute: String, defaultValue: Int): Int {
    val value = this.extendedAttributes[attribute]
    when (value) {
        null -> return defaultValue
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
