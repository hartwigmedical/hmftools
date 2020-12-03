package com.hartwig.hmftools.extensions

import htsjdk.variant.variantcontext.Genotype

fun Genotype.qual(isSingleBreakEnd: Boolean): Double {
    return if (isSingleBreakEnd) requireAttributeAsDouble("BQ") else requireAttributeAsDouble("QUAL")
}

fun Genotype.allelicFrequency(isSingleBreakEnd: Boolean, isShort: Boolean): Double {
    val fragmentSupport = fragmentSupport(isSingleBreakEnd)
    val readPairSupport = if (isSingleBreakEnd || !isShort) this.refSupportReadPair() else 0
    val totalSupport = fragmentSupport + this.refSupportRead() + readPairSupport

    return fragmentSupport.toDouble() / totalSupport.toDouble()
}

fun Genotype.assemblyReadPairs(): Int = attributeAsInt("ASRP", 0)

fun Genotype.readPairs(): Int = attributeAsInt("RP", 0)

fun Genotype.splitRead(): Int = requireAttributeAsInt("SR")

fun Genotype.indelCount(): Int = requireAttributeAsInt("IC")

fun Genotype.refSupportRead(): Int = requireAttributeAsInt("REF")

fun Genotype.refSupportReadPair(): Int = requireAttributeAsInt("REFPAIR")

fun Genotype.fragmentSupport(isSingleBreakEnd: Boolean): Int {
    val attribute = when (isSingleBreakEnd) {
        true -> "BVF"
        false -> "VF"
    }
    return requireAttributeAsInt(attribute)
}

fun Genotype.requireAttributeAsDouble(attribute: String): Double {
    this.checkAttribute(attribute)
    return attributeAsDouble(attribute, 0.0)
}

fun Genotype.requireAttributeAsInt(attribute: String): Int {
    this.checkAttribute(attribute)
    return attributeAsInt(attribute, 0)
}

fun Genotype.attributeAsInt(attribute: String, defaultValue: Int): Int {
    val value = this.extendedAttributes[attribute]
    when (value) {
        null -> return defaultValue
        is Int -> return value
        is String -> return value.toInt()
    }

    throw IllegalStateException("Unable to interpret $value as Int")
}

fun Genotype.attributeAsDouble(attribute: String, defaultValue: Double): Double {
    val value = this.extendedAttributes[attribute]
    when (value) {
        null -> return defaultValue
        is Double -> return value
        is Int -> return value.toDouble()
        is String -> return value.toDouble()
    }

    throw IllegalStateException("Unable to interpret $value as Int")
}

private fun Genotype.checkAttribute(attribute: String) {
    if (!this.hasExtendedAttribute(attribute)) {
        throw IllegalStateException("Genotype must contain attribute $attribute")
    }
}
