package com.hartwig.hmftools.gripss

import htsjdk.variant.variantcontext.Genotype

fun Genotype.qual(isSingleBreakEnd: Boolean): Double {
    if(isSingleBreakEnd)
    {
        var bqQual = requireAttributeAsDouble("BQ")
        var bumQual = attributeAsDouble("BUMQ", 0.0)
        return bqQual - bumQual
    }
    else
    {
        return requireAttributeAsDouble("QUAL")
    }
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

    if(isSingleBreakEnd)
    {
        var bumCount = attributeAsInt("BUM", 0)
        var bvf = requireAttributeAsInt("BVF")

        return if (bvf == bumCount) 0 else bvf
    }
    else
    {
        return requireAttributeAsInt("VF")
    }
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
