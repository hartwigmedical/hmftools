package com.hartwig.hmftools.gripsskt

import htsjdk.variant.variantcontext.Genotype

fun Genotype.qual(isSingleBreakEnd: Boolean, isLineInsertion: Boolean): Double {
    if(isSingleBreakEnd)
    {
        return if (isLineInsertion) attributeAsDouble("BQ", 0.0) else attributeAsDouble("BAQ", 0.0)
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

    return if (totalSupport == 0) 0.0 else fragmentSupport.toDouble() / totalSupport.toDouble()
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
        var bscCount = attributeAsInt("BSC", 0)
        var basrpCount = attributeAsInt("BASRP", 0)
        var bassrCount = attributeAsInt("BASSR", 0)

        var bvf = requireAttributeAsInt("BVF")

        return if (bscCount == 0 && basrpCount == 0 && bassrCount == 0) 0 else bvf
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
