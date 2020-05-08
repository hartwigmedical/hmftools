package com.hartwig.hmftools.extensions

import com.hartwig.hmftools.gripss.VariantType
import htsjdk.variant.variantcontext.VariantContext
import kotlin.math.max

fun VariantContext.strandBias(): Double {
    val strandBias = this.getAttributeAsDouble("SB", 0.5)
    return max(strandBias, 1 - strandBias)
}

fun VariantContext.imprecise(): Boolean {
    return this.getAttributeAsBoolean("IMPRECISE", false)
}

fun VariantContext.homologyLength(): Int {
    return this.getAttributeAsInt("HOMLEN", 0)
}

fun VariantContext.inexactHomologyLength(): Int {
    if (!this.hasAttribute("IHOMPOS")) {
        return 0
    }

    val (iHom1, iHom2) = this.getAttributeAsIntList("IHOMPOS", 0)
    return max(0, iHom2 - iHom1)
}

fun VariantContext.toVariantType(): VariantType {
    return VariantType.create(this.contig, this.start, this.alleles[0].displayString, this.alleles[1].displayString)
}