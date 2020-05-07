package com.hartwig.hmftools.extensions

import htsjdk.variant.variantcontext.VariantContext
import kotlin.math.max

fun VariantContext.strandBias(): Double {
    val strandBias = this.getAttributeAsDouble("SB", 0.5)
    return max(strandBias, 1 - strandBias)
}

fun VariantContext.imprecise(): Boolean {
    return this.getAttributeAsBoolean("IMPRECISE", false)
}

fun VariantContext.homLength(): Int {
    return this.getAttributeAsInt("HOMLEN", 0)
}