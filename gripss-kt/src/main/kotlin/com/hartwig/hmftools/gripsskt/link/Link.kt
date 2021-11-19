package com.hartwig.hmftools.gripsskt.link

import com.hartwig.hmftools.gripsskt.StructuralVariantContext
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

data class Link(val link: String, val vcfId: String, val otherVcfId: String, val minDistance: Int, val maxDistance: Int) {

    companion object {
        operator fun invoke(variant: StructuralVariantContext): Link {
            val distance = variant.duplicationLength + variant.insertSequenceLength
            return Link("PAIR", variant.vcfId, variant.mateId!!, distance, distance)
        }

        operator fun invoke(link: String, variants: Pair<StructuralVariantContext, StructuralVariantContext>): Link {
            val (minDistance, maxDistance) = distance(variants.first, variants.second)
            return Link(link, variants.first.vcfId, variants.second.vcfId, minDistance, maxDistance)
        }


        private fun distance(first: StructuralVariantContext, second: StructuralVariantContext): Pair<Int, Int> {
            val minDistance = abs(first.maxStart - second.minStart)
            val maxDistance = abs(first.minStart - second.maxStart)
            return Pair(min(minDistance, maxDistance), max(minDistance, maxDistance))
        }
    }

    fun reverse(): Link {
        return Link(link, otherVcfId, vcfId, minDistance, maxDistance)
    }

    override fun toString(): String {
        return "$vcfId<$link>$otherVcfId"
    }
}