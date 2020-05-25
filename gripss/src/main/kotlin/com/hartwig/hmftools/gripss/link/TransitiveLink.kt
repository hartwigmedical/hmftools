package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import java.util.*

class TransitiveLink(private val assemblyLinkStore: LinkStore, private val variantStore: VariantStore) {

    companion object {
        const val MAX_TRANSITIVE_JUMPS = 2
        const val MAX_TRANSITIVE_DISTANCE = 1000
    }

    fun transitiveLink(variant: StructuralVariantContext, maxTransitiveJumps: Int = MAX_TRANSITIVE_JUMPS): List<Link> {

        if (!variant.isSingle) {
            val target = variantStore.select(variant.mateId!!)

            val alternativeStart = variantStore.selectAlternatives(variant)
            for (alternative in alternativeStart) {
                val assemblyLinks = findLinks("trs_${variant.vcfId}_", alternative, target, maxTransitiveJumps, mutableListOf())
                if (assemblyLinks.isNotEmpty()) {
                    return assemblyLinks
                }
            }
        }

        return Collections.emptyList()
    }

    private fun findLinks(transLinkPrefix: String, current: StructuralVariantContext, target: StructuralVariantContext, maxTransitiveJumps: Int, path: List<Link>)
            : List<Link> {
        val mate = variantStore.select(current.mateId!!)
        val newPath: List<Link> = path + Link(current)
        if (isMatch(mate, target)) {
            if (!target.imprecise) {
                val (minTotalDistance, maxTotalDistance) = newPath.totalDistance()
                if (target.insertSequenceLength < minTotalDistance || target.insertSequenceLength > maxTotalDistance) {
                    return Collections.emptyList()
                }
            }
            return newPath
        }

        // Always try assembled links first!
        val assemblyLinkedVariants = assemblyLinkStore.linkedVariants(current.mateId)
        for (link in assemblyLinkedVariants) {
            val linkedVariant = variantStore.select(link.otherVcfId)
            if (!linkedVariant.isSingle && !linkedVariant.imprecise) {
                val nextJump = Link(link.link, Pair(mate, linkedVariant))
                val newAssemblyLinks = findLinks(transLinkPrefix, linkedVariant, target, maxTransitiveJumps, newPath + nextJump)
                if (newAssemblyLinks.isNotEmpty()) {
                    return newAssemblyLinks
                }
            }
        }

        if (maxTransitiveJumps > 0) {
            val transitiveLinkedVariants = variantStore.selectTransitive(mate).filter { x -> !x.imprecise && !x.isSingle }
            for (linkedVariant in transitiveLinkedVariants) {
                val nextJump = Link("$transLinkPrefix${MAX_TRANSITIVE_JUMPS - maxTransitiveJumps}", Pair(mate, linkedVariant))
                val newAssemblyLinks = findLinks(transLinkPrefix, linkedVariant, target, maxTransitiveJumps - 1, newPath + nextJump)
                if (newAssemblyLinks.isNotEmpty()) {
                    return newAssemblyLinks
                }
            }
        }

        return Collections.emptyList()
    }


    private fun isMatch(current: StructuralVariantContext, target: StructuralVariantContext): Boolean {
        return current.vcfId != target.vcfId
                && current.mateId!! != target.vcfId
                && current.orientation == target.orientation
                && target.confidenceIntervalsOverlap(current)
    }


    private fun List<Link>.totalDistance(): Pair<Int, Int> {
        var minDistance = 0
        var maxDistance = 0

        for (variant in this) {
            minDistance += variant.minDistance
            maxDistance += variant.maxDistance
        }

        return Pair(minDistance, maxDistance)
    }

    private fun VariantStore.selectAlternatives(variant: StructuralVariantContext): List<StructuralVariantContext> {
        val maxDistance = if (variant.orientation == 1.toByte()) {
            Pair(0, variant.insertSequenceLength)
        } else {
            Pair(variant.insertSequenceLength, 0)
        }

        val alternativeFilter = { other: StructuralVariantContext -> other.orientation == variant.orientation && !other.imprecise && !other.isSingle }
        return selectOthersNearby(variant, maxDistance) { x -> alternativeFilter(x) }
    }

    private fun VariantStore.selectTransitive(variant: StructuralVariantContext): List<StructuralVariantContext> {
        val leftFilter = { other: StructuralVariantContext -> other.minStart <= variant.maxStart }
        val rightFilter = { other: StructuralVariantContext -> other.maxStart >= variant.minStart }
        val directionFilter = if (variant.orientation == 1.toByte()) leftFilter else rightFilter
        val transitiveFilter = { other: StructuralVariantContext -> other.orientation != variant.orientation && !other.imprecise && !other.isSingle }

        return selectOthersNearby(variant, MAX_TRANSITIVE_DISTANCE) { x -> directionFilter(x) && transitiveFilter(x) }
    }

}


