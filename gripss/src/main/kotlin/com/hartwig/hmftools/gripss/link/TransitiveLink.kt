package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import org.apache.logging.log4j.LogManager
import java.util.*

class TransitiveLink(private val assemblyLinkStore: LinkStore, private val variantStore: VariantStore) {

    companion object {
        private const val MAX_ALTERNATIVES = 25
        private const val MAX_ALTERNATIVES_SEEK_DISTANCE = 1000
        private const val MAX_ALTERNATIVES_ADDITIONAL_DISTANCE = MAX_ALTERNATIVES_SEEK_DISTANCE // This needs to be high to account for insert sequences

        private const val MAX_ASSEMBLY_JUMPS = 5
        private const val MAX_TRANSITIVE_JUMPS = 2
        private const val MAX_TRANSITIVE_SEEK_DISTANCE = 2000
        private const val MAX_TRANSITIVE_ADDITIONAL_DISTANCE = 1000

        private val logger = LogManager.getLogger(this::class.java)
    }

    fun transitiveLink(variant: StructuralVariantContext, maxAssemblyJumps: Int = MAX_ASSEMBLY_JUMPS, maxTransitiveJumps: Int = MAX_TRANSITIVE_JUMPS): List<Link> {

        if (!variant.isSingle) {
            val target = variantStore.select(variant.mateId!!)
            val alternativeStart = variantStore.selectAlternatives(variant)
            // Note: we really shouldn't expect to see that many variants within the CIPOS. If we do it is likely that it is a large
            // poly-G region or something equally messy.
            if (alternativeStart.size in 1..MAX_ALTERNATIVES) {
                logger.debug("Examining ${alternativeStart.size} alternative(s) to variant $target")
                for (alternative in alternativeStart) {
                    val assemblyLinks = findLinks("trs_${variant.vcfId}_", alternative, target, maxAssemblyJumps, maxTransitiveJumps, mutableListOf())
                    if (assemblyLinks.isNotEmpty()) {
                        return assemblyLinks
                    }
                }
            }
        }

        return Collections.emptyList()
    }

    private fun findLinks(transLinkPrefix: String, current: StructuralVariantContext, target: StructuralVariantContext, maxAssemblyJumps: Int, maxTransitiveJumps: Int, path: List<Link>)
            : List<Link> {
        val mate = variantStore.select(current.mateId!!)
        val newPath: List<Link> = path + Link(current)
        if (isAlternative(target, mate)) {
            if (!target.imprecise) {
                val targetDistance = target.insertSequenceLength + target.duplicationLength
                val (minTotalDistance, maxTotalDistance) = newPath.totalDistance()
                if (targetDistance < minTotalDistance || targetDistance > maxTotalDistance) {
                    return Collections.emptyList()
                }
            }
            return newPath
        }

        // Always try assembled links first!
        if (maxAssemblyJumps > 0) {
            val assemblyLinkedVariants = assemblyLinkStore.linkedVariants(current.mateId).filter { x -> !path.contains(x) }
            for (link in assemblyLinkedVariants) {
                val linkedVariant = variantStore.select(link.otherVcfId)
                if (!linkedVariant.isSingle && !linkedVariant.imprecise) {
                    val newAssemblyLinks = findLinks(transLinkPrefix, linkedVariant, target, maxAssemblyJumps - 1, maxTransitiveJumps, newPath + link)
                    if (newAssemblyLinks.isNotEmpty()) {
                        return newAssemblyLinks
                    }
                }
            }
        }

        if (maxTransitiveJumps > 0) {
            val transitiveLinkedVariants = variantStore.selectTransitive(mate).filter { x -> !x.imprecise && !x.isSingle && assemblyLinkStore.linkedVariants(x.vcfId).isEmpty() }
            for (linkedVariant in transitiveLinkedVariants) {
                val nextJump = Link("$transLinkPrefix${MAX_TRANSITIVE_JUMPS - maxTransitiveJumps}", Pair(mate, linkedVariant))
                val newAssemblyLinks = findLinks(transLinkPrefix, linkedVariant, target, maxAssemblyJumps, maxTransitiveJumps - 1, newPath + nextJump)
                if (newAssemblyLinks.isNotEmpty()) {
                    return newAssemblyLinks
                }
            }
        }

        return Collections.emptyList()
    }


    private fun insertSequenceAdditionalDistance(variant: StructuralVariantContext): Pair<Int, Int> {
        return if (variant.orientation == 1.toByte()) {
            Pair(0, variant.insertSequenceLength)
        } else {
            Pair(variant.insertSequenceLength, 0)
        }
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

    private fun isAlternative(target: StructuralVariantContext, other: StructuralVariantContext): Boolean {
        return target.orientation == other.orientation && isMatchingPosition(target, other)
    }

    private fun VariantStore.selectAlternatives(variant: StructuralVariantContext): Collection<StructuralVariantContext> {
        val isMatchingPositionAndOrientation = { other: StructuralVariantContext -> isAlternative(variant, other) }
        val alternativeFilter = { other: StructuralVariantContext -> !other.imprecise && !other.isSingle }
        return selectOthersNearby(variant, MAX_ALTERNATIVES_ADDITIONAL_DISTANCE, MAX_ALTERNATIVES_SEEK_DISTANCE) { other -> isMatchingPositionAndOrientation(other) && alternativeFilter(other) }
    }

    private fun VariantStore.selectTransitive(variant: StructuralVariantContext): Collection<StructuralVariantContext> {
        val leftFilter = { other: StructuralVariantContext -> other.minStart <= variant.maxStart }
        val rightFilter = { other: StructuralVariantContext -> other.maxStart >= variant.minStart }
        val directionFilter = if (variant.orientation == 1.toByte()) leftFilter else rightFilter
        val transitiveFilter = { other: StructuralVariantContext -> other.orientation != variant.orientation && !other.imprecise && !other.isSingle }

        return selectOthersNearby(variant, MAX_TRANSITIVE_ADDITIONAL_DISTANCE, MAX_TRANSITIVE_SEEK_DISTANCE) { x -> directionFilter(x) && transitiveFilter(x) }
    }

    private fun isMatchingPosition(target: StructuralVariantContext, other: StructuralVariantContext): Boolean {
        val targetInsDistance = insertSequenceAdditionalDistance(target)
        val minStart = target.minStart - targetInsDistance.first
        val maxStart = target.maxStart + targetInsDistance.second

        val otherInsDistance = insertSequenceAdditionalDistance(other);
        val otherMinStart = other.minStart - otherInsDistance.first
        val otherMaxStart = other.maxStart + otherInsDistance.second

        return target.vcfId != other.vcfId
                && target.mateId!! != other.vcfId
                && otherMinStart <= maxStart
                && otherMaxStart >= minStart
    }

}


