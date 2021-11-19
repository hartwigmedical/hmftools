package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import org.apache.logging.log4j.LogManager
import java.util.*

private typealias SvFilter = (StructuralVariantContext) -> Boolean

class TransitiveLink(private val assemblyLinkStore: LinkStore, private val variantStore: VariantStore) {

    companion object {
        private const val MAX_VARIANTS = 500_000
        private const val MAX_ALTERNATIVES = 25
        private const val MAX_ALTERNATIVES_SEEK_DISTANCE = 1000
        private const val MAX_ALTERNATIVES_ADDITIONAL_DISTANCE = MAX_ALTERNATIVES_SEEK_DISTANCE // This needs to be high to account for insert sequences

        private const val MAX_ASSEMBLY_JUMPS = 5
        private const val MAX_TRANSITIVE_JUMPS = 2

        private val logger = LogManager.getLogger(this::class.java)
    }

    fun transitiveLink(variant: StructuralVariantContext, maxAssemblyJumps: Int = MAX_ASSEMBLY_JUMPS, maxTransitiveJumps: Int = MAX_TRANSITIVE_JUMPS): List<Link> {
        if (variantStore.size() > MAX_VARIANTS) {
            return Collections.emptyList()
        }

        if (!variant.isSingle) {
            val target = variantStore.select(variant.mateId!!)
            val alternativeStarts = variantStore.selectAlternatives(variant)
            // Note: we really shouldn't expect to see that many variants within the CIPOS. If we do it is likely that it is a large
            // poly-G region or something equally messy.
            if (alternativeStarts.size in 1..MAX_ALTERNATIVES) {
                logger.debug("Examining ${alternativeStarts.size} alternative(s) to variant $target")
                val transLinkPrefix = "trs_${variant.vcfId}_"
                val assemblyNodes = ArrayDeque<Node>()

                for (alternative in alternativeStarts) {
                    val alternativeMate = variantStore.select(alternative.mateId!!)
                    val node = Node(transLinkPrefix, maxTransitiveJumps, maxAssemblyJumps, maxTransitiveJumps, alternative, alternativeMate, listOf(Link(alternative)))
                    assemblyNodes.add(node)
                }

                val assemblyLinks = findLinks(target, assemblyNodes, ArrayDeque(), ArrayDeque())
                if (assemblyLinks.isNotEmpty()) {
                    return assemblyLinks
                }
            }
        }

        return Collections.emptyList()
    }

    private tailrec fun findLinks(target: StructuralVariantContext, assemblyNodes: ArrayDeque<Node>, transitiveNodes: ArrayDeque<Node>, matchedTransitive: ArrayDeque<Node>): List<Link> {
        if (transitiveNodes.size > 1) {
            // No result if we there is more than one transitive path (and no assembly path)
            return Collections.emptyList()
        }

        if (assemblyNodes.isEmpty() && transitiveNodes.isEmpty()) {
            if (matchedTransitive.size == 1) {
                return matchedTransitive.pop().links
            }

            return Collections.emptyList()
        }

        if (assemblyNodes.isNotEmpty()) {
            val node = assemblyNodes.removeFirst()
            if (node.matchesTarget(target)) {
                // Return first (breath-wise) completely assembled link we see
                return node.links
            }

            node.assemblyNodes(assemblyLinkStore, variantStore).forEach { assemblyNodes.addLast(it) }
            node.transitiveNodes(assemblyLinkStore, variantStore).forEach { transitiveNodes.addLast(it) }
        } else {
            val node = transitiveNodes.removeFirst()
            if (node.matchesTarget(target)) {
                matchedTransitive.add(node)
            }
            node.assemblyNodes(assemblyLinkStore, variantStore).forEach { transitiveNodes.addLast(it) }
            node.transitiveNodes(assemblyLinkStore, variantStore).forEach { transitiveNodes.addLast(it) }
        }

        return findLinks(target, assemblyNodes, transitiveNodes, matchedTransitive)
    }

    private fun VariantStore.selectAlternatives(variant: StructuralVariantContext): Collection<StructuralVariantContext> {
        val isMatchingPositionAndOrientation: SvFilter = { other -> Node.isAlternative(variant, other) }
        val alternativeFilter: SvFilter = { other -> !other.imprecise && !other.isSingle }
        return selectOthersNearby(variant, MAX_ALTERNATIVES_ADDITIONAL_DISTANCE, MAX_ALTERNATIVES_SEEK_DISTANCE) { other -> isMatchingPositionAndOrientation(other) && alternativeFilter(other) }.sortByQualDesc()
    }
}

private data class Node(val transLinkPrefix: String, val maxTransitiveJumps: Int, val remainingAssemblyJumps: Int, val remainingTransitiveJumps: Int, val start: StructuralVariantContext, val end: StructuralVariantContext, val links: List<Link>) {

    companion object {
        private const val MIN_TRANSITIVE_DISTANCE = 30
        private const val MAX_TRANSITIVE_SEEK_DISTANCE = 2000
        private const val MAX_TRANSITIVE_ADDITIONAL_DISTANCE = 1000

        fun isAlternative(target: StructuralVariantContext, other: StructuralVariantContext, additionalAllowance: Int = 1): Boolean {
            if (target.orientation != other.orientation) {
                return false
            }

            val targetInsDistance = insertSequenceAdditionalDistance(target)
            val minStart = target.minStart - targetInsDistance.first - additionalAllowance
            val maxStart = target.maxStart + targetInsDistance.second + additionalAllowance

            val otherInsDistance = insertSequenceAdditionalDistance(other)
            val otherMinStart = other.minStart - otherInsDistance.first
            val otherMaxStart = other.maxStart + otherInsDistance.second

            return target.vcfId != other.vcfId
                    && target.mateId!! != other.vcfId
                    && otherMinStart <= maxStart
                    && otherMaxStart >= minStart
        }

        fun insertSequenceAdditionalDistance(variant: StructuralVariantContext): Pair<Int, Int> {
            return if (variant.orientation == 1.toByte()) {
                Pair(0, variant.insertSequenceLength)
            } else {
                Pair(variant.insertSequenceLength, 0)
            }
        }
    }

    fun isTransitive(): Boolean {
        return links.any { x -> x.link.contains("trs") }
    }

    fun matchesTarget(targetEnd: StructuralVariantContext): Boolean {
        if (isAlternative(targetEnd, end)) {
            if (!targetEnd.imprecise) {
                val targetDistance = targetEnd.insertSequenceLength + targetEnd.duplicationLength
                val (minTotalDistance, maxTotalDistance) = links.totalDistance()
                if (targetDistance < minTotalDistance || targetDistance > maxTotalDistance) {
                    return false
                }
            }
            return true
        }

        return false
    }

    fun assemblyNodes(assemblyLinkStore: LinkStore, variantStore: VariantStore): List<Node> {
        val result = mutableListOf<Node>()

        // Always try assembled links first!
        val unfilteredAssemblyLinks = assemblyLinkStore.linkedVariants(end.vcfId)
        if (remainingAssemblyJumps > 0) {
            val assemblyLinkedVariants = unfilteredAssemblyLinks
                    .filter { x -> !links.contains(x) }
                    .map { x -> Pair(x, variantStore.select(x.otherVcfId)) }
                    .sortedByDescending { x -> x.second.tumorQual }

            for (linkPair in assemblyLinkedVariants) {
                val link = linkPair.first
                val linkedVariant = linkPair.second
                if (!linkedVariant.isSingle && !linkedVariant.imprecise) {
                    val linkedVariantMate = variantStore.select(linkedVariant.mateId!!)
                    result.add(Node(transLinkPrefix, maxTransitiveJumps, remainingAssemblyJumps - 1, remainingTransitiveJumps, linkedVariant, linkedVariantMate, links + link + Link(linkedVariant)))
                }
            }
        }

        return result
    }

    fun transitiveNodes(assemblyLinkStore: LinkStore, variantStore: VariantStore): List<Node> {
        if (assemblyLinkStore.linkedVariants(end.vcfId).isNotEmpty()) {
            // Cannot have transitive links if there are assembly links
            return Collections.emptyList()
        }

        val result = mutableListOf<Node>()

        // unfilteredAssemblyLinks.isEmpty() &&
        if (remainingTransitiveJumps > 0) {
            val unlinkedFilter: SvFilter = { assemblyLinkStore.linkedVariants(it.vcfId).isEmpty() }
            val transitiveLinkedVariants = variantStore.selectTransitive(end).filter { x -> unlinkedFilter(x) }
            for (linkedVariant in transitiveLinkedVariants) {
                val linkedVariantMate = variantStore.select(linkedVariant.mateId!!)
                val link = Link("$transLinkPrefix${maxTransitiveJumps - remainingTransitiveJumps}", Pair(end, linkedVariant))
                result.add(Node(transLinkPrefix, maxTransitiveJumps, remainingAssemblyJumps, remainingTransitiveJumps - 1, linkedVariant, linkedVariantMate, links + link + Link(linkedVariant)))
            }
        }

        return result
    }

    private fun VariantStore.selectTransitive(variant: StructuralVariantContext): Collection<StructuralVariantContext> {
        val leftFilter: SvFilter = { other -> other.maxStart <= variant.minStart - MIN_TRANSITIVE_DISTANCE }
        val rightFilter: SvFilter = { other -> other.minStart >= variant.maxStart + MIN_TRANSITIVE_DISTANCE }
        val directionFilter: SvFilter = if (variant.orientation == 1.toByte()) leftFilter else rightFilter
        val transitiveFilter: SvFilter = { other -> other.orientation != variant.orientation && !other.imprecise && !other.isSingle }

        return selectOthersNearby(variant, MAX_TRANSITIVE_ADDITIONAL_DISTANCE, MAX_TRANSITIVE_SEEK_DISTANCE) { x -> directionFilter(x) && transitiveFilter(x) }.sortByQualDesc()
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
}

private fun Collection<StructuralVariantContext>.sortByQualDesc(): List<StructuralVariantContext> {
    return this.sortedByDescending { x -> x.tumorQual }
}
