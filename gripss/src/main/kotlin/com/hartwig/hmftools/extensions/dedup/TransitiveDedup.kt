package com.hartwig.hmftools.extensions.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.link.Link
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import org.apache.logging.log4j.LogManager
import java.util.*

class TransitiveDedup(private val assemblyLinkStore: LinkStore, private val variantStore: VariantStore) {

    companion object {
        const val MAX_TRANSITIVE_JUMPS = 2
        private val logger = LogManager.getLogger(this::class.java)
    }

    fun dedup(variant: StructuralVariantContext): Boolean {

        if (!variant.isSingle) {
            val target = variantStore.select(variant.mateId!!)

            val alternativeStart = variantStore.selectAlternatives(variant).filter { x -> !x.imprecise && !x.isSingle }
            for (alternative in alternativeStart) {
                val assemblyLinks = assemblyLinks("trs_${variant.vcfId}_", alternative, target, MAX_TRANSITIVE_JUMPS, mutableListOf())
                if (assemblyLinks.isNotEmpty()) {
                    logger.info("Found alternate mapping of $variant CIPOS:${variant.confidenceInterval} IMPRECISE:${variant.imprecise} -> ${assemblyLinks.alternatePath()}")
                    return true
                }
            }
        }

        return false
    }


    private fun assemblyLinks(transLinkPrefix: String, current: StructuralVariantContext, target: StructuralVariantContext, maxTransitiveJumps: Int, path: List<Link>): List<Link> {
        val mate = variantStore.select(current.mateId!!)
        val newPath: List<Link> = path + Link(current)
        if (matchTarget(mate, target)) {
            if (!target.imprecise) {
                val (minTotalDistance, maxTotalDistance) = newPath.totalDistance()
                if (target.insertSequenceLength < minTotalDistance || target.insertSequenceLength > maxTotalDistance) {
                    return Collections.emptyList()
                }
            }
            return newPath
        }

        // Always try assembled links first!
        val assemblyLinkedVariants = assemblyLinkStore.linkedVariants(current.mateId).map { vcfId -> variantStore.select(vcfId) }.filter { x -> !x.imprecise && !x.isSingle }
        for (linkedVariant in assemblyLinkedVariants) {
            val nextJump = Link("ASM", Pair(mate, linkedVariant))
            val newAssemblyLinks = assemblyLinks(transLinkPrefix, linkedVariant, target, maxTransitiveJumps, newPath + nextJump)
            if (newAssemblyLinks.isNotEmpty()) {
                return newAssemblyLinks
            }
        }

        if (maxTransitiveJumps > 0) {
            val transitiveLinkedVariants = variantStore.selectNearbyFacingVariants(mate).filter { x -> !x.imprecise && !x.isSingle }
            for (linkedVariant in transitiveLinkedVariants) {
                val nextJump = Link("$transLinkPrefix${MAX_TRANSITIVE_JUMPS - maxTransitiveJumps}", Pair(mate, linkedVariant))
                val newAssemblyLinks = assemblyLinks(transLinkPrefix, linkedVariant, target, maxTransitiveJumps - 1, newPath + nextJump)
                if (newAssemblyLinks.isNotEmpty()) {
                    return newAssemblyLinks
                }
            }
        }

        return Collections.emptyList()
    }


    private fun matchTarget(current: StructuralVariantContext, target: StructuralVariantContext): Boolean {
        return current.vcfId != target.vcfId && current.mateId!! != target.vcfId && current.orientation == target.orientation && target.confidenceIntervalsOverlap(current)
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

    private fun List<Link>.alternatePath(): String {
        val stringJoiner = StringJoiner("")

        for (i in this.indices) {
            if (i == 0) {
                stringJoiner.add(this[i].vcfId)
            }

            val link = this[i]
            if (link.link == "PAIR") {
                stringJoiner.add("-")
            } else {
                stringJoiner.add("<${link.link}>")
            }

            stringJoiner.add(link.otherVcfId)
        }

        return stringJoiner.toString()
    }
}


