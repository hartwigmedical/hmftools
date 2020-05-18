package com.hartwig.hmftools.extensions.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import java.util.*
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

class TransitiveDedup(private val linkStore: LinkStore, private val variantStore: VariantStore) {

    fun dedup(variant: StructuralVariantContext): Boolean {

        if (!variant.isSingle) {
            val target = variantStore.select(variant.mateId!!)

            val alternativeStart = variantStore.selectAlternatives(variant)
            for (alternative in alternativeStart) {
                if (!alternative.isSingle) {
                    val (minDistance, maxDistance) = distance(variant, alternative)
                    val nextJump = LinkedVariantStart(alternative.vcfId, alternative.mateId!!, minDistance, maxDistance, alternative.insertSequenceLength)
                    val assemblyLinks = assemblyLinks(nextJump, target, 2, mutableListOf())
                    if (assemblyLinks.isNotEmpty()) {
                        println("$variant CIPOS:${variant.confidenceInterval} IMPRECISE:${variant.imprecise} -> ${assemblyLinks.joinToString("")}")
                        return true
                    }
                }
            }
        }

        return false
    }


    private fun assemblyLinks(current: LinkedVariant, target: StructuralVariantContext, maxTransitiveJumps: Int, path: List<LinkedVariant>): List<LinkedVariant> {
        val mate = variantStore.select(current.mateId!!)
        val newPath = path + current
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
        val assemblyLinkedVariants = linkStore.followLinks(current.mateId).map { vcfId -> variantStore.select(vcfId) }
        for (linkedVariant in assemblyLinkedVariants) {
            if (!linkedVariant.isSingle && !linkedVariant.imprecise) {
                val (minDistance, maxDistance) = distance(linkedVariant, mate)
                val nextJump = AssemblyLinkedVariant(linkedVariant.vcfId, linkedVariant.mateId!!, minDistance, maxDistance, linkedVariant.insertSequenceLength)
                val newAssemblyLinks = assemblyLinks(nextJump, target, maxTransitiveJumps, newPath)
                if (newAssemblyLinks.isNotEmpty()) {
                    return newAssemblyLinks
                }
            }
        }

        if (maxTransitiveJumps > 0) {
            val transitiveLinkedVariants = variantStore.selectTransitivelyLinkedVariants(mate)
            for (linkedVariant in transitiveLinkedVariants) {
                if (!linkedVariant.isSingle && !linkedVariant.imprecise) {
                    val (minDistance, maxDistance) = distance(linkedVariant, mate)
                    val nextJump = TransitiveLinkedVariant(linkedVariant.vcfId, linkedVariant.mateId!!, minDistance, maxDistance, linkedVariant.insertSequenceLength)
                    val newAssemblyLinks = assemblyLinks(nextJump, target, maxTransitiveJumps - 1, newPath)
                    if (newAssemblyLinks.isNotEmpty()) {
                        return newAssemblyLinks
                    }
                }
            }
        }

        return Collections.emptyList()
    }


    private fun matchTarget(current: StructuralVariantContext, target: StructuralVariantContext): Boolean {
        return current.vcfId != target.vcfId && current.mateId!! != target.vcfId && current.orientation == target.orientation && target.confidenceIntervalsOverlap(current)
    }

    private fun distance(first: StructuralVariantContext, second: StructuralVariantContext): Pair<Int, Int> {
        val minDistance = abs(first.maxStart - second.minStart)
        val maxDistance = abs(first.minStart - second.maxStart)
        return Pair(min(minDistance, maxDistance), max(minDistance, maxDistance))
    }

    private fun List<LinkedVariant>.totalDistance(): Pair<Int, Int> {

        var minDistance = 0
        var maxDistance = 0

        for (variant in this) {
            minDistance += variant.minDistance + variant.insertSequenceLength
            maxDistance += variant.maxDistance + variant.insertSequenceLength
        }

        return Pair(minDistance, maxDistance)
    }
}


private sealed class LinkedVariant {
    abstract val vcfId: String
    abstract val mateId: String
    abstract val insertSequenceLength: Int
    abstract val minDistance: Int
    abstract val maxDistance: Int
    override fun toString(): String {
        return "$vcfId-$mateId"
    }
}

private data class LinkedVariantStart(
        override val vcfId: String,
        override val mateId: String,
        override val minDistance: Int,
        override val maxDistance: Int,
        override val insertSequenceLength: Int) : LinkedVariant() {
    override fun toString(): String {
        return "$vcfId-$mateId"
    }
}

private data class AssemblyLinkedVariant(
        override val vcfId: String,
        override val mateId: String,
        override val minDistance: Int,
        override val maxDistance: Int,
        override val insertSequenceLength: Int) : LinkedVariant() {
    override fun toString(): String {
        return "<ASM>$vcfId-$mateId"
    }
}

private data class TransitiveLinkedVariant(
        override val vcfId: String,
        override val mateId: String,
        override val minDistance: Int,
        override val maxDistance: Int,
        override val insertSequenceLength: Int) : LinkedVariant() {
    override fun toString(): String {
        return "<TRS>$vcfId-$mateId"
    }
}
