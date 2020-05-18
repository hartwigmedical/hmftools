package com.hartwig.hmftools.extensions.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import java.util.*

class TransitiveDedup(private val linkStore: LinkStore, private val variantStore: VariantStore) {

    fun dedup(variant: StructuralVariantContext): Boolean {

        if (!variant.isSingle) {
            val target = variantStore.select(variant.mateId!!)

            val alternativeStart = variantStore.selectAlternatives(variant)
            for (alternative in alternativeStart) {
                if (!alternative.isSingle) {
                    val assemblyLinks = assemblyLinks(LinkedVariantStart(alternative.vcfId, alternative.mateId!!), target, 2, mutableListOf())
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
            return newPath
        }

        // Always try assembled links first!
        val assemblyLinkedVariants = linkStore.followLinks(current.mateId).map { vcfId -> variantStore.select(vcfId) }
        for (linkedVariant in assemblyLinkedVariants) {
            if (!linkedVariant.isSingle && !linkedVariant.imprecise) {
                val nextJump = AssemblyLinkedVariant(linkedVariant.vcfId, linkedVariant.mateId!!)
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
                    val nextJump = TransitiveLinkedVariant(linkedVariant.vcfId, linkedVariant.mateId!!)
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

}

sealed class LinkedVariant {
    abstract val vcfId: String
    abstract val mateId: String
    override fun toString(): String {
        return "$vcfId-$mateId"
    }
}
data class LinkedVariantStart(override val vcfId: String, override val mateId: String) : LinkedVariant() {
    override fun toString(): String {
        return "$vcfId-$mateId"
    }
}
data class AssemblyLinkedVariant(override val vcfId: String, override val mateId: String) : LinkedVariant() {
    override fun toString(): String {
        return "<ASM>$vcfId-$mateId"
    }
}
data class TransitiveLinkedVariant(override val vcfId: String, override val mateId: String) : LinkedVariant() {
    override fun toString(): String {
        return "<TRS>$vcfId-$mateId"
    }
}
