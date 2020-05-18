package com.hartwig.hmftools.extensions.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import java.util.*

class AssemblyDedup(private val linkStore: LinkStore, private val variantStore: VariantStore) {

    fun dedup(variant: StructuralVariantContext): Boolean {

        if (!variant.isSingle && variant.imprecise) {
            val target = variantStore.select(variant.mateId!!)

            val alternativeStart = variantStore.selectAlternatives(variant)
            for (alternative in alternativeStart) {
                val assemblyLinks = assemblyLinks(alternative, target, mutableListOf())
                if (assemblyLinks.isNotEmpty()) {
                    println("$variant CIPOS:${variant.confidenceInterval} IMPRECISE:${variant.imprecise} -> ${assemblyLinks.joinToString(",")}")
                    return true
                }
            }

        }

        return false
    }

    private fun assemblyLinks(current: StructuralVariantContext, target: StructuralVariantContext, path: List<String>): List<String> {
        if (current.isSingle) {
            return Collections.emptyList()
        }

        val mate = variantStore.select(current.mateId!!)
        val newPath = path + current.vcfId + mate.vcfId
        if (matchTarget(mate, target)) {
            return newPath
        }

        val linkedVariants = linkStore.followLinks(current.mateId).map { vcfId -> variantStore.select(vcfId) }
        for (linkedVariant in linkedVariants) {
            val newAssemblyLinks = assemblyLinks(linkedVariant, target, newPath)
            if (newAssemblyLinks.isNotEmpty()) {
                return newAssemblyLinks
            }
        }

        return Collections.emptyList()
    }

    private fun matchTarget(current: StructuralVariantContext, target: StructuralVariantContext): Boolean {
        return current.vcfId != target.vcfId && current.mateId!! != target.vcfId && current.orientation == target.orientation && target.confidenceIntervalsOverlap(current)
    }

}