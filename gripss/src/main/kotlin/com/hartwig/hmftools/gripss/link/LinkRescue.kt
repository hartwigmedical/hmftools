package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.GripssFilterConfig
import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

data class LinkRescue(val rescues: Set<String>) {
    companion object {

        operator fun invoke(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore, rescueShort: Boolean): LinkRescue {
            val rescues = mutableSetOf<String>()
            val isShort = { x: StructuralVariantContext -> x.isShortDel || x.isShortIns || x.isShortDup }
            val isRescueCandidate = { x: StructuralVariantContext -> softFilterStore.isRescueCandidate(x.vcfId, x.mateId) && (rescueShort || !isShort(x)) }

            for (vcfId in links.linkedVariants()) {
                if (rescues.contains(vcfId)) {
                    continue
                }

                val variant = variantStore.select(vcfId)
                if (isRescueCandidate(variant)) {
                    val allLinkedVariants = allLinkedVariants(variant.vcfId, links, variantStore)
                    val anyPassing = allLinkedVariants.any { softFilterStore.isPassing(it) }
                    if (anyPassing) {
                        for (linkedVariant in allLinkedVariants.map { variantStore.select(it) }) {
                            if (isRescueCandidate(linkedVariant)) {
                                rescues.add(linkedVariant.vcfId)
                                linkedVariant.mateId?.let { rescues.add(it) }
                            }
                        }
                    }
                }
            }

            return LinkRescue(rescues)
        }

        fun rescueSingles(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore, config: GripssFilterConfig): LinkRescue {
            val rescues = mutableSetOf<String>()
            val isRescueCandidate = { x: StructuralVariantContext -> x.isSingle && softFilterStore.isExclusivelyMinQualFiltered(x.vcfId)}

            for (vcfId in links.linkedVariants()) {
                if (rescues.contains(vcfId)) {
                    continue
                }
                val variant = variantStore.select(vcfId)
                if (isRescueCandidate(variant)) {
                    val allLinkedSingles = allLinkedVariants(variant.vcfId, links, variantStore)
                            .map { variantStore.select(it) }
                            .filter { it.isSingle && !softFilterStore.containsDuplicateFilter(it.vcfId) }
                    val combinedQual = variant.qual + allLinkedSingles.map { it.qual }.sum()
                    if (combinedQual > config.minQualBreakEnd) {
                        rescues.add(variant.vcfId)
                        allLinkedSingles.forEach { rescues.add(it.vcfId) }
                    }
                }
            }

            return LinkRescue(rescues)
        }

        private fun allLinkedVariants(vcfId: String, links: LinkStore, variantStore: VariantStore): Set<String> {
            val result = mutableSetOf<String>()
            allLinkedVariantsInner(vcfId, links, variantStore, result)
            return result
        }

        private fun allLinkedVariantsInner(vcfId: String, links: LinkStore, variantStore: VariantStore, result: MutableSet<String>) {
            if (result.contains(vcfId)) {
                return
            }

            result.add(vcfId)
            val linkedVariants = links.linkedVariants(vcfId).map { variantStore.select(it.otherVcfId) }
            for (linkedVariant in linkedVariants) {
                allLinkedVariantsInner(linkedVariant.vcfId, links, variantStore, result)
                linkedVariant.mateId?.let { allLinkedVariantsInner(it, links, variantStore, result) }
            }
        }

        fun SoftFilterStore.isRescueCandidate(vcfId: String, mateId: String?): Boolean {
            return isRescueCandidate(vcfId) && mateId?.let { !containsDuplicateFilter(it) } ?: true
        }

        fun SoftFilterStore.isRescueCandidate(vcfId: String): Boolean {
            return isFiltered(vcfId) && !containsDuplicateFilter(vcfId)
        }

    }
}