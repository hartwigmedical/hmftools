package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

data class LinkRescue(val rescues: Set<String>) {
    companion object {

        operator fun invoke(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore, rescueShort: Boolean): LinkRescue {
            val rescues = mutableSetOf<String>()
            val isShort = { x: StructuralVariantContext -> x.isShortDel || x.isShortIns || x.isShortDup }
            val isEligibleForRescue = { x: StructuralVariantContext -> softFilterStore.isEligibleForRescue(x.vcfId, x.mateId) && (rescueShort || !isShort(x)) }

            for (vcfId in links.linkedVariants()) {
                if (rescues.contains(vcfId)) {
                    continue
                }

                val variant = variantStore.select(vcfId)
                if (isEligibleForRescue(variant)) {
                    val allLinkedVariants = allLinkedVariants(variant.vcfId, links, variantStore)
                    val anyPassing = allLinkedVariants.any { softFilterStore.isPassing(it) }
                    if (anyPassing) {
                        for (linkedVariant in allLinkedVariants.map { variantStore.select(it) }) {
                            if (isEligibleForRescue(linkedVariant)) {
                                rescues.add(linkedVariant.vcfId)
                                linkedVariant.mateId?.let { rescues.add(it) }
                            }
                        }
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
    }
}