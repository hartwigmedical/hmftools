package com.hartwig.hmftools.gripsskt.link

import com.hartwig.hmftools.gripsskt.GripssFilterConfig
import com.hartwig.hmftools.gripsskt.StructuralVariantContext
import com.hartwig.hmftools.gripsskt.store.LinkStore
import com.hartwig.hmftools.gripsskt.store.SoftFilterStore
import com.hartwig.hmftools.gripsskt.store.VariantStore

data class LinkRescue(val rescues: Set<String>) {
    companion object {

        fun rescueDsb(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
            return rescue(links, softFilterStore, variantStore, false)
        }

        fun rescueAssembly(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
            return rescue(links, softFilterStore, variantStore, true)
        }

        fun rescueTransitive(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
            return rescue(links, softFilterStore, variantStore, true)
        }

        private fun rescue(links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore, rescueShort: Boolean): LinkRescue {
            val rescues = mutableSetOf<String>()
            val isValid = { x: StructuralVariantContext -> (rescueShort || !x.isTooShortToRescue) && !softFilterStore.containsDuplicateFilter(x.vcfId, x.mateId) }
            val isRescueCandidate = { x: StructuralVariantContext -> isValid(x) && softFilterStore.isFiltered(x.vcfId) }

            for (vcfId in links.linkedVariants()) {
                if (rescues.contains(vcfId)) {
                    continue
                }

                val variant = variantStore.select(vcfId)
                if (isRescueCandidate(variant)) {
                    val allLinkedVariants = allLinkedVariants(variant.vcfId, links, variantStore).map { variantStore.select(it) }
                    val anyPassing = allLinkedVariants.any { x -> isValid(x) && softFilterStore.isPassing(x.vcfId) }
                    if (anyPassing) {
                        for (linkedVariant in allLinkedVariants ){
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

        fun rescueDsbMobileElementInsertion(config: GripssFilterConfig, links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
            val rescues = mutableSetOf<String>()
            val isValid = { x: StructuralVariantContext -> !softFilterStore.containsPONFilter(x.vcfId, x.mateId) && !softFilterStore.containsDuplicateFilter(x.vcfId, x.mateId) && !x.isTooShortToRescue }
            val isRescueCandidate = { x: StructuralVariantContext -> softFilterStore.isFiltered(x.vcfId) && isValid(x) }

            for (vcfId in links.linkedVariants()) {
                if (rescues.contains(vcfId)) {
                    continue
                }

                for (link in links.linkedVariants(vcfId)) {
                    val variant = variantStore.select(vcfId)
                    if (isRescueCandidate(variant)) {
                        val other = variantStore.select(link.otherVcfId)
                        val combinedQual = variant.tumorQual + other.tumorQual
                        if (combinedQual >= config.minQualRescueMobileElementInsertion && isValid(other) && (variant.isMobileElementInsertion || other.isMobileElementInsertion)) {
                            rescues.add(variant.vcfId)
                            variant.mateId?.let { rescues.add(it) }
                            rescues.add(other.vcfId)
                            other.mateId?.let { rescues.add(it) }
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