package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

data class LinkRescue(val rescues: Set<String>) {
    companion object {

        operator fun invoke(minRescueQual: Int, links: LinkStore, softFilterStore: SoftFilterStore, variantStore: VariantStore): LinkRescue {
            val rescues = mutableSetOf<String>()

            for (vcfId in links.linkedVariants()) {
                val variant = variantStore.select(vcfId)

                if (softFilterStore.isEligibleForRescue(vcfId, variant.mateId) && variant.qual >= minRescueQual) {
                    val hasPassingLink = links.linkedVariants(vcfId).any { x -> softFilterStore.isPassing(x.otherVcfId) }
                    if (hasPassingLink) {
                        rescues.add(vcfId)
                        variant.mateId?.let { rescues.add(it) }
                    }
                }
            }

            return LinkRescue(rescues)
        }
    }
}