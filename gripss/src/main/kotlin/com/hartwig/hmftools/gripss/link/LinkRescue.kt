package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore

data class LinkRescue(val rescues: Set<String>) {
    companion object {

        operator fun invoke(links: LinkStore, softFilterStore: SoftFilterStore): LinkRescue {
            val rescues = mutableSetOf<String>()

            for (vcfId in links.linkedVariants()) {
                if (softFilterStore.isEligibleForRescue(vcfId)) {
                    val hasPassingLink = links.linkedVariants(vcfId).any { x -> softFilterStore[x.otherVcfId].isEmpty() }
                    if (hasPassingLink) {
                        rescues.add(vcfId)
                    }
                }
            }

            return LinkRescue(rescues)
        }
    }
}