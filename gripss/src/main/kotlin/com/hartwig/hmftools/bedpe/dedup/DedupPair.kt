package com.hartwig.hmftools.bedpe.dedup

import com.hartwig.hmftools.gripss.link.AlternatePath
import com.hartwig.hmftools.gripss.store.SoftFilterStore

data class DedupPair(val duplicates: Set<String>, val rescue: Set<String>) {

    companion object {
        operator fun invoke(filterStore: SoftFilterStore, alternatePaths: Collection<AlternatePath>): DedupPair {
            val rescues = mutableSetOf<String>()
            val duplicates = mutableSetOf<String>()

            for (alt in alternatePaths) {
                val originalPasses = filterStore[alt.vcfId].isEmpty()
                val anyAltPasses = alt.pathVcfIds().map { x -> filterStore[x] }.any { x -> x.isEmpty() }

                if (alt.size() == 1) {


                } else {
                    duplicates.add(alt.vcfId)
                    duplicates.add(alt.mateId)

                    if (originalPasses || anyAltPasses) {
                        rescues.addAll(alt.pathVcfIds())
                    }
                }
            }

            // duplicates supersede rescues
            rescues.removeIf { x -> duplicates.contains(x) }
            return DedupPair(duplicates, rescues)
        }


    }


}