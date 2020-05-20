package com.hartwig.hmftools.bedpe.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.link.AlternatePath
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

data class DedupPair(val duplicates: Set<String>, val rescue: Set<String>) {

    companion object {
        operator fun invoke(filterStore: SoftFilterStore, alternatePaths: Collection<AlternatePath>, variantStore: VariantStore): DedupPair {
            val rescues = mutableSetOf<String>()
            val duplicates = mutableSetOf<String>()

            for (alt in alternatePaths) {
                val originalPasses = filterStore[alt.vcfId].isEmpty()
                val anyInAltPathPasses = alt.pathVcfIds().map { x -> filterStore[x] }.any { x -> x.isEmpty() }

                if (alt.size() == 1) {
                    val variant = variantStore.select(alt.vcfId)
                    val other = variantStore.select(alt.path[0].otherVcfId)

                    // Favour PRECISE, PASSING, then QUAL
                    if (!keepOriginal(variant, other, originalPasses, anyInAltPathPasses)) {
                        duplicates.add(alt.vcfId)
                    }
                } else {
                    duplicates.add(alt.vcfId)
                    duplicates.add(alt.mateId)

                    if (originalPasses || anyInAltPathPasses) {
                        rescues.addAll(alt.pathVcfIds())
                    }
                }
            }

            // duplicates supersede rescues
            rescues.removeIf { x -> duplicates.contains(x) }
            return DedupPair(duplicates, rescues)
        }

        private fun keepOriginal(variant: StructuralVariantContext, other: StructuralVariantContext, variantPass: Boolean, otherPass: Boolean): Boolean {
            if (!variant.imprecise && other.imprecise) {
                return true;
            }

            if (variantPass && !otherPass) {
                return true
            }

            return variant.qual > other.qual
        }
    }

}