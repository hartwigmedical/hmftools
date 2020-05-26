package com.hartwig.hmftools.gripss.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

class DedupSingle(val duplicates: Set<String>) {

    companion object {
        operator fun invoke(variantStore: VariantStore, softFilterStore: SoftFilterStore): DedupSingle {
            val duplicates = mutableSetOf<String>()

            for (sgl in variantStore.selectAll().filter { x -> x.isSingle }) {
                val sglPasses = softFilterStore.isPassing(sgl.vcfId)

                val exactPositionFilter = { other: StructuralVariantContext -> other.start >= sgl.minStart && other.start <= sgl.maxStart }
                val duplicateFilter = { other: StructuralVariantContext -> other.orientation == sgl.orientation && (other.precise || exactPositionFilter(other)) }
                val others = variantStore.selectOthersNearby(sgl, 0, duplicateFilter)
                if (!others.all { x -> keepOriginal(sglPasses, sgl, x, softFilterStore) }) {
                    duplicates.add(sgl.vcfId)
                }
            }
            return DedupSingle(duplicates)
        }

        private fun keepOriginal(originalPass: Boolean, original: StructuralVariantContext, alternative: StructuralVariantContext, softFilterStore: SoftFilterStore): Boolean {
            if (!alternative.isSingle) {
                return false
            }

            val alternativePass = softFilterStore.isPassing(alternative.vcfId)
            if (originalPass != alternativePass) {
                return originalPass
            }

            return original.qual > alternative.qual
        }
    }

}