package com.hartwig.hmftools.bedpe.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

class DedupSingle(val duplicates: Set<String>) {

    companion object {
        operator fun invoke(variantStore: VariantStore, softFilterStore: SoftFilterStore): DedupSingle {
            val duplicates = mutableSetOf<String>()

            for (sgl in variantStore.selectAll().filter { x -> x.isSingle && softFilterStore.isPassing(x.vcfId) }) {
                val duplicateFilter = { other: StructuralVariantContext -> other.orientation == sgl.orientation && softFilterStore[other.vcfId].isEmpty() }
                val others = variantStore.selectOthersNearby(sgl, 0, duplicateFilter)
                if (isDuplicate(sgl, others)) {
                    duplicates.add(sgl.vcfId)
                }
            }
            return DedupSingle(duplicates)
        }

        private fun isDuplicate(sgl: StructuralVariantContext, others: List<StructuralVariantContext>): Boolean {
            return when (others.size) {
                0 -> false
                1 -> !others[0].isSingle || others[0].qual > sgl.qual
                else -> true
            }
        }
    }

}