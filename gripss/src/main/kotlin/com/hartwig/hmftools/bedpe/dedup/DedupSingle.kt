package com.hartwig.hmftools.bedpe.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

class DedupSingle(val duplicates: Set<String>) {

    companion object {
        fun create(variantStore: VariantStore,  softFilterStore: SoftFilterStore, variants: List<StructuralVariantContext>): DedupSingle {
            val duplicates = mutableSetOf<String>()


            // Any passing single breakend which matches the position and orientation of another passing breakend (within CIPOS bounds)
            // is also filtered as DEDUP (if the matching variant is also a single breakend then the highest scoring breakend will be kept).

            for (sgl in variants.filter { x -> x.isSingle && softFilterStore[x.vcfId].isEmpty() }) {
                val duplicateFilter = { other: StructuralVariantContext -> other.orientation == sgl.orientation && softFilterStore[other.vcfId].isEmpty()}
                val others = variantStore.selectOthersNearby(sgl, 0, duplicateFilter)
                if (others.isNotEmpty()) {
                    println("Single: $sgl");
                    for (other in others) {
                        println("Nearby: $other")
                    }
                }
            }

            return DedupSingle(duplicates)
        }

        private fun VariantStore.selectDuplicates(variant: StructuralVariantContext): List<StructuralVariantContext> {
            val duplicate = { other: StructuralVariantContext -> other.orientation == variant.orientation}
            return selectOthersNearby(variant, 0, duplicate)
        }

    }

}