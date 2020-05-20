package com.hartwig.hmftools.bedpe.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.VariantStore

class DedupSingle(val duplicates: Set<String>) {

    companion object {
        fun create(variantStore: VariantStore, variants: List<StructuralVariantContext>): DedupSingle {
            val duplicates = mutableSetOf<String>()

            // Any single breakend which matches the position and orientation of another breakend (within CIPOS bounds) is also filtered
            // as DEDUP (if the matching variant is also a single breakend then the highest scoring passing breakend will be kept).

            for (sgl in variants.filter { x -> x.isSingle }) {

                val others = variantStore.selectOthersInConfidenceIntervals(sgl) {x -> x.orientation == sgl.orientation}
                if (others.isNotEmpty()) {
                    println("Single: $sgl");
                    for (other in others) {
                        println("Nearby: $other")
                    }
                }
            }

            return DedupSingle(duplicates)
        }
    }

}