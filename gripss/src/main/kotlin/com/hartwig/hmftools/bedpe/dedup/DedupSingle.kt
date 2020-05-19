package com.hartwig.hmftools.bedpe.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.VariantStore

class DedupSingle {

    companion object {
        fun create(variantStore: VariantStore, variants: List<StructuralVariantContext>) {

            for (sgl in variants.filter { x -> x.isSingle }) {

                val others = variantStore.selectOthersInConfidenceIntervals(sgl) {x -> x.orientation == sgl.orientation}
                if (others.isNotEmpty()) {
                    println("Single: $sgl");
                    for (other in others) {
                        println("Nearby: $other")
                    }
                }
            }
        }
    }

}