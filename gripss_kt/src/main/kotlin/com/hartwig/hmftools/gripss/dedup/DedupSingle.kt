package com.hartwig.hmftools.gripss.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore

class DedupSingle(val duplicates: Set<String>) {

    companion object {
        private const val MAX_DEDUP_SGL_SEEK_DISTANCE = 1000
        private const val MAX_DEDUP_SGL_ADDITIONAL_DISTANCE = 0

        operator fun invoke(variantStore: VariantStore, softFilterStore: SoftFilterStore, linkStore: LinkStore): DedupSingle {
            val duplicates = mutableSetOf<String>()

            for (sgl in variantStore.selectAll().filter { x -> x.isSingle }) {
                val sglPasses = softFilterStore.isPassing(sgl.vcfId)

                val exactPositionFilter = { other: StructuralVariantContext -> other.start >= sgl.minStart && other.start <= sgl.maxStart }
                val duplicateFilter = { other: StructuralVariantContext -> other.orientation == sgl.orientation && (other.precise || exactPositionFilter(other)) }
                val others = variantStore.selectOthersNearby(sgl, MAX_DEDUP_SGL_ADDITIONAL_DISTANCE, MAX_DEDUP_SGL_SEEK_DISTANCE, duplicateFilter)
                if (!others.all { x -> keepSingle(sglPasses, sgl, x, softFilterStore, linkStore) }) {
                    duplicates.add(sgl.vcfId)
                } else {
                    others.forEach { x ->
                        x.vcfId.let { duplicates.add(it) }
                        x.mateId?.let { duplicates.add(it) }
                    }
                }
            }
            return DedupSingle(duplicates)
        }

        private fun keepSingle(originalPass: Boolean, original: StructuralVariantContext, alternative: StructuralVariantContext, softFilterStore: SoftFilterStore, linkStore: LinkStore): Boolean {
            if (linkStore.linkedVariants(alternative.vcfId).isNotEmpty()) {
                return false
            }

            val alternativePass = softFilterStore.isPassing(alternative.vcfId)
            if (originalPass != alternativePass) {
                return originalPass
            }

            return original.tumorQual > alternative.tumorQual
        }
    }
}