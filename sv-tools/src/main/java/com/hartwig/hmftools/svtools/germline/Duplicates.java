package com.hartwig.hmftools.svtools.germline;

import java.util.List;

import com.google.common.collect.Lists;

public final class Duplicates
{
    private static final int MAX_DEDUP_SGL_SEEK_DISTANCE = 1000;
    private static final int MAX_DEDUP_SGL_ADDITIONAL_DISTANCE = 0;

    public static List<String> dedupSingles()
    {

        return Lists.newArrayList();
    }

    /*
    class DedupSingle(val duplicates: Set<String>) {

    companion object {

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
     */
}
