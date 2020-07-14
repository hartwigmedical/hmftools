package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.common.gripss.GripssFilters.DEDUP
import com.hartwig.hmftools.common.gripss.GripssFilters.MIN_QUAL
import com.hartwig.hmftools.gripss.GripssFilterConfig
import com.hartwig.hmftools.gripss.MATE
import com.hartwig.hmftools.gripss.PON
import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*

class SoftFilterStore(private val filters: Map<String, Set<String>>) {

    companion object {
        private val mateFiltered = setOf(MATE)

        operator fun invoke(config: GripssFilterConfig, variants: List<StructuralVariantContext>, ponFiltered: Set<String>, hotspots: Set<String>): SoftFilterStore {
            val filters = mutableMapOf<String, Set<String>>()
            for (variant in variants) {
                if (!hotspots.contains(variant.vcfId)) {
                    val variantFilters = mutableSetOf<String>()
                    if (ponFiltered.contains(variant.vcfId)) {
                        variantFilters.add(PON)
                    }
                    variantFilters.addAll(variant.softFilters(config))
                    if (variantFilters.isNotEmpty()) {
                        filters[variant.vcfId] = variantFilters
                    }
                }
            }

            return SoftFilterStore(filters)
        }
    }

    operator fun get(vcfId: String): Set<String> {
        return filters.getOrDefault(vcfId, Collections.emptySet())
    }

    fun filters(vcfId: String, mateId: String?): Set<String> {
        val result = get(vcfId)
        val isMateFiltered = mateId?.let { isPassing(it) } == false

        if (result.isEmpty() && isMateFiltered) {
            return mateFiltered
        }
        return result
    }

    fun duplicates(): Set<String> {
        return filters.entries.filter { (_, filters) -> filters.contains(DEDUP) }.map { x -> x.key }.toSet()
    }

    fun isPassing(vcfId: String) = get(vcfId).isEmpty()

    fun isFiltered(vcfId: String) = get(vcfId).isNotEmpty()

    fun containsDuplicateFilter(vcfId: String): Boolean = filters[vcfId]?.contains(DEDUP) == true

    fun isExclusivelyMinQualFiltered(vcfId: String): Boolean {
        return filters[vcfId]?.let { it.size == 1 && it.contains(MIN_QUAL) } == true
    }

    fun update(duplicates: Set<String>, rescues: Set<String>): SoftFilterStore {
        val result = mutableMapOf<String, MutableSet<String>>()
        for ((vcfId, softFilters) in filters) {
            if (!rescues.contains(vcfId)) {
                result[vcfId] = softFilters.toMutableSet()
            }
        }

        for (vcfId in duplicates) {
            result.computeIfAbsent(vcfId) { mutableSetOf() }.add(DEDUP)
        }

        return SoftFilterStore(result)
    }
}