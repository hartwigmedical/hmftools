package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.DEDUP
import com.hartwig.hmftools.gripss.GripssFilterConfig
import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*

class SoftFilterStore(private val filters: Map<String, Set<String>>) {

    companion object {
        operator fun invoke(config: GripssFilterConfig, variants: List<StructuralVariantContext>): SoftFilterStore {
            val filters = mutableMapOf<String, Set<String>>()
            for (variant in variants) {
                filters[variant.vcfId] = variant.softFilters(config)
            }

            return SoftFilterStore(filters)
        }
    }

    operator fun get(vcfId: String): Set<String> {
        return filters.getOrDefault(vcfId, Collections.emptySet())
    }

    fun duplicates(): Set<String> {
        return filters.entries.filter { (_, filters) -> filters.contains(DEDUP) }.map { x -> x.key }.toSet()
    }

    fun isPassing(vcfId: String) = get(vcfId).isEmpty()

    fun isEligibleForRescue(vcfId: String): Boolean = filters[vcfId]?.contains(DEDUP) == false

    fun isDuplicate(vcfId: String): Boolean = filters[vcfId]?.contains(DEDUP) == true

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