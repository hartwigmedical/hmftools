package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.DEDUP
import com.hartwig.hmftools.gripss.GripssFilterConfig
import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*

class SoftFilterStore(private val filters: Map<String, List<String>>) {

    companion object {
        operator fun invoke(config: GripssFilterConfig, variants: List<StructuralVariantContext>): SoftFilterStore {
            val filters = mutableMapOf<String, List<String>>()
            for (variant in variants) {
                filters[variant.vcfId] = variant.softFilters(config)
            }

            return SoftFilterStore(filters)
        }
    }

    operator fun get(vcfId: String): List<String> {
        return filters.getOrDefault(vcfId, Collections.emptyList())
    }

    fun update(duplicates: Set<String>, rescues: Set<String>): SoftFilterStore {
        val result = mutableMapOf<String, MutableList<String>>()
        for ((vcfId, softFilters) in filters) {
            if (!rescues.contains(vcfId)) {
                result[vcfId] = softFilters.toMutableList()
            }
        }

        for (vcfId in duplicates) {
            result.computeIfAbsent(vcfId) {mutableListOf()}.add(DEDUP)
        }


        return SoftFilterStore(result)
    }
}