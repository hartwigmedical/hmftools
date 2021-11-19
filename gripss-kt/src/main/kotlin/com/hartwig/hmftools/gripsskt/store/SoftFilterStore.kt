package com.hartwig.hmftools.gripsskt.store

import com.hartwig.hmftools.common.gripss.GripssFilters.DEDUP
import com.hartwig.hmftools.gripsskt.ContigComparator
import com.hartwig.hmftools.gripsskt.GripssFilterConfig
import com.hartwig.hmftools.gripsskt.PON
import com.hartwig.hmftools.gripsskt.StructuralVariantContext
import java.util.*

class SoftFilterStore(private val filters: Map<String, Set<String>>) {

    companion object {
        operator fun invoke(config: GripssFilterConfig, comparator: ContigComparator, variants: List<StructuralVariantContext>, ponFiltered: Set<String>, hotspots: Set<String>): SoftFilterStore {
            val filters = mutableMapOf<String, Set<String>>()
            for (variant in variants) {
                if (!hotspots.contains(variant.vcfId)) {
                    val variantFilters = mutableSetOf<String>()
                    if (ponFiltered.contains(variant.vcfId)) {
                        variantFilters.add(PON)
                    }
                    variantFilters.addAll(variant.softFilters(config, comparator))
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
        return get(vcfId) + (mateId?.let { get(it) } ?: setOf())
    }

    fun duplicates(): Set<String> {
        return filters.entries.filter { (_, filters) -> filters.contains(DEDUP) }.map { x -> x.key }.toSet()
    }

    fun isPassing(vcfId: String) = get(vcfId).isEmpty()

    fun isFiltered(vcfId: String) = get(vcfId).isNotEmpty()

    fun containsDuplicateFilter(vcfId: String, mateId: String?): Boolean {
        return containsFilter(DEDUP, vcfId, mateId)
    }

    fun containsPONFilter(vcfId: String, mateId: String?): Boolean {
        return containsFilter(PON, vcfId, mateId)
    }

    fun containsFilter(filter: String, vcfId: String, mateId: String?): Boolean {
        return (filters[vcfId]?.contains(filter) ?: false) || (mateId?.let { filters[it]?.contains(filter) } ?: false)
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