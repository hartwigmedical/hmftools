package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.*
import org.apache.logging.log4j.LogManager
import java.util.*

class SoftFilterStore(private val filters: Map<String, Set<String>>) {

    companion object {
        private val logger = LogManager.getLogger(this::class.java)
        private val mateFiltered = setOf(MATE)

        operator fun invoke(config: GripssFilterConfig, variants: List<StructuralVariantContext>, ponStore: LocationStore, hotspotStore: LocationStore): SoftFilterStore {
            val filters = mutableMapOf<String, Set<String>>()
            for (variant in variants) {
                if (!hotspotStore.contains(variant)) {
                    val variantFilters = mutableSetOf<String>()
                    if (ponStore.contains(variant)) {
                        variantFilters.add(PON)
                    }
                    variantFilters.addAll(variant.softFilters(config))
                    if (variantFilters.isNotEmpty()) {
                        filters[variant.vcfId] = variantFilters
                    }
                } else {
                    logger.debug("Found hotspot: $variant")
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

    fun isEligibleForRescue(vcfId: String, mateId: String?): Boolean {
        return isEligibleForRescue(vcfId) && mateId?.let { isDuplicate(it) } != true
    }

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