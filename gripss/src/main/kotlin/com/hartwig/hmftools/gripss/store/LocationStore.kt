package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.bedpe.BreakendLocation
import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.util.*

class LocationStore(private val singlesMap: Map<String, List<BreakendLocation>>, private val pairedMap: Map<String, List<Pair<BreakendLocation, BreakendLocation>>>) {

    companion object {
        private const val MAX_DISTANCE = 4

        operator fun invoke(single: List<BreakendLocation>, paired: List<Pair<BreakendLocation, BreakendLocation>>): LocationStore {

            val singlesMap = mutableMapOf<String, MutableList<BreakendLocation>>()
            for (bedEntry in single) {
                singlesMap.computeIfAbsent(bedEntry.contig) { mutableListOf() }.add(bedEntry.expand(MAX_DISTANCE))
            }

            val pairedMap = mutableMapOf<String, MutableList<Pair<BreakendLocation, BreakendLocation>>>()
            for (bedpePair in paired) {
                val expandedFirst = bedpePair.first.expand(MAX_DISTANCE)
                val expandedSecond = bedpePair.second.expand(MAX_DISTANCE)

                pairedMap.computeIfAbsent(bedpePair.first.contig) { mutableListOf() }.add(Pair(expandedFirst, expandedSecond))
                pairedMap.computeIfAbsent(bedpePair.second.contig) { mutableListOf() }.add(Pair(expandedSecond, expandedFirst))
            }

            return LocationStore(singlesMap, pairedMap)
        }
    }

    fun isFound(variant: StructuralVariantContext): Boolean {
        return if (variant.isSingle) {
            isSingleFound(variant)
        } else {
            isPairFiltered(variant)
        }
    }

    private fun isSingleFound(sgl: StructuralVariantContext): Boolean {
        return singlesMap.getOrDefault(sgl.vcfId, Collections.emptyList()).any { x -> sgl.startLocation.isEquivalent(x) }
    }

    private fun isPairFiltered(variant: StructuralVariantContext): Boolean {
        val startLocation = variant.startLocation
        val endLocation = variant.endLocation!!
        return pairedMap.getOrDefault(variant.vcfId, Collections.emptyList()).any { (otherFirst, otherSecond) -> otherFirst.isEquivalent(startLocation) && otherSecond.isEquivalent(endLocation) }
    }
}