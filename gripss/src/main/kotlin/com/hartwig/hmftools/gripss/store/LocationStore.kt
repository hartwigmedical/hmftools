package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.bedpe.Breakend
import com.hartwig.hmftools.bedpe.Breakpoint
import com.hartwig.hmftools.bedpe.Location
import com.hartwig.hmftools.gripss.ContigComparator
import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.io.Serializable

class LocationStore private constructor(private val compare: ContigComparator, private val singlesMap: Map<String, LocationSeek<Breakend>>, private val pairedMap: Map<String, LocationSeek<Breakpoint>>) : Serializable {

    companion object {
        operator fun invoke(compare: ContigComparator, single: List<Breakend>, paired: List<Breakpoint>, additionalBuffer: Int = 0): LocationStore {

            val singlesMap = mutableMapOf<String, MutableList<Breakend>>()
            for (bedEntry in single) {
                val expanded = bedEntry.expand(additionalBuffer)
                bedEntry.locationKey().forEach { key -> singlesMap.computeIfAbsent(key) { mutableListOf() }.add(expanded) }
            }

            val pairedMap = mutableMapOf<String, MutableList<Breakpoint>>()
            for (bedpePair in paired) {
                val expanded = bedpePair.expand(additionalBuffer)
                bedpePair.locationKey().forEach { key -> pairedMap.computeIfAbsent(key) { mutableListOf() }.add(expanded) }
            }

            val breakendSeeks: Map<String, LocationSeek<Breakend>> = singlesMap.entries.associate { (key: String, list: List<Breakend>) -> Pair(key, LocationSeek(list)) }
            val breakpointSeeks: Map<String, LocationSeek<Breakpoint>> = pairedMap.entries.associate { (key: String, list: List<Breakpoint>) -> Pair(key, LocationSeek(list)) }

            return LocationStore(compare, breakendSeeks, breakpointSeeks)
        }

        private fun Breakend.locationKey(): List<String> {
            val startPositionKey = roundedPosition(this.start)
            val endPositionKey = roundedPosition(this.start)

            val result = mutableListOf<String>()
            result.add("${this.contig}:${this.orientation}:${startPositionKey}")
            if (startPositionKey != endPositionKey) {
                result.add("${this.contig}:${this.orientation}:${endPositionKey}")
            }

            return result
        }

        private fun Breakpoint.locationKey(): List<String> {
            val startKeys = startBreakend.locationKey()
            val endKeys = endBreakend.locationKey()

            val result = mutableListOf<String>()
            for (startKey in startKeys) {
                for (endKey in endKeys) {
                    result.add("$startKey:$endKey")
                }
            }

            return result
        }

        private fun roundedPosition(position: Int) = position / 1_000_000

        private fun Location.overlaps(other: Location): Boolean {
            // No need to check for contig and orientation as the maps are already keyed by this
            return other.start <= end && other.end >= start
        }
    }

    fun contains(variant: StructuralVariantContext): Boolean {
        return if (variant.isSingle) {
            contains(variant.startBreakend)
        } else {
            contains(variant.breakpoint!!)
        }
    }

    fun contains(start: Breakend): Boolean {
        val keys = start.locationKey()
        return keys.any { singlesMap[it]?.any { x -> x.overlaps(start) } == true }
    }

    fun contains(breakpoint: Breakpoint): Boolean {
        if (compare.compare(breakpoint.startBreakend, breakpoint.endBreakend) > 0) {
            return contains(Breakpoint(breakpoint.endBreakend, breakpoint.startBreakend))
        }

        val keys = breakpoint.locationKey()
        return keys.any { pairedMap[it]?.any { x -> x.startBreakend.overlaps(breakpoint.startBreakend) && x.endBreakend.overlaps(breakpoint.endBreakend) } == true }
    }

    private fun ContigComparator.compare(breakend1: Breakend, breakend2: Breakend): Int {
        return compare(breakend1.contig, breakend1.start, breakend2.contig, breakend2.start)
    }
}