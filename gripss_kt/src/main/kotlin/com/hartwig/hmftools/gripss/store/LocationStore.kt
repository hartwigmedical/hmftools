package com.hartwig.hmftools.gripss.store

import com.hartwig.hmftools.gripss.ContigComparator
import com.hartwig.hmftools.gripss.StructuralVariantContext
import java.io.Serializable

open class LocationStore(private val compare: ContigComparator, private val singlesMap: Map<String, LocationSeek<Breakend>>, private val pairedMap: Map<String, LocationSeek<Breakpoint>>) : Serializable {

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
            val endPositionKey = roundedPosition(this.end)

            val result = mutableListOf<String>()
            for (key in startPositionKey..endPositionKey) {
                result.add("${this.contig}:${this.orientation}:${key}")
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

        private fun Location.overlaps(other: Location, buffer: Int): Boolean {
            // No need to check for contig and orientation as the maps are already keyed by this
            // return other.start <= end && other.end >= start
            return start <= other.end + buffer && end >= other.start - buffer
        }
    }

    fun contains(variant: StructuralVariantContext): Boolean {
        return if (variant.isSingle) {
            contains(variant.startBreakend, variant.inexactHomologyStart())
        } else {
            contains(variant.breakpoint!!, variant.inexactHomologyStart(), variant.inexactHomologyEnd())
        }
    }

    fun contains(start: Breakend): Boolean {
        return contains(start, 0);
    }

    fun contains(start: Breakend, inexactHom: Int): Boolean {
        val keys = start.locationKey()
        return keys.any { singlesMap[it]?.any { x -> x.overlaps(start, inexactHom) } == true }
    }

    fun contains(breakpoint: Breakpoint): Boolean {
        return contains(breakpoint, 0, 0);
    }

    fun contains(breakpoint: Breakpoint, inexactHomStart: Int, inexactHomEnd: Int): Boolean {
        if (compare.compare(breakpoint.startBreakend, breakpoint.endBreakend) > 0) {
            return contains(Breakpoint(breakpoint.endBreakend, breakpoint.startBreakend), inexactHomStart, inexactHomEnd)
        }

        val keys = breakpoint.locationKey()
        return keys.any { pairedMap[it]?.any { x -> x.startBreakend.overlaps(breakpoint.startBreakend, inexactHomStart)
                && x.endBreakend.overlaps(breakpoint.endBreakend, inexactHomEnd) } == true }
    }

    private fun ContigComparator.compare(breakend1: Breakend, breakend2: Breakend): Int {
        return compare(breakend1.contig, breakend1.start, breakend2.contig, breakend2.start)
    }
}