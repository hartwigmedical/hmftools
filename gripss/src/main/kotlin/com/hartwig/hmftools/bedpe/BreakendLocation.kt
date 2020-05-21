package com.hartwig.hmftools.bedpe

import java.io.File
import java.nio.file.Files

data class BreakendLocation(val contig: String, val start: Int, val end: Int, val orientation: Byte) {

    companion object {
        fun fromBedFile(file: String): List<BreakendLocation> {
            return Files.readAllLines(File(file).toPath()).map { fromBed(it) }
        }

        fun fromBed(line: String): BreakendLocation {
            fun String.toOrientation(): Byte = if (this == "+") 1 else -1

            val (contig1, start1, end1, _, _, strand1) = line.split("\t")
            return BreakendLocation(contig1, start1.toInt() + 1, end1.toInt(), strand1.toOrientation())
        }

        fun fromBedpeFile(file: String): List<Pair<BreakendLocation, BreakendLocation>> {
            return Files.readAllLines(File(file).toPath()).map { fromBedpe(it) }
        }

        fun fromBedpe(line: String): Pair<BreakendLocation, BreakendLocation> {
            fun String.toOrientation(): Byte = if (this == "+") 1 else -1

            val (contig1, start1, end1, contig2, start2, end2, _, _, strand1, strand2) = line.split("\t")
            val start = BreakendLocation(contig1, start1.toInt() + 1, end1.toInt(), strand1.toOrientation())
            val end = BreakendLocation(contig2, start2.toInt() + 1, end2.toInt(), strand2.toOrientation())
            return Pair(start, end);
        }

        operator fun <T> List<T>.component6(): T = get(5)
        operator fun <T> List<T>.component7(): T = get(6)
        operator fun <T> List<T>.component8(): T = get(7)
        operator fun <T> List<T>.component9(): T = get(8)
        operator fun <T> List<T>.component10(): T = get(9)
    }

    fun isEquivalent(other: BreakendLocation): Boolean {
        return contig == other.contig && orientation == other.orientation && other.start <= end && other.end >= start
    }

    fun expand(distance: Int): BreakendLocation {
        return this.copy(start = start - distance, end = end + distance)
    }


}
