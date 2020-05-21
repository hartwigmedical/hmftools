package com.hartwig.hmftools.bedpe

import java.io.File
import java.nio.file.Files

interface Location {
    val contig: String
    val start: Int
    val end: Int
    val orientation: Byte

    fun isEquivalent(other: Location): Boolean {
        return contig == other.contig && orientation == other.orientation && other.start <= end && other.end >= start
    }
}

data class Breakend(override val contig: String, override val start: Int, override val end: Int, override val orientation: Byte) : Location {

    companion object {
        fun fromBedFile(file: String): List<Breakend> {
            return Files.readAllLines(File(file).toPath()).map { fromBed(it) }
        }

        fun fromBed(line: String): Breakend {
            val (contig1, start1, end1, _, _, strand1) = line.split("\t")
            return Breakend(contig1, start1.toInt() + 1, end1.toInt(), strand1.toOrientation())
        }
    }


    fun expand(distance: Int): Breakend {
        return this.copy(start = start - distance, end = end + distance)
    }

}

data class Breakpoint(val startBreakend: Breakend, val endBreakend: Breakend) : Location {
    override val contig: String = startBreakend.contig
    override val start: Int = startBreakend.start
    override val end: Int = startBreakend.end
    override val orientation = startBreakend.orientation

    companion object {
        fun fromBedpeFile(file: String): List<Breakpoint> {
            return Files.readAllLines(File(file).toPath()).map { fromBedpe(it) }
        }

        fun fromBedpe(line: String): Breakpoint {
            val (contig1, start1, end1, contig2, start2, end2, _, _, strand1, strand2) = line.split("\t")
            val start = Breakend(contig1, start1.toInt() + 1, end1.toInt(), strand1.toOrientation())
            val end = Breakend(contig2, start2.toInt() + 1, end2.toInt(), strand2.toOrientation())
            return Breakpoint(start, end);
        }
    }
}


private fun String.toOrientation(): Byte = if (this == "+") 1 else -1

private operator fun <T> List<T>.component6(): T = get(5)
private operator fun <T> List<T>.component7(): T = get(6)
private operator fun <T> List<T>.component8(): T = get(7)
private operator fun <T> List<T>.component9(): T = get(8)
private operator fun <T> List<T>.component10(): T = get(9)