package com.hartwig.hmftools.paddle.dnds

import java.io.File
import java.nio.file.Files

data class DndsMutation(
        val sample: String,
        val contig: String, val position: Int, val ref: String, val alt: String,
        val worstCodingEffect: String, val canonicalCodingEffect: String,
        val repeatCount: Int, val biallelic: Boolean, val hotspot: Boolean,
        val gene: String, val dndsImpact: String) : Comparable<DndsMutation> {

    companion object {
        fun fromFile(file: String): List<DndsMutation> {
            return Files.readAllLines(File(file).toPath()).map { fromString(it) }
        }

        fun fromString(line: String): DndsMutation {
            val array = line.split("\t")
            return DndsMutation(
                    array[0],
                    array[1],
                    array[2].toInt(),
                    array[3],
                    array[4],
                    array[5],
                    array[6],
                    array[7].toInt(),
                    array[8] == "1",
                    array[9] == "HOTSPOT",
                    array[10],
                    array[11])
        }
    }

    val impact = impact(dndsImpact, canonicalCodingEffect, worstCodingEffect)

    fun impact(dndsImpact: String, canonicalCodingEffect: String, worstCodingEffect: String): Impact {
        if (dndsImpact == "no-SNV" || dndsImpact == "NA") {
            when (canonicalCodingEffect) {
                "MISSENSE" -> return Impact.INFRAME
                "NONSENSE_OR_FRAMESHIFT" -> return Impact.FRAMESHIFT
                "SYNONYMOUS" -> return Impact.SYNONYMOUS
            }

            when (worstCodingEffect) {
                "MISSENSE" -> return Impact.INFRAME
                "NONSENSE_OR_FRAMESHIFT" -> return Impact.FRAMESHIFT
                "SYNONYMOUS" -> return Impact.SYNONYMOUS
            }

            return Impact.UNKNOWN
        }

        return when(dndsImpact) {
            "Missense" -> Impact.MISSENSE
            "Nonsense" -> Impact.NONSENSE
            "Stop_loss" -> Impact.NONSENSE
            "Synonymous" -> Impact.SYNONYMOUS
            "Essential_Splice" -> Impact.SPLICE
            else -> throw IllegalArgumentException("Unknown impact: $dndsImpact $canonicalCodingEffect $worstCodingEffect")
        }
    }

    override fun compareTo(other: DndsMutation): Int {
        fun Boolean.toInt() = if (this) 1 else 0

        if (this.hotspot.xor(other.hotspot)) {
            return other.hotspot.toInt() - this.hotspot.toInt()
        }

        if (this.biallelic.xor(other.biallelic)) {
            return other.biallelic.toInt() - this.biallelic.toInt()
        }

        return this.impact.ordinal - other.impact.ordinal
    }
}
