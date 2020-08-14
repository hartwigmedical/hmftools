package com.hartwig.hmftools.paddle.dnds

import com.hartwig.hmftools.paddle.gene.Gene
import java.io.File
import java.nio.file.Files

data class DndsMutation(
        val sample: String,
        val contig: String,
        val position: Int,
        val ref: String,
        val alt: String,
        val worstCodingEffect: String,
        val canonicalCodingEffect: String,
        val repeatCount: Int,
        private val biallelic: Boolean,
        val isHotspot: Boolean,
        val gene: Gene,
        val dndsImpact: String) {

    val impact = impact(dndsImpact, canonicalCodingEffect, worstCodingEffect)
    val isBiallelic = biallelic && impact != Impact.MISSENSE

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

        return when (dndsImpact) {
            "Missense" -> Impact.MISSENSE
            "Nonsense" -> Impact.NONSENSE
            "Stop_loss" -> Impact.NONSENSE
            "Synonymous" -> Impact.SYNONYMOUS
            "Essential_Splice" -> Impact.SPLICE
            else -> throw IllegalArgumentException("Unknown impact: $dndsImpact $canonicalCodingEffect $worstCodingEffect")
        }
    }

    companion object {
        fun fromFile(file: String): List<DndsMutation> {
            return Files.readAllLines(File(file).toPath()).filter { x -> !x.startsWith("sample") }.map { fromString(it) }
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
}
