package com.hartwig.hmftools.paddle.dnds

import com.hartwig.hmftools.paddle.Gene
import com.hartwig.hmftools.paddle.Impact
import java.io.File
import java.nio.file.Files

private val EXCLUDE_HOTSPOT = setOf(Impact.SYNONYMOUS, Impact.UNKNOWN)
private val EXCLUDE_BIALLELIC = setOf(Impact.MISSENSE, Impact.SYNONYMOUS, Impact.UNKNOWN)

data class DndsMutation(
        val sample: String,
        val contig: String,
        val position: Int,
        val ref: String,
        val alt: String,
        val worstCodingEffect: String,
        val canonicalCodingEffect: String,
        val repeatCount: Int,
        val biallelic: Boolean,
        private val hotspot: Boolean,
        val gene: Gene,
        val dndsImpact: String) {

    val impact = impact(dndsImpact, canonicalCodingEffect, worstCodingEffect)
    val isBiallelic = biallelic && impact !in EXCLUDE_BIALLELIC
    val isHotspot = hotspot && impact !in EXCLUDE_HOTSPOT
    val isKnownOncoDriver = isHotspot || impact == Impact.INFRAME && repeatCount < 8
    val isKnownTsgDriver = isHotspot || isBiallelic

    private fun impact(dndsImpact: String, canonicalCodingEffect: String, worstCodingEffect: String): Impact {
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
                    array[8].toBoolean(),
                    array[9].toBoolean(),
                    array[10],
                    array[11])
        }
    }
}
