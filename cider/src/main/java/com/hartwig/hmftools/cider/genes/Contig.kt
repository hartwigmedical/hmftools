package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix

enum class ContigType {
    HUMAN_CHROMOSOME,   // E.g. 1, 2, ..., X, Y
    MITOCHONDRIAL,      // E.g. M
    ALTERNATE_LOCUS,    // E.g. xxx_alt
    UNLOCALISED,        // E.g. xxx_random
    UNPLACED,           // E.g. un_xxx
    OTHER;

    fun inPrimaryAssembly() = this == HUMAN_CHROMOSOME || this == MITOCHONDRIAL || this == UNLOCALISED || this == UNPLACED
}

data class Contig(
    val name: String,
    val type: ContigType,
    // Only present if it's a main human chromosome contig. I.e. 1, 2, ..., X, Y
    val chromosome: HumanChromosome?
) {
    override fun toString() = name

    companion object {
        fun fromName(contigName: String): Contig {
            if (contigName.isEmpty()) {
                throw IllegalArgumentException("contigName must not be empty")
            }

            // Example possibilities
            // chr1
            // chrM
            // chr1_KI270763v1_alt
            // chr15_KI270727v1_random
            // chrUn_KI270423v1

            val normalisedContig = stripChrPrefix(contigName)
            val humanChromosome = try {
                HumanChromosome.fromString(normalisedContig)
            } catch (e: IllegalArgumentException) {
                null
            }
            val parts = normalisedContig.split("_", limit = 1)
            val baseChromosome = parts[0]
            val extension = parts.getOrNull(1)
            val isAltLocus = extension?.endsWith("_alt") ?: false
            val isUnlocalised = extension?.endsWith("_random") ?: false
            val baseChromosomeLower = baseChromosome.lowercase()
            val isUnplaced = baseChromosomeLower.startsWith("un_")
            val isMitochondrial = baseChromosomeLower == "mt" || baseChromosomeLower == "m"

            val contigType = if (isUnplaced) {
                ContigType.UNPLACED
            } else if (isAltLocus) {
                ContigType.ALTERNATE_LOCUS
            } else if (isUnlocalised) {
                ContigType.UNLOCALISED
            } else if (isMitochondrial) {
                ContigType.MITOCHONDRIAL
            } else if (humanChromosome != null) {
                ContigType.HUMAN_CHROMOSOME
            } else {
                ContigType.OTHER
            }

            return Contig(
                contigName,
                contigType,
                if (contigType == ContigType.HUMAN_CHROMOSOME) humanChromosome else null
            )
        }
    }
}
