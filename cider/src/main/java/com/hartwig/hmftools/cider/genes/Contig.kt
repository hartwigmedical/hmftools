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

    val inPrimaryAssembly: Boolean
        get() = this == HUMAN_CHROMOSOME || this == MITOCHONDRIAL || this == UNLOCALISED || this == UNPLACED
}

data class Contig(
    // Example possibilities:
    // chr1
    // chrM
    // chr1_KI270763v1_alt
    // chr15_KI270727v1_random
    // chrUn_KI270423v1
    // The chr prefix may or may not be present and does affect the identity of the contig.
    val name: String,
) {
    val type: ContigType
        get() =
            if (chromosome != null) {
                ContigType.HUMAN_CHROMOSOME
            }
            else if (isUnplaced) {
                ContigType.UNPLACED
            } else if (isAltLocus) {
                ContigType.ALTERNATE_LOCUS
            } else if (isUnlocalised) {
                ContigType.UNLOCALISED
            } else if (isMitochondrial) {
                ContigType.MITOCHONDRIAL
            } else {
                ContigType.OTHER
            }

    val chromosome: HumanChromosome?
        get() = try {
            HumanChromosome.fromString(normalisedContig)
        } catch (e: IllegalArgumentException) {
            null
        }

    override fun toString() = name

    private val normalisedContig: String
        get() = stripChrPrefix(name)

    private val baseChromosome: String
        get() = normalisedContig.split("_", limit = 2)[0].lowercase()

    private val extension: String?
        get() = normalisedContig.split("_", limit = 2).getOrNull(1)

    private val isAltLocus: Boolean
        get() = extension?.endsWith("_alt") ?: false

    private val isUnlocalised: Boolean
        get() = extension?.endsWith("_random") ?: false

    private val isUnplaced: Boolean
        get() = baseChromosome.startsWith("un_")

    private val isMitochondrial: Boolean
        get() = baseChromosome == "mt" || baseChromosome == "m"
}
