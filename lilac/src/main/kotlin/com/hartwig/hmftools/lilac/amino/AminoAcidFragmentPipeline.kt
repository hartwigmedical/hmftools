package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.SequenceCountDiff
import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideQualEnrichment
import com.hartwig.hmftools.lilac.nuc.NucleotideSpliceEnrichment
import org.apache.logging.log4j.LogManager
import java.util.function.Predicate

class AminoAcidFragmentPipeline(private val config: LilacConfig, private val referenceFragments: List<NucleotideFragment>, tumorFragments: List<NucleotideFragment>) {
    private val minBaseQuality = config.minBaseQual
    private val minEvidence = config.minEvidence
    private val aminoAcidEnricher = AminoAcidQualEnrichment(minEvidence)
    private val nucleotideQualEnrichment = NucleotideQualEnrichment(minBaseQuality, minEvidence)
    private val highQualityTumorFragments = qualFiltered(minBaseQuality, tumorFragments)



    companion object {
        const val MAX_AMINO_ACID_BOUNDARY = 298

        val logger = LogManager.getLogger(this::class.java)
    }


    fun referencePhasingFragments(context: HlaContext): List<AminoAcidFragment> {
        val gene = "HLA-${context.gene}"
        val geneReferenceFragments = referenceFragments.filter { it.genes.contains(gene) }

        val referenceAminoAcids = process(context.aminoAcidBoundaries, geneReferenceFragments)
        val referenceNucleotideCounts = SequenceCount.nucleotides(minEvidence, referenceAminoAcids)
        val referenceAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, referenceAminoAcids)
        referenceAminoAcidCounts.writeVertically("${config.outputFilePrefix}.${gene}.aminoacids.txt")
        referenceNucleotideCounts.writeVertically("${config.outputFilePrefix}.${gene}.nucleotides.txt")

        return referenceAminoAcids
    }

    fun referenceCoverageFragments(): List<AminoAcidFragment> {
        return qualFiltered(minBaseQuality, referenceFragments)
    }

    fun tumorCoverageFragments(): List<AminoAcidFragment> {
        if (highQualityTumorFragments.isEmpty()) {
            return listOf()
        }

        val referenceAminoAcids = referenceCoverageFragments()
        val referenceNucleotideCounts = SequenceCount.nucleotides(minEvidence, referenceAminoAcids)
        val referenceAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, referenceAminoAcids)
        val tumorNucleotideCounts = SequenceCount.nucleotides(minEvidence, highQualityTumorFragments)
        val tumorAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, highQualityTumorFragments)
        val nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts).filter { it.tumorCount > 0 }
        val aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts).filter { it.tumorCount > 0 }

        val variantFilteredTumorAminoAcids = highQualityTumorFragments.filter { !it.containsVariant(nucleotideDifferences, aminoAcidDifferences) }
        return variantFilteredTumorAminoAcids
    }


    private fun AminoAcidFragment.containsVariant(nucelotideVariants: List<SequenceCountDiff>, aminoAcidVariants: List<SequenceCountDiff>): Boolean {
        return nucelotideVariants.any { this.containsNucleotideVariant(it) } || aminoAcidVariants.any { this.containsAminoAcidVariant(it) }
    }

    private fun AminoAcidFragment.containsNucleotideVariant(variant: SequenceCountDiff): Boolean {
        return this.containsNucleotide(variant.loci) && this.nucleotide(variant.loci) == variant.sequence
    }

    private fun AminoAcidFragment.containsAminoAcidVariant(variant: SequenceCountDiff): Boolean {
        return this.containsAminoAcid(variant.loci) && this.aminoAcid(variant.loci) == variant.sequence
    }


    private fun qualFiltered(minBaseQuality: Int, fragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val qualityFilteredFragments = fragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        return qualityFilteredFragments.map { it.toAminoAcidFragment() }
    }

    private fun process(boundaries: Set<Int>, fragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        if (fragments.isEmpty()) {
            return listOf()
        }

        val qualEnriched = nucleotideQualEnrichment.enrich(fragments)
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minEvidence, boundaries.filter { it <= MAX_AMINO_ACID_BOUNDARY }.toSet())
        val spliceEnriched = spliceEnricher.enrich(qualEnriched)
        val result = aminoAcidEnricher.enrich(spliceEnriched)

        return result
    }
}