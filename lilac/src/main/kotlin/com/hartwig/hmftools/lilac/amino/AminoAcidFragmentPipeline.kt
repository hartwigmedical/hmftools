package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.SequenceCountDiff
import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideQualEnrichment
import com.hartwig.hmftools.lilac.nuc.NucleotideSpliceEnrichment

class AminoAcidFragmentPipeline(private val config: LilacConfig, private val referenceFragments: List<NucleotideFragment>, tumorFragments: List<NucleotideFragment>) {
    private val tumorAminoAcids = qualFiltered(tumorFragments)

    private val minBaseQuality = config.minBaseQual
    private val minEvidence = config.minEvidence
    private val aminoAcidEnricher = AminoAcidQualEnrichment(minEvidence)
    private val nucleotideQualEnrichment = NucleotideQualEnrichment(minBaseQuality, minEvidence)

    fun phasingFragments(context: HlaContext): List<AminoAcidFragment> {
        val gene = "HLA-${context.gene}"
        val geneReferenceFragments = referenceFragments.filter { it.genes.contains(gene) }
        val referenceAminoAcids = process(context.aminoAcidBoundaries, geneReferenceFragments)
        val referenceNucleotideCounts = SequenceCount.nucleotides(minEvidence, referenceAminoAcids)
        val referenceAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, referenceAminoAcids)
        referenceAminoAcidCounts.writeVertically("${config.outputFilePrefix}.reference.aminoacids.${gene}.count.txt")
        referenceNucleotideCounts.writeVertically("${config.outputFilePrefix}.reference.nucleotides.${gene}.count.txt")

        if (tumorAminoAcids.isEmpty()) {
            return referenceAminoAcids
        }

        val tumorAminoAcids = tumorAminoAcids.filter { it.genes.contains(gene) }
        val tumorNucleotideCounts = SequenceCount.nucleotides(minEvidence, tumorAminoAcids)
        val tumorAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, tumorAminoAcids)
        tumorAminoAcidCounts.writeVertically("${config.outputFilePrefix}.tumor.aminoacids.${gene}.count.txt")
        tumorNucleotideCounts.writeVertically("${config.outputFilePrefix}.tumor.nucleotides.${gene}.count.txt")

        val nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts).filter { it.tumorCount > 0 }
        val aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts).filter { it.tumorCount > 0 }

        val variantFilteredTumorAminoAcids = tumorAminoAcids.filter { !it.containsVariant(nucleotideDifferences, aminoAcidDifferences) }
        return referenceAminoAcids //TODO: CHANGE?
//        return referenceAminoAcids + variantFilteredTumorAminoAcids
    }

    private fun AminoAcidFragment.containsVariant(nucelotideVariants: List<SequenceCountDiff>, aminoAcidVariants: List<SequenceCountDiff>): Boolean {
        return nucelotideVariants.any { this.containsNucleotideVariant(it) } || aminoAcidVariants.any { this.containsAminoAcidVariant(it) }
    }

    private fun AminoAcidFragment.containsNucleotideVariant(variant: SequenceCountDiff): Boolean {
        return this.containsNucleotide(variant.loci) && this.nucleotide(variant.loci) == variant.sequence
    }

    fun AminoAcidFragment.containsAminoAcidVariant(variant: SequenceCountDiff): Boolean {
        return this.containsAminoAcid(variant.loci) && this.aminoAcid(variant.loci) == variant.sequence
    }

    fun referenceCoverageFragments(): List<AminoAcidFragment> {
        return qualFiltered(referenceFragments)
    }

    fun tumorCoverageFragments(): List<AminoAcidFragment> {
        if (tumorAminoAcids.isEmpty()) {
            return listOf()
        }

        val referenceAminoAcids = referenceCoverageFragments()
        val referenceNucleotideCounts = SequenceCount.nucleotides(minEvidence, referenceAminoAcids)
        val referenceAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, referenceAminoAcids)
        val tumorNucleotideCounts = SequenceCount.nucleotides(minEvidence, tumorAminoAcids)
        val tumorAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, tumorAminoAcids)
        val nucleotideDifferences = SequenceCountDiff.create(referenceNucleotideCounts, tumorNucleotideCounts).filter { it.tumorCount > 0 }
        val aminoAcidDifferences = SequenceCountDiff.create(referenceAminoAcidCounts, tumorAminoAcidCounts).filter { it.tumorCount > 0 }

        val variantFilteredTumorAminoAcids = tumorAminoAcids.filter { !it.containsVariant(nucleotideDifferences, aminoAcidDifferences) }
//        return referenceAminoAcids
        return variantFilteredTumorAminoAcids
    }

    private fun qualFiltered(fragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val qualityFilteredFragments = fragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        return qualityFilteredFragments.map { it.toAminoAcidFragment() }
    }



    fun process(boundaries: Set<Int>, fragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        if (fragments.isEmpty()) {
            return listOf()
        }

        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minEvidence, boundaries)

        val qualEnriched = nucleotideQualEnrichment.enrich(fragments)
        val spliceEnriched = spliceEnricher.enrich(qualEnriched)
        val result = aminoAcidEnricher.enrich(spliceEnriched)

        return result
    }
}