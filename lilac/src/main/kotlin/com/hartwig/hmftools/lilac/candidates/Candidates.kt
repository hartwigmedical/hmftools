package com.hartwig.hmftools.lilac.candidates

import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence
import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import org.apache.logging.log4j.LogManager

class Candidates(private val config: LilacConfig,
                 private val nucleotideSequences: List<HlaSequenceLoci>,
                 private val aminoAcidSequences: List<HlaSequenceLoci>) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    fun candidates(context: HlaContext, fragments: List<AminoAcidFragment>, phasedEvidence: List<PhasedEvidence>): List<HlaSequenceLoci> {
        val gene = context.gene
        val aminoAcidBoundary = context.aminoAcidBoundaries

        logger.info("Determining initial candidate set for gene HLA-$gene")
        val aminoAcidCounts = SequenceCount.aminoAcids(config.minEvidence, fragments)
        val nucleotideCounts = SequenceCount.nucleotides(config.minEvidence, fragments)
        aminoAcidCounts.writeVertically("${config.outputFilePrefix}.aminoacids.${gene}.count.txt")
        nucleotideCounts.writeVertically("${config.outputFilePrefix}.nucleotides.${gene}.count.txt")

        val geneCandidates = aminoAcidSequences.filter { it.allele.gene == gene }
        logger.info(" ... ${geneCandidates.size} candidates before filtering")

        // Amino acid filtering
        val aminoAcidFilter = AminoAcidFiltering(aminoAcidBoundary)
        val aminoAcidCandidates = aminoAcidFilter.aminoAcidCandidates(geneCandidates, aminoAcidCounts)
        val aminoAcidCandidateAlleles = aminoAcidCandidates.map { it.allele }.toSet()
        val aminoAcidSpecificAllelesCandidate = aminoAcidCandidateAlleles.map { it.asFourDigit() }.toSet()
        logger.info(" ... ${aminoAcidCandidates.size} candidates after amino acid filtering")

        // Nucleotide filtering
        val nucleotideFiltering = NucleotideFiltering(config.minEvidence, aminoAcidBoundary)
        val nucleotideCandidatesAfterAminoAcidFiltering = nucleotideSequences
                .filter { it.allele.asFourDigit() in aminoAcidSpecificAllelesCandidate }
        val nucleotideSpecificAllelesCandidate = nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(nucleotideCandidatesAfterAminoAcidFiltering, fragments)
                .map { it.allele.asFourDigit() }
                .toSet()

        val nucleotideCandidates = aminoAcidCandidates.filter { it.allele.asFourDigit() in nucleotideSpecificAllelesCandidate }
        logger.info(" ... ${nucleotideCandidates.size} candidates after exon boundary filtering")

        val phasedCandidates = filterCandidates(nucleotideCandidates, phasedEvidence)
        logger.info(" ... ${phasedCandidates.size} candidates after phasing: " + phasedCandidates.map { it.allele }.joinToString(", "))

        return phasedCandidates

    }

    private fun filterCandidates(initialCandidates: List<HlaSequenceLoci>, evidence: List<PhasedEvidence>): List<HlaSequenceLoci> {
        var candidates = initialCandidates
        for (i in evidence.indices) {
            val newEvidence = evidence[i]
            candidates = candidates.filter { it.consistentWithAny(newEvidence.evidence.keys, *newEvidence.aminoAcidIndices) }
        }

        return candidates
    }
}