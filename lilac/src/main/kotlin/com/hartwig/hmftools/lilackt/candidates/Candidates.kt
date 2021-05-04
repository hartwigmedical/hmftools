package com.hartwig.hmftools.lilackt.candidates

import com.hartwig.hmftools.lilackt.LilacConfig
import com.hartwig.hmftools.lilackt.SequenceCount
import com.hartwig.hmftools.lilackt.amino.AminoAcidFragment
import com.hartwig.hmftools.lilackt.evidence.PhasedEvidence
import com.hartwig.hmftools.lilackt.hla.HlaAllele
import com.hartwig.hmftools.lilackt.hla.HlaContext
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci
import org.apache.logging.log4j.LogManager

class Candidates(private val config: LilacConfig,
                 private val nucleotideSequences: List<HlaSequenceLoci>,
                 private val aminoAcidSequences: List<HlaSequenceLoci>) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    fun unphasedCandidates(context: HlaContext, fragments: List<AminoAcidFragment>): List<HlaAllele> {
        val gene = context.gene
        val aminoAcidBoundary = context.aminoAcidBoundaries

        logger.info("Determining un-phased candidate set for gene HLA-$gene")
        val aminoAcidCounts = SequenceCount.aminoAcids(config.minEvidence, fragments)
        val geneCandidates = aminoAcidSequences.filter { it.allele.gene == gene }
        logger.info("    ${geneCandidates.size} candidates before filtering")

        // Amino acid filtering
        val aminoAcidFilter = AminoAcidFiltering(aminoAcidBoundary)
        val aminoAcidCandidates = aminoAcidFilter.aminoAcidCandidates(geneCandidates, aminoAcidCounts)
        val aminoAcidCandidateAlleles = aminoAcidCandidates.map { it.allele }.toSet()
        val aminoAcidSpecificAllelesCandidate = aminoAcidCandidateAlleles.map { it.asFourDigit() }.toSet()
        if (aminoAcidSpecificAllelesCandidate.isEmpty()) {
            logger.warn("    0 candidates after amino acid filtering - reverting to all gene candidates")
            return geneCandidates.map { it.allele }
        }

        logger.info("    ${aminoAcidCandidates.size} candidates after amino acid filtering")

        // Nucleotide filtering
        val nucleotideFiltering = NucleotideFiltering(config.minEvidence, aminoAcidBoundary)
        val nucleotideCandidatesAfterAminoAcidFiltering = nucleotideSequences
                .filter { it.allele.asFourDigit() in aminoAcidSpecificAllelesCandidate }
        val nucleotideSpecificAllelesCandidate = nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(nucleotideCandidatesAfterAminoAcidFiltering, fragments)
                .map { it.allele.asFourDigit() }
                .distinct()

        if (nucleotideSpecificAllelesCandidate.isEmpty()) {
            logger.warn("    0 candidates after exon boundary filtering - reverting to amino acid candidates")
            return aminoAcidCandidateAlleles.toList()
        }

        logger.info("    ${nucleotideSpecificAllelesCandidate.size} candidates after exon boundary filtering")
        return nucleotideSpecificAllelesCandidate
    }

    fun phasedCandidates(context: HlaContext, unphasedCandidateAlleles: List<HlaAllele>, phasedEvidence: List<PhasedEvidence>): List<HlaAllele> {
        val gene = context.gene
        logger.info("Determining phased candidate set for gene HLA-$gene")

        val unphasedCandidates = aminoAcidSequences.filter { it.allele.asFourDigit() in unphasedCandidateAlleles }
        val phasedCandidates = filterCandidates(unphasedCandidates, phasedEvidence)
        logger.info("    ${phasedCandidates.size} candidates after phasing: " + phasedCandidates.map { it.allele }.joinToString(", "))

        return phasedCandidates.map { it.allele }
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