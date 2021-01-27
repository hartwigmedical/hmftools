package com.hartwig.hmftools.lilac.candidates

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceFactory
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.nuc.NucleotideFiltering
import com.hartwig.hmftools.lilac.seq.HlaSequence
import org.apache.logging.log4j.LogManager

class Candidates(
        private val minFragmentsPerAllele: Int,
        private val minFragmentsToRemoveSingles: Int,
        private val minBaseCount: Int,
        private val nucleotideSequences: List<HlaSequence>,
        private val aminoAcidSequences: List<HlaSequence>) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    fun candidates(context: HlaContext, aminoAcidFragments: List<AminoAcidFragment>): List<HlaSequence> {
        val gene = context.gene
        val aminoAcidBoundary = context.aminoAcidBoundaries
        val expectedAlleles = context.expectedAlleles

        logger.info("Determining initial candidate set for gene HLA-$gene")
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseCount, aminoAcidFragments)
//        aminoAcidCounts.writeVertically("/Users/jon/hmf/analysis/hla/output/aminoacids.${gene}.count.txt")
//        val nucleotideCounts = SequenceCount.nucleotides(minBaseCount, aminoAcidFragments)

        val geneCandidates = aminoAcidSequences.filter { it.allele.gene == gene }
        logger.info(" ... ${geneCandidates.size} candidates before filtering")

        // Amino acid filtering
        val aminoAcidCandidates = aminoAcidCandidates(aminoAcidBoundary, aminoAcidCounts, geneCandidates)
        val aminoAcidCandidateAlleles = aminoAcidCandidates.map { it.allele }.toSet()
        val aminoAcidSpecificAllelesCandidate = aminoAcidCandidateAlleles.map { it.specificProtein() }.toSet()

        logger.info(" ... ${aminoAcidCandidates.size} candidates after amino acid filtering")

        // Nucleotide filtering
        val nucleotideFiltering = NucleotideFiltering(minBaseCount, aminoAcidBoundary)
        val nucleotideCandidatesAfterAminoAcidFiltering = nucleotideSequences
                .filter { it.allele.specificProtein() in aminoAcidSpecificAllelesCandidate }
        val nucleotideSpecificAllelesCandidate = nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(nucleotideCandidatesAfterAminoAcidFiltering, aminoAcidFragments)
                .map { it.allele.specificProtein() }
                .toSet()

        val nucleotideCandidates = aminoAcidCandidates.filter { it.allele.specificProtein() in nucleotideSpecificAllelesCandidate }
        logger.info(" ... ${nucleotideCandidates.size} candidates after exon boundary filtering")

        val phasedEvidenceFactory = PhasedEvidenceFactory(minFragmentsPerAllele, minFragmentsToRemoveSingles, minBaseCount)
        val phasedEvidence = phasedEvidenceFactory.evidence(expectedAlleles, aminoAcidFragments)

        val phasedCandidates = filterCandidates(nucleotideCandidates, phasedEvidence)
        logger.info(" ... ${phasedCandidates.size} candidates after phasing: " + phasedCandidates.map { it.allele }.joinToString(", "))

        return phasedCandidates

    }

    private fun checkCandidates(candidates: Collection<HlaSequence>): Int {
        var count = 0

        if (candidates.any { it.allele == HlaAllele("A*01:01:01:01") }) {
            count++;
        }

        if (candidates.any { it.allele == HlaAllele("A*11:01:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("B*08:01:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("B*56:01:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("C*01:02:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("C*07:01:01:01") }) {
            count++;
        }


        return count;
    }

    private fun checkColo8289Candidates(candidates: Collection<HlaSequence>): Int {
        var count = 0

        if (candidates.any { it.allele == HlaAllele("A*01:01:01:01") }) {
            count++;
        }

        if (candidates.any { it.allele == HlaAllele("C*03:04:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("C*07:01:01:01") }) {
            count++;
        }

        if (candidates.any { it.allele == HlaAllele("B*08:01:01:01") }) {
            count++;
        }

        if (candidates.any { it.allele == HlaAllele("B*40:02:01:01") }) {
            count++;
        }

        return count;
    }


    private fun aminoAcidCandidates(boundaries: Set<Int>, aminoAcidCount: SequenceCount, candidates: List<HlaSequence>): List<HlaSequence> {
        var result = candidates
//        val locations = (0 until aminoAcidCount.length).filter { aminoAcidCount.depth(it) >= 12 }
        val locations = (0 until aminoAcidCount.length).toSet() subtract boundaries
        for (location in locations) {
            result = filterCandidates(location, aminoAcidCount.sequenceAt(location), result)
        }
        return result
    }

    private fun filterCandidates(index: Int, expectedCharacters: Collection<Char>, candidates: Collection<HlaSequence>): List<HlaSequence> {
        return candidates.filter { it.length <= index || it.sequence[index] == '*' || it.sequence[index] in expectedCharacters }
    }

    private fun filterCandidates(initialCandidates: List<HlaSequence>, evidence: List<PhasedEvidence>): List<HlaSequence> {
        var candidates = initialCandidates
        for (i in evidence.indices) {
            val newEvidence = evidence[i]
            logger.info(" ... $newEvidence")
            candidates = matchingCandidates(newEvidence, candidates)
//            logger.info(" ... contains ${(checkColo8289Candidates(candidates))} COLO829 ")

        }

        return candidates
    }

    private fun matchingCandidates(evidence: PhasedEvidence, candidates: Collection<HlaSequence>): List<HlaSequence> {
        return candidates.filter { it.consistentWith(evidence) }
    }


}