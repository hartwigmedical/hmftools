package com.hartwig.hmftools.lilac.candidates

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.nuc.NucleotideFiltering
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.phase.PhasedEvidence
import com.hartwig.hmftools.lilac.phase.TypeEvidence
import com.hartwig.hmftools.lilac.seq.HlaSequence
import org.apache.logging.log4j.LogManager

class Candidates(private val minBaseCount: Int, private val minFragmentCount: Int, private val nucleotideSequences: List<HlaSequence>, private val aminoAcidSequences: List<HlaSequence>) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    fun candidates(gene: String, aminoAcidBoundary: Set<Int>, nucleotideFragments: List<NucleotideFragment>): List<HlaSequence> {
        val aminoAcidFragments = nucleotideFragments.map { it.toAminoAcidFragment() }

        logger.info("Determining initial candidate set for gene HLA-$gene")
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseCount, aminoAcidFragments)
        val nucleotideCounts = SequenceCount.nucleotides(minBaseCount, aminoAcidFragments)

        val geneCandidates = aminoAcidSequences.filter { it.allele.gene == gene }
        logger.info(" ... ${geneCandidates.size} candidates before filtering")

        // Nucleotide filtering
        val nucleotideFiltering = NucleotideFiltering(minBaseCount, aminoAcidBoundary)
        val nucleotideSpecificProteinCandidates = nucleotideFiltering.filterCandidatesOnAminoAcidBoundaries(nucleotideSequences, aminoAcidFragments).map { it.allele.specificProtein() }.toSet()
        val nucleotideCandidates = geneCandidates.filter { it.allele.specificProtein() in nucleotideSpecificProteinCandidates }
        logger.info(" ... ${nucleotideCandidates.size} candidates after nucleotide exon boundary filtering")

        // Amino acid filtering
        val aminoAcidCandidates = initialCandidates(aminoAcidCounts, nucleotideCandidates)
        logger.info(" ... ${aminoAcidCandidates.size} candidates after amino acid filtering")

        val typeEvidenceFactory = TypeEvidence(minBaseCount, minFragmentCount)
        val typeEvidence = typeEvidenceFactory.evidence(aminoAcidFragments)

        val result = filterCandidates(aminoAcidCandidates, typeEvidence)
        logger.info(" ... ${result.size} candidates after phasing: " + result.map { it.allele }.joinToString(", "))
        return result

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


    private fun initialCandidates(aminoAcidCount: SequenceCount, candidates: List<HlaSequence>): List<HlaSequence> {
        var result = candidates
        val locations = (0 until aminoAcidCount.length).toSet()
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
            candidates = matchingCandidates(newEvidence, candidates)
        }

        return candidates
    }

    private fun matchingCandidates(evidence: PhasedEvidence, candidates: Collection<HlaSequence>): List<HlaSequence> {
        return candidates.filter { it.consistentWith(evidence) }
    }


}