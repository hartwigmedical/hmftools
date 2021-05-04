package com.hartwig.hmftools.lilackt.evidence

import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci
import org.apache.logging.log4j.LogManager

object PhasedEvidenceValidation {

    val logger = LogManager.getLogger(this::class.java)

    fun validateExpected(gene: String, evidence: List<PhasedEvidence>, candidates: List<HlaSequenceLoci>) {
        val expectedSequences = candidates.filter { it.allele.gene == gene }
        for (sequence in expectedSequences) {
            for (phasedEvidence in evidence) {
                if (!sequence.consistentWith(phasedEvidence)) {
                    logger.warn("Expected allele ${sequence.allele} filtered by $phasedEvidence")
                }
            }
        }
    }

    fun validateAgainstFinalCandidates(gene: String, evidence: List<PhasedEvidence>, candidates: List<HlaSequenceLoci>) {
        for (inconsistentEvidence in unmatchedEvidence(evidence, candidates)) {
            logger.warn("HLA-$gene phased evidence not found in candidates: $inconsistentEvidence")
        }
    }


    private fun unmatchedEvidence(evidence: List<PhasedEvidence>, candidates: List<HlaSequenceLoci>): List<PhasedEvidence> {
        return evidence.map { it.inconsistentEvidence(candidates) }.filter { it.evidence.isNotEmpty() }
    }


}