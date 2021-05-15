package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
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
}