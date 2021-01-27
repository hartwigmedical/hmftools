package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.seq.HlaSequence
import org.apache.logging.log4j.LogManager

class PhasedEvidenceValidation(val candidates: List<HlaSequence>) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    fun validateEvidence(evidence: List<PhasedEvidence>) {
        for (inconsistentEvidence in unmatchedEvidence(evidence)) {
            logger.warn("Phased evidence not found in candidates: $inconsistentEvidence")
        }

    }

    fun unmatchedEvidence(evidence: List<PhasedEvidence>): List<PhasedEvidence> {
        return evidence.map { it.inconsistentEvidence(candidates) }.filter { it.evidence.isNotEmpty() }
    }


}