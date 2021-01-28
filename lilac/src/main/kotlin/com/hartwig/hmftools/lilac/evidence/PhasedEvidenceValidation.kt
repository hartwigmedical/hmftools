package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.seq.HlaSequence
import org.apache.logging.log4j.LogManager

class PhasedEvidenceValidation(val candidates: List<HlaSequence>, val expectedAlleles: List<HlaAllele>) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    fun validateEvidence(gene: String, evidence: List<PhasedEvidence>) {
        if (expectedAlleles.isNotEmpty()) {
            val expectedSequences = candidates.filter { it.allele.gene == gene && it.allele in expectedAlleles }
            for (sequence in expectedSequences) {
                for (phasedEvidence in evidence) {
                    if (!sequence.consistentWith(phasedEvidence)) {
                        logger.warn("Expected allele ${sequence.allele} filtered by $phasedEvidence")
                    }
                }
            }
        }

        for (inconsistentEvidence in unmatchedEvidence(evidence)) {
            logger.warn("Phased evidence not found in candidates: $inconsistentEvidence")
        }

    }

    fun unmatchedEvidence(evidence: List<PhasedEvidence>): List<PhasedEvidence> {
        return evidence.map { it.inconsistentEvidence(candidates) }.filter { it.evidence.isNotEmpty() }
    }


}