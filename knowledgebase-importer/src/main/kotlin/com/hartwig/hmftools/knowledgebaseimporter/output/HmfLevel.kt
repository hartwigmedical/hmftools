package com.hartwig.hmftools.knowledgebaseimporter.output

enum class HmfLevel {
    A, B, C, D, E, UNKNOWN;

    companion object {
        private val evidenceLevelMap = mapOf("FDA guidelines" to A,
                                             "NCCN guidelines" to A,
                                             "NCCN/CAP guidelines" to A,
                                             "CPIC guidelines" to A,
                                             "European LeukemiaNet guidelines" to A,
                                             "Clinical trials" to B,
                                             "Late trials" to B,
                                             "Late trials,Pre-clinical" to B,
                                             "Early trials" to B,
                                             "Case report" to C,
                                             "Early Trials,Case Report" to C,
                                             "Pre-clinical" to D,
                                             "1" to A,
                                             "2A" to A,
                                             "2B" to A,
                                             "3A" to B,
                                             "3B" to B,
                                             "4" to D,
                                             "R1" to A,
                                             "R2" to B,
                                             "R3" to D,
                                             "A" to A,
                                             "B" to B,
                                             "C" to C,
                                             "D" to D,
                                             "E" to E)

        operator fun invoke(evidenceLevel: String): HmfLevel {
            return evidenceLevelMap[evidenceLevel.trim()] ?: UNKNOWN
        }
    }
}
