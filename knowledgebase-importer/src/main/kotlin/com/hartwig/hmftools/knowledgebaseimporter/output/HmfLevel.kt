package com.hartwig.hmftools.knowledgebaseimporter.output

enum class HmfLevel {
    A, B, C, D, E, UNKNOWN;

    companion object {
        private val evidenceLevelMap = mapOf("FDA guidelines" to A,
                                             "NCCN guidelines" to B,
                                             "NCCN/CAP guidelines" to B,
                                             "CPIC guidelines" to B,
                                             "European LeukemiaNet guidelines" to B,
                                             "Clinical trials" to B,
                                             "Late trials" to B,
                                             "Early trials" to C,
                                             "Case report" to C,
                                             "Pre-clinical" to D,
                                             "1" to A,
                                             "2A" to B,
                                             "2B" to B,
                                             "3A" to C,
                                             "3B" to C,
                                             "4" to D,
                                             "R1" to B,
                                             "R2" to C,
                                             "R3" to D,
                                             "A" to A,
                                             "B" to B,
                                             "C" to C,
                                             "D" to D,
                                             "E" to E)

        operator fun invoke(evidenceLevel: String): HmfLevel {
            return evidenceLevelMap[evidenceLevel] ?: UNKNOWN
        }
    }
}
