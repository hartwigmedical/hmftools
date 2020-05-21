package com.hartwig.hmftools.gripss

data class GripssConfig(
        val inputVcf: String,
        val outputVcf: String,
        val singleBreakendPonFile: String,
        val pairedBreakpointPonFile: String,
        val filterConfig: GripssFilterConfig)

data class GripssFilterConfig(
        val maxNormalSupport: Double,
        val minNormalCoverage: Int,
        val minTumorAF: Double,
        val maxShortStrandBias: Double,
        val minQualBreakEnd: Int,
        val minQualBreakPoint: Int,
        val maxHomLength: Int,
        val maxHomLengthShortInversion: Int,
        val maxInexactHomLength: Int,
        val maxInexactHomLengthShortDel: Int,
        val minSize: Int) {

    companion object {
        fun default(): GripssFilterConfig {
            return GripssFilterConfig(
                    0.03,
                    8,
                    0.005,
                    0.95,
                    1000,
                    350,
                    50,
                    6,
                    50,
                    5,
                    32)
        }
    }

}