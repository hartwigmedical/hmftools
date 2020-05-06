package com.hartwig.hmftools.gripss

data class GripssConfig(val inputVcf: String, val outputVcf: String, val filterConfig: GripssFilterConfig)


data class GripssFilterConfig(val maxNormalSupport: Double, val minNormalCoverage: Int, val minTumorAF: Double, val maxShortStrandBias: Double)