package com.hartwig.hmftools.knowledgebaseimporter.output

enum class HmfResponse {
    Resistant, Responsive, OTHER;

    companion object {
        private val responseMap = mapOf("Resistance or Non-Response" to Resistant,
                                        "Sensitivity" to Responsive,
                                        "Resistant" to Resistant,
                                        "Responsive" to Responsive)

        operator fun invoke(response: String): HmfResponse {
            return responseMap[response] ?: OTHER
        }
    }
}
