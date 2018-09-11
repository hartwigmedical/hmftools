package com.hartwig.hmftools.knowledgebaseimporter.output

enum class HmfResponse {
    Resistant, Responsive, OTHER;

    companion object {
        private val responseMap = mapOf("Resistance or Non-Response" to Resistant,
                                        "Resistant" to Resistant,
                                        "Resistance" to Resistant,
                                        "Sensitivity" to Responsive,
                                        "Sensitivity/Response" to Responsive,
                                        "Responsive" to Responsive)

        operator fun invoke(response: String): HmfResponse {
            return responseMap[response] ?: OTHER
        }
    }
}
