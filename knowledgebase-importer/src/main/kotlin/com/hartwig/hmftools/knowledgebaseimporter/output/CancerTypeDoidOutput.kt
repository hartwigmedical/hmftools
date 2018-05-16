package com.hartwig.hmftools.knowledgebaseimporter.output

data class CancerTypeDoidOutput(val cancerType: String, val doidSet: String) {
    companion object {
        val header = listOf("cancerType", "doids")
    }

    val record: List<String> = listOf(cancerType, doidSet)
}
