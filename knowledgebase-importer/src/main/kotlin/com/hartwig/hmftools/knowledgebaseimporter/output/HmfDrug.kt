package com.hartwig.hmftools.knowledgebaseimporter.output

data class HmfDrug(val name: String, val type: String) {
    companion object {
        val header = listOf("name", "type")
    }

    val record: List<String> = listOf(name, type)
}
