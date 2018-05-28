package com.hartwig.hmftools.knowledgebaseimporter.output

data class Actionability(val source: String, val cancerType: String, val drug: HmfDrug, val level: String, val significance: String,
                         val evidenceType: String, val hmfLevel: HmfLevel, val hmfResponse: HmfResponse) {
    companion object {
        val header = listOf("source", "drug", "drugType", "cancerType", "level", "hmfLevel", "evidenceType", "significance", "hmfResponse")

        operator fun invoke(source: String, cancerTypes: List<String>, drugs: List<HmfDrug>, level: String, significance: String,
                            evidenceType: String, hmfLevel: HmfLevel, hmfResponse: HmfResponse): List<Actionability> {
            return cancerTypes.flatMap { cancerType ->
                drugs.map { drug ->
                    Actionability(source, cancerType, drug, level, significance, evidenceType, hmfLevel, hmfResponse)
                }
            }
        }
    }

    val record: List<String> = listOf(source, drug.name, drug.type, cancerType, level, hmfLevel.name, evidenceType, significance,
                                      hmfResponse.name)
}
