package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.extensions.csv.CsvData

data class Actionability(val source: String, val reference: String, val drug: HmfDrug, val cancerType: String, val level: String,
                         val hmfLevel: String, val evidenceType: String, val significance: String, val hmfResponse: String) : CsvData {
    companion object {
        operator fun invoke(source: String, reference: String, cancerTypes: List<String>, drugs: List<HmfDrug>, level: String,
                            significance: String, evidenceType: String, hmfLevel: HmfLevel, hmfResponse: HmfResponse): List<Actionability> {
            return cancerTypes.flatMap { cancerType ->
                drugs.map { drug ->
                    Actionability(source, reference, drug, cancerType, level, hmfLevel.toString(), evidenceType, significance,
                                  hmfResponse.toString())
                }
            }
        }
    }
}
