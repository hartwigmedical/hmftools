package com.hartwig.hmftools.knowledgebaseimporter.civic.input

import com.hartwig.hmftools.extensions.csv.CsvData

data class CivicEvidenceInput(val evidence_id: String, val variant_id: String, val evidence_type: String, val evidence_direction: String,
                              val evidence_level: String, val drugs: String, private val disease: String, private val doid: String,
                              val clinical_significance: String) : CsvData {

    //MIVO: map melanoma to skin melanoma
    val cancerType = if (disease == "Melanoma") "Skin Melanoma" else disease
    val cancerDoid = if (doid == "1909") "8923" else doid
}
