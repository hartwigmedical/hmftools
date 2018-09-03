package com.hartwig.hmftools.knowledgebaseimporter.civic.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class CivicEvidenceInput(val evidence_id: String, val variant_id: String, val evidence_type: String, val evidence_direction: String,
                              val evidence_level: String, val drugs: String, val disease: String, val doid: String,
                              val clinical_significance: String) : CsvData, CorrectedInput<CivicEvidenceInput> {

    //MIVO: map melanoma to skin melanoma
    override fun correct(): CivicEvidenceInput? {
        if (evidence_id == "1481") return null
        val disease = if (disease == "Melanoma") "Skin Melanoma" else disease
        val doid = if (doid == "1909") "8923" else doid
        return copy(disease = disease, doid = doid)
    }
}
