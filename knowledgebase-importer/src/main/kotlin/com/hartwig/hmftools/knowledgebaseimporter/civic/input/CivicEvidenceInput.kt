package com.hartwig.hmftools.knowledgebaseimporter.civic.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class CivicEvidenceInput(val evidence_id: String, val variant_id: String, val evidence_type: String, val evidence_direction: String,
                              val evidence_level: String, val drugs: String, val disease: String, val doid: String,
                              val clinical_significance: String) : CsvData, CorrectedInput<CivicEvidenceInput> {

    override fun correct(): CivicEvidenceInput? {
        // MIVO: map melanoma to skin melanoma
        val disease = if (disease == "Melanoma") "Skin Melanoma" else disease
        // KODU: Also map melanoma DOID to skin melanoma DOID
        val doid = if (doid == "1909") "8923" else doid
        return copy(disease = disease, doid = doid)
    }
}
