package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import org.apache.commons.csv.CSVRecord

data class OncoAnnotatedVariantRecord(val transcript: String, val gene: String, val alteration: String, val oncogenicity: String) {
    companion object {
        operator fun invoke(csvRecord: CSVRecord): OncoAnnotatedVariantRecord {
            return OncoAnnotatedVariantRecord(csvRecord["Isoform"], csvRecord["Gene"], csvRecord["Alteration"], csvRecord["Oncogenicity"])
        }
    }
}
