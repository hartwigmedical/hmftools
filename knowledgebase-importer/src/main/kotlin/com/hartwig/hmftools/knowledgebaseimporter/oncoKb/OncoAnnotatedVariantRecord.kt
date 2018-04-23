package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import org.apache.commons.csv.CSVRecord

data class OncoAnnotatedVariantRecord(private val csvRecord: CSVRecord) {
    val transcript: String = csvRecord["Isoform"]
    val gene: String = csvRecord["Gene"]
    val alteration: String = csvRecord["Alteration"]
    val oncogenicity: String = csvRecord["Oncogenicity"]
}
