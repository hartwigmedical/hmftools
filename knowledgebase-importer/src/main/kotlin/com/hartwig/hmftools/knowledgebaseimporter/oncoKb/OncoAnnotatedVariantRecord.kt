package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import org.apache.commons.csv.CSVRecord

data class OncoAnnotatedVariantRecord(private val csvRecord: CSVRecord) {
    val transcript: String = csvRecord[0]
    val gene: String = csvRecord[3]
    val impact: String = csvRecord[4]
    val oncogenicity: String = csvRecord[5]
}
