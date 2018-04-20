package com.hartwig.hmftools.knowledgebaseimporter.cgi

import org.apache.commons.csv.CSVRecord

data class CgiRecord(private val csvRecord: CSVRecord) {
    val gene: String = csvRecord["gene"]
    val gdna: String = csvRecord["gdna"]
    val impact: String = csvRecord["protein"]
    val transcript: String = csvRecord["transcript"]
    val context: String = csvRecord["context"]
}
