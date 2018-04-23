package com.hartwig.hmftools.knowledgebaseimporter.civic

import org.apache.commons.csv.CSVRecord

data class CivicEvidence(private val csvRecord: CSVRecord) {
    val cancerType: String = csvRecord["disease"]
    val doid: String = csvRecord["doid"]
    val drugs: List<String> = csvRecord["drugs"].split(",").map { it.trim() }
    val type: String = csvRecord["evidence_type"]
    val direction: String = csvRecord["evidence_direction"]
    val level: String = csvRecord["evidence_level"]
    val significance: String = csvRecord["clinical_significance"]
}