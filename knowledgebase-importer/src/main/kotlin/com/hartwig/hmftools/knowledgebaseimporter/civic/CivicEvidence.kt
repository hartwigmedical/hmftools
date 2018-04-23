package com.hartwig.hmftools.knowledgebaseimporter.civic

import org.apache.commons.csv.CSVRecord

data class CivicEvidence(private val csvRecord: CSVRecord) {
    val cancerType: String = csvRecord["disease"]
    val doid: String = csvRecord["doid"]
    val drugs: List<String> = csvRecord["drugs"].split(",").map { it.trim() }
    val evidenceType: String = csvRecord["evidence_type"]
    val evidenceDirection: String = csvRecord["evidence_direction"]
    val evidenceLevel: String = csvRecord["evidence_level"]
    val significance: String = csvRecord["clinical_significance"]
}