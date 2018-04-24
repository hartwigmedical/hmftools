package com.hartwig.hmftools.knowledgebaseimporter.cgi

import org.apache.commons.csv.CSVRecord

data class CgiBiomarkersRecord(private val csvRecord: CSVRecord) {
    val gene: String = csvRecord["Gene"]
    val gdna: String = csvRecord["gDNA"] ?: "na"
    val cdna: String = csvRecord["cDNA"] ?: "na"
    val transcript: String = csvRecord["transcript"] ?: "na"
    val protein: String = csvRecord["individual_mutation"]?.substringAfter(':', "") ?: "na"
    val drug: String = csvRecord["Drug full name"]
    val level: String = csvRecord["Evidence level"]
    val association: String = csvRecord["Association"]
    val cancerTypes = csvRecord["Primary Tumor type"].split(";").map { it.trim() }
}
