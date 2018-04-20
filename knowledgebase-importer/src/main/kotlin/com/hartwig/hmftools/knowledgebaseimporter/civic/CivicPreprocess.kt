package com.hartwig.hmftools.knowledgebaseimporter.civic

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

private val format = CSVFormat.TDF.withFirstRecordAsHeader()

fun preProcessCivic(variantFileLocation: String, evidenceFileLocation: String): List<CivicRecord> {
    val variantEvidenceMap = readEvidenceMap(evidenceFileLocation)
    val parser = CSVParser.parse(File(variantFileLocation), Charset.defaultCharset(), format)
    return parser.iterator().asSequence().map { CivicRecord(it, variantEvidenceMap) }
            .filter { (!it.chromosome.isBlank() && !it.start.isBlank() && !it.stop.isBlank() && (!it.ref.isBlank() || !it.alt.isBlank())) || (!it.hgvs.isBlank()) }
            .toList()
}

private fun highestEvidence(
        variantEvidencePairs: List<Pair<String, String>>): String = variantEvidencePairs.map { it.second }.sorted().first()

private fun readEvidenceMap(evidenceFileLocation: String): Map<String, String> {
    val evidenceParser = CSVParser.parse(File(evidenceFileLocation), Charset.defaultCharset(),
                                         format)
    return evidenceParser.iterator().asSequence()
            .map { csvRecord -> csvRecord.get("variant_id") to csvRecord.get("evidence_level") }
            .groupBy { it.first }
            .mapValues { highestEvidence(it.value) }
            .toMap()
}
