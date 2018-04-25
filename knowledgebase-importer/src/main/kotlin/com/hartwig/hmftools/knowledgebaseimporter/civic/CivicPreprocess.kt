package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

private val format = CSVFormat.TDF.withFirstRecordAsHeader()

fun preProcessCivicRecords(variantFileLocation: String, evidenceFileLocation: String): List<CivicRecord> {
    val variantEvidenceMap = readEvidenceMap(evidenceFileLocation)
    val parser = CSVParser.parse(File(variantFileLocation), Charset.defaultCharset(), format)
    return parser.iterator().asSequence().map { CivicRecord(it, variantEvidenceMap) }
            .toList()
}

fun preProcessCivicVariants(variantFileLocation: String, evidenceFileLocation: String): List<CivicRecord> {
    return preProcessCivicRecords(variantFileLocation, evidenceFileLocation)
            .filter { (!it.chromosome.isBlank() && !it.start.isBlank() && !it.stop.isBlank() && (!it.ref.isBlank() || !it.alt.isBlank())) || (!it.hgvs.isBlank()) }
}

private fun readEvidenceMap(evidenceFileLocation: String): Multimap<String, CivicEvidence> {
    val evidenceMap = ArrayListMultimap.create<String, CivicEvidence>()
    val evidenceParser = CSVParser.parse(File(evidenceFileLocation), Charset.defaultCharset(), format)
    evidenceParser.iterator().asSequence()
            .forEach { csvRecord -> evidenceMap.put(csvRecord["variant_id"], CivicEvidence(csvRecord)) }
    return evidenceMap
}
