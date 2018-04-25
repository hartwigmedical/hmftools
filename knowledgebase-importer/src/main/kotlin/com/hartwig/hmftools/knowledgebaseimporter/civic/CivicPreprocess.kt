package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords

fun preProcessCivicRecords(variantFileLocation: String, evidenceFileLocation: String): List<CivicRecord> {
    val variantEvidenceMap = readEvidenceMap(evidenceFileLocation)
    return readCSVRecords(variantFileLocation) { CivicRecord(it, variantEvidenceMap) }
}

fun preProcessCivicVariants(variantFileLocation: String, evidenceFileLocation: String): List<CivicRecord> {
    return preProcessCivicRecords(variantFileLocation, evidenceFileLocation)
            .filter { (!it.chromosome.isBlank() && !it.start.isBlank() && !it.stop.isBlank() && (!it.ref.isBlank() || !it.alt.isBlank())) || (!it.hgvs.isBlank()) }
}

private fun readEvidenceMap(evidenceFileLocation: String): Multimap<String, CivicEvidence> {
    val evidenceMap = ArrayListMultimap.create<String, CivicEvidence>()
    readCSVRecords(evidenceFileLocation) { csvRecord -> evidenceMap.put(csvRecord["variant_id"], CivicEvidence(csvRecord)) }
    return evidenceMap
}
