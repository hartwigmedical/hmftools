package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import org.apache.commons.csv.CSVRecord

data class CgiKnownKbRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                            override val events: List<SomaticEvent>) : RecordMetadata by metadata, KnownRecord {
    companion object {
        operator fun invoke(csvRecord: CSVRecord): CgiKnownKbRecord? {
            if (csvRecord["context"] != "somatic") return null
            val metadata = CgiMetadata(csvRecord["gene"], csvRecord["transcript"])
            val events = readSomaticEvents(csvRecord)
            return CgiKnownKbRecord(metadata, "", events)
        }

        private fun readSomaticEvents(csvRecord: CSVRecord): List<SomaticEvent> {
            val gDnaVariants = csvRecord["gdna"].orEmpty().split("__").map { it.trim() }.filterNot { it.isBlank() }.map { GDnaVariant(it) }
            val proteinAnnotations = listOfNotNull(csvRecord["protein"]).filterNot { it.isBlank() }.map {
                ProteinAnnotation(csvRecord["transcript"], it)
            }
            return gDnaVariants + proteinAnnotations
        }
    }
}
