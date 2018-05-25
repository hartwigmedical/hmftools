package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GDnaVariant
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import org.apache.commons.csv.CSVRecord

data class CgiActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                               override val actionability: List<Actionability>, val cgiDrugs: List<CgiDrug>) : RecordMetadata by metadata,
        ActionableRecord {
    companion object {

        private val FUSION_SEPARATORS = setOf("__")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ABL1", "BCR"),
                                            FusionPair("PDGFRA", "FIP1L1"),
                                            FusionPair("PDGFB", "COL1A1"))
        private val FUSIONS_TO_FILTER = setOf(FusionPair("RET", "TPCN1"))
        private val fusionReader = FusionReader(separators = FUSION_SEPARATORS, filterSet = FUSIONS_TO_FILTER, flipSet = FUSIONS_TO_FLIP)

        operator fun invoke(record: CSVRecord): CgiActionableRecord? {
            val metadata = CgiMetadata(record["Gene"], record["transcript"] ?: "na")
            val events = readSomaticEvents(record)
            return CgiActionableRecord(metadata, events, readActionability(record), readCgiDrugs(record))
        }

        private fun readSomaticEvents(record: CSVRecord): List<SomaticEvent> {
            return readGdnaVariants(record) + listOfNotNull(readProteinAnnotation(record), readCdnaVariant(record), readCNV(record),
                                                            readFusion(record))
        }

        private fun readGdnaVariants(record: CSVRecord): List<GDnaVariant> {
            return record["gDNA"].orEmpty().split("__").map { it.trim() }.filterNot { it.isBlank() }.map { GDnaVariant(it) }
        }

        private fun readProteinAnnotation(record: CSVRecord): ProteinAnnotation? {
            val transcript = record["transcript"]
            val proteinAnnotation = record["individual_mutation"]?.substringAfter(':', "")
            return if (proteinAnnotation.isNullOrBlank() || transcript.isNullOrBlank()) {
                null
            } else {
                ProteinAnnotation(transcript, proteinAnnotation!!)
            }
        }

        private fun readCdnaVariant(record: CSVRecord): CDnaAnnotation? {
            val transcript = record["transcript"]
            val cdnaAnnotation = record["cDNA"]
            return if (cdnaAnnotation.isNullOrBlank() || transcript.isNullOrBlank()) {
                null
            } else {
                CDnaAnnotation(transcript, cdnaAnnotation!!)
            }
        }

        private fun readCNV(record: CSVRecord): CnvEvent? = when {
            record["Alteration type"] != "CNA"   -> null
            record["Alteration"].contains("amp") -> CnvEvent(record["Gene"], "Amplification")
            else                                 -> CnvEvent(record["Gene"], "Deletion")
        }

        private fun readFusion(record: CSVRecord): FusionEvent? = when {
            record["Alteration type"] == "FUS" -> fusionReader.read(record["Gene"], record["Alteration"])
            else                               -> null
        }

        private fun readActionability(record: CSVRecord): List<Actionability> {
            val cancerTypes = record["Primary Tumor type"].split(";").map { it.trim() }
            val level = record["Evidence level"]
            val association = record["Association"]
            return Actionability("cgi", cancerTypes, readDrugs(record), level, association, "Predictive",
                                 highestLevel(level), HmfResponse(association))
        }

        private fun readDrugs(record: CSVRecord): List<String> {
            val drugNames = readDrugsField(record["Drug"].orEmpty())
            return if (drugNames.isEmpty()) {
                readDrugsField(record["Drug family"].orEmpty())
            } else {
                drugNames
            }
        }

        private fun readDrugsField(drugField: String): List<String> {
            return drugField.replace(";", " + ")
                    .removeSurrounding("[", "]")
                    .split(",")
                    .filterNot { it.isBlank() }
        }

        private fun highestLevel(levelField: String): HmfLevel {
            return levelField.split(",").map { it.trim() }.map { HmfLevel(it) }.sorted().firstOrNull() ?: HmfLevel.UNKNOWN
        }
    }
}
