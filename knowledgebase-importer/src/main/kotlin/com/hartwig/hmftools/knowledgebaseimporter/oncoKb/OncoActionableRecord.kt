package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.extractFusion
import com.hartwig.hmftools.knowledgebaseimporter.flipFusion
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import org.apache.commons.csv.CSVRecord

data class OncoActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                                override val actionability: List<Actionability>) : RecordMetadata by metadata, ActionableRecord {
    companion object {
        private val FUSION_SEPARATORS = listOf("-", " - ", "?")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ROS1", "CD74"),
                                            FusionPair("EP300", "MLL"),
                                            FusionPair("EP300", "MOZ"),
                                            FusionPair("RET", "CCDC6"))
        private val FUSIONS_TO_FILTER = listOf<FusionPair>()

        operator fun invoke(csvRecord: CSVRecord): OncoActionableRecord {
            val metadata = OncoKbMetadata(csvRecord["Gene"], csvRecord["Isoform"])
            val level: String = readLevel(csvRecord["Level"])
            val significance = if (csvRecord["Level"].startsWith("R")) HmfResponse.Resistant else HmfResponse.Responsive
            val drugs = csvRecord["Drugs(s)"].split(",").map { it.trim() }
            val cancerType: String = csvRecord["Cancer Type"]
            val actionability = Actionability("oncoKb", listOf(cancerType), drugs, level, significance.name,
                                              "Predictive", HmfLevel(csvRecord["Level"]), significance)
            return OncoActionableRecord(metadata, readSomaticEvents(csvRecord), actionability)
        }

        private fun readLevel(levelField: String): String {
            return if (levelField.startsWith("R")) levelField.drop(1) else levelField
        }

        private fun readSomaticEvents(record: CSVRecord): List<SomaticEvent> {
            val fusionsAndCnvs = listOfNotNull(readFusions(record), readCnv(record))
            return if (fusionsAndCnvs.isEmpty()) {
                listOfNotNull(ProteinAnnotation(record["Isoform"], record["Alteration"]))
            } else {
                fusionsAndCnvs
            }
        }

        private fun readFusions(record: CSVRecord): FusionEvent? {
            val alteration = record["Alteration"]
            val fusion = when {
                alteration == null                   -> null
                alteration.contains(Regex("Fusion")) -> extractFusion(record["Gene"], alteration.trim(), FUSION_SEPARATORS)
                else                                 -> null
            }
            return if (FUSIONS_TO_FILTER.contains(fusion) || fusion == null) null else flipFusion(fusion, FUSIONS_TO_FLIP)
        }

        private fun readCnv(record: CSVRecord): CnvEvent? {
            return when {
                record["Alteration"] == "Amplification" -> CnvEvent(record["Gene"], "Amplification")
                record["Alteration"] == "Deletion"      -> CnvEvent(record["Gene"], "Deletion")
                else                                    -> null
            }
        }
    }
}
