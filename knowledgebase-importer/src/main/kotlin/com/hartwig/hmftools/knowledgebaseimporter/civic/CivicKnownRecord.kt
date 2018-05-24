package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.Multimap
import com.hartwig.hmftools.knowledgebaseimporter.extractFusion
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import org.apache.commons.csv.CSVRecord
import java.util.regex.Pattern

data class CivicKnownRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                            override val events: List<SomaticEvent>, override val actionability: List<Actionability>,
                            val cancerDoids: Map<String, String>) :
        RecordMetadata by metadata, KnownRecord, ActionableRecord {
    companion object {
        private val FUSION_SEPARATORS = listOf("-")
        private val FUSIONS_TO_FILTER = setOf(FusionPair("BRAF", "CUL1"))

        operator fun invoke(csvRecord: CSVRecord, variantEvidenceMap: Multimap<String, CivicEvidence>): CivicKnownRecord {
            val metadata = CivicMetadata(csvRecord["gene"], csvRecord["representative_transcript"] ?: "-")
            val evidence = variantEvidenceMap[csvRecord["variant_id"]]
            val additionalInfo = additionalInfo(evidence)
            val actionability = evidence.filter { it.direction == "Supports" }.flatMap { it.actionabilityItems }
            val doids = evidence.associateBy({ it.cancerType }, { it.doid })
            return CivicKnownRecord(metadata, additionalInfo, readSomaticEvents(csvRecord), actionability, doids)
        }

        private fun readSomaticEvents(record: CSVRecord): List<SomaticEvent> {
            val events = mutableListOf<SomaticEvent>()
            if (hasVariant(record)) {
                events.add(KnowledgebaseVariant(record["gene"], record["chromosome"], record["start"].toLong(), record["reference_bases"],
                                                record["variant_bases"]))
            }
            if (hasHgvs(record)) {
                val hgvsParts = extractHgvs(record).split(":")
                events.add(CDnaAnnotation(hgvsParts[0], hgvsParts[1]))
            }
            return events + listOfNotNull(readCNV(record), readFusion(record))
        }

        private fun hasVariant(record: CSVRecord): Boolean {
            return (!record["chromosome"].isNullOrBlank() && !record["start"].isNullOrBlank() && !record["stop"].isNullOrBlank() &&
                    (!record["reference_bases"].isNullOrBlank() || !record["variant_bases"].isNullOrBlank()))
        }

        private fun hasHgvs(record: CSVRecord): Boolean {
            return extractHgvs(record).isNotBlank()
        }

        private fun extractHgvs(record: CSVRecord): String {
            val hgvsField = record["hgvs_expressions"].orEmpty()
            val matcher = Pattern.compile("(ENST[0-9]+\\.[0-9+]:c\\.[0-9][^,\\t\\s\\n]+)").matcher(hgvsField)
            return if (matcher.find()) matcher.group(1) else ""
        }

        private fun additionalInfo(evidence: Collection<CivicEvidence>): String {
            val highestEvidenceLevel = evidence.map { it.level }.sorted().firstOrNull() ?: "N"
            return (highestEvidenceLevel == "A" || highestEvidenceLevel == "B" || highestEvidenceLevel == "C").toString()
        }

        private fun readCNV(record: CSVRecord): CnvEvent? {
            return when {
                record["variant"] == "AMPLIFICATION" -> CnvEvent(record["gene"], "Amplification")
                record["variant"] == "DELETION"      -> CnvEvent(record["gene"], "Deletion")
                record["variant"] == "LOH"           -> CnvEvent(record["gene"], "Deletion")
                else                                 -> null
            }
        }

        private fun readFusion(record: CSVRecord): FusionEvent? {
            val variantTypes = record["variant_types"].orEmpty()
            val fusion = when {
                variantTypes.contains("fusion") -> extractCorrectFusion(record["gene"], record["variant"].trim())
                else                            -> null
            }
            return if (FUSIONS_TO_FILTER.contains(fusion) || fusion == null) null else fusion
        }

        private fun extractCorrectFusion(gene: String, variant: String): FusionEvent {
            return if (variant.contains(Regex("MLL-MLLT3")) && gene == "KMT2A") {
                extractFusion(gene, variant.replace("MLL-MLLT3", "KMT2A-MLLT3"), FUSION_SEPARATORS)
            } else {
                extractFusion(gene, variant, FUSION_SEPARATORS)
            }
        }
    }
}
