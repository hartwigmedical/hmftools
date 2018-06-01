package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.Multimap
import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import org.apache.commons.csv.CSVRecord
import java.util.regex.Pattern

data class CivicRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                       override val events: List<SomaticEvent>, override val actionability: List<Actionability>,
                       val cancerDoids: Map<String, String>) :
        RecordMetadata by metadata, KnownRecord, ActionableRecord {
    companion object {
        private val FUSION_SEPARATORS = setOf("-")
        private val FUSIONS_TO_FILTER = setOf(FusionPair("BRAF", "CUL1"))
        private val fusionReader = FusionReader(separators = FUSION_SEPARATORS, filterSet = FUSIONS_TO_FILTER)

        operator fun invoke(record: CSVRecord, variantEvidenceMap: Multimap<String, CivicEvidence>): CivicRecord {
            val metadata = CivicMetadata(record["gene"], record["representative_transcript"] ?: "-")
            val evidence = variantEvidenceMap[record["variant_id"]]
            val additionalInfo = additionalInfo(evidence)
            val actionability = evidence.filter { it.direction == "Supports" }.flatMap { it.actionabilityItems }
            val doids = evidence.associateBy({ it.cancerType }, { it.doid })
            val (gene, variant) = correctRecord(record["gene"], record["variant"])
            return CivicRecord(metadata, additionalInfo, readSomaticEvents(gene, variant, record), actionability, doids)
        }

        private fun readSomaticEvents(gene: String, variant: String, record: CSVRecord): List<SomaticEvent> {
            val events = mutableListOf<SomaticEvent>()
            if (hasVariant(record)) {
                events.add(KnowledgebaseVariant(record["gene"], record["chromosome"], record["start"].toLong(), record["reference_bases"],
                                                record["variant_bases"]))
            }
            if (hasHgvs(record)) {
                val hgvsParts = extractHgvs(record).split(":")
                events.add(CDnaAnnotation(hgvsParts[0], hgvsParts[1]))
            }
            return events + listOfNotNull(readCNV(gene, variant), readFusion(gene, variant, record))
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

        private fun readCNV(gene: String, variant: String): CnvEvent? = when (variant) {
            "AMPLIFICATION" -> CnvEvent(gene, "Amplification")
            "DELETION"      -> CnvEvent(gene, "Deletion")
            else            -> null
        }

        private fun readFusion(gene: String, variant: String, record: CSVRecord): FusionEvent? {
            val variantTypes = record["variant_types"].orEmpty()
            return when {
                variantTypes.contains("fusion") -> fusionReader.read(gene, variant)
                else                            -> null
            }
        }

        private fun correctRecord(gene: String, variant: String): Pair<String, String> = when {
            variant.contains(Regex("MLL-MLLT3")) && gene == "KMT2A" ->
                Pair(gene, variant.replace("MLL-MLLT3", "KMT2A-MLLT3"))
            else                                                    -> Pair(gene, variant)
        }
    }
}
