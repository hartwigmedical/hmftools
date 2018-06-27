package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.Multimap
import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager
import java.util.regex.Pattern

data class CivicRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                       override val events: List<SomaticEvent>, override val actionability: List<Actionability>,
                       val cancerDoids: Map<String, String>) :
        RecordMetadata by metadata, KnownRecord, ActionableRecord {
    companion object {
        private val logger = LogManager.getLogger("CivicRecord")
        private const val RANGE_VARIANTS = "gene_variant|transcript_variant|exon_variant|coding_sequence_variant|protein_altering_variant"
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
            val somaticEvents = readSomaticEvents(gene, variant, record)
            if (actionability.isNotEmpty() && somaticEvents.isEmpty())
                logger.warn("Could not extract any somatic event from ${record["gene"]}: ${record["variant"]} (${variantTypes(record)})")
            return CivicRecord(metadata, additionalInfo, somaticEvents, actionability, doids)
        }

        private fun readSomaticEvents(gene: String, variant: String, record: CSVRecord): List<SomaticEvent> {
            return when {
                isFusion(record)               -> listOfNotNull(fusionReader.read(gene, variant))
                isCNV(variant)                 -> listOfNotNull(readCNV(gene, variant))
                isVariantRecord(record)        -> readVariants(record)
                isGenomicRangeMutation(record) -> listOf(readGenomicRange(record))
                else                           -> emptyList()
            }
        }

        private fun isCNV(variant: String) = variant == "AMPLIFICATION" || variant == "DELETION"
        private fun isFusion(record: CSVRecord) = variantTypes(record).isNotEmpty() && variantTypes(record).all { it.contains("fusion") }
        private fun isVariantRecord(record: CSVRecord) = variantTypes(record).none { it.contains("fusion") } &&
                (hasKnownVariant(record) || hasHgvs(record))

        private fun isGenomicRangeMutation(record: CSVRecord) = hasPosition(record) && !hasRefOrAlt(record) &&
                variantTypes(record).isNotEmpty() && (isGenericMutation(record) || isGenericMissense(record))

        private fun isGenericMutation(record: CSVRecord) = record["variant"].toLowerCase() == "mutation" ||
                variantTypes(record).any { it.contains(RANGE_VARIANTS.toRegex()) }

        private fun isGenericMissense(record: CSVRecord) = !record["variant"].contains("+") && variantTypes(record).size == 1 &&
                !record["variant"].toLowerCase().contains(" and ") && variantTypes(record).first() == "missense_variant"

        private fun hasPosition(record: CSVRecord) = !record["chromosome"].isNullOrBlank()
                && !record["start"].isNullOrBlank() && !record["stop"].isNullOrBlank()

        private fun hasRefOrAlt(record: CSVRecord) = !record["reference_bases"].isNullOrBlank() || !record["variant_bases"].isNullOrBlank()

        private fun hasKnownVariant(record: CSVRecord): Boolean = hasPosition(record) && hasRefOrAlt(record)

        private fun hasHgvs(record: CSVRecord) = extractHgvs(record).isNotBlank()

        private fun readVariants(record: CSVRecord): List<SomaticEvent> {
            val events = mutableListOf<SomaticEvent>()
            if (hasKnownVariant(record)) {
                events.add(KnowledgebaseVariant(record["gene"], record["chromosome"], record["start"].toLong(), record["reference_bases"],
                                                record["variant_bases"]))
            }
            if (hasHgvs(record)) {
                val hgvsParts = extractHgvs(record).split(":")
                events.add(CDnaAnnotation(hgvsParts[0], hgvsParts[1]))
            }
            return events
        }

        private fun extractHgvs(record: CSVRecord): String {
            val hgvsField = record["hgvs_expressions"].orEmpty()
            val matcher = Pattern.compile("(ENST[0-9]+\\.[0-9+]:c\\.[0-9][^,\\t\\s\\n]+)").matcher(hgvsField)
            return if (matcher.find()) matcher.group(1) else ""
        }

        private fun readCNV(gene: String, variant: String): CnvEvent? = when (variant) {
            "AMPLIFICATION" -> CnvEvent(gene, "Amplification")
            "DELETION"      -> CnvEvent(gene, "Deletion")
            else            -> null
        }

        private fun readGenomicRange(record: CSVRecord): GenericRangeMutations {
            return GenericRangeMutations(record["gene"], record["representative_transcript"].orEmpty().substringBefore("."),
                                         record["start"].toInt(), record["stop"].toInt())
        }

        private fun variantTypes(record: CSVRecord): List<String> {
            val variantTypes = record["variant_types"].orEmpty()
            return if (variantTypes.matches("N/A".toRegex(RegexOption.IGNORE_CASE))) {
                emptyList()
            } else {
                variantTypes.split(",")
            }
        }

        private fun additionalInfo(evidence: Collection<CivicEvidence>): String {
            val highestEvidenceLevel = evidence.map { it.level }.sorted().firstOrNull() ?: "N"
            return (highestEvidenceLevel == "A" || highestEvidenceLevel == "B" || highestEvidenceLevel == "C").toString()
        }

        private fun correctRecord(gene: String, variant: String): Pair<String, String> = when {
            variant.contains(Regex("MLL-MLLT3")) && gene == "KMT2A" ->
                Pair(gene, variant.replace("MLL-MLLT3", "KMT2A-MLLT3"))
            else                                                    -> Pair(gene, variant)
        }
    }
}
