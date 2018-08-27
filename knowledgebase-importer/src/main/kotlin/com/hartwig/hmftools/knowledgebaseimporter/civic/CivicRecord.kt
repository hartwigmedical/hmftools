package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.Doid
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import org.apache.logging.log4j.LogManager

data class CivicRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                       override val events: List<SomaticEvent>, override val actionability: List<Actionability>,
                       val cancerDoids: Map<String, Doid>) :
        RecordMetadata by metadata, KnownRecord, ActionableRecord {
    companion object {
        private val logger = LogManager.getLogger("CivicRecord")
        private val FUSION_SEPARATORS = setOf("-")
        private val FUSIONS_TO_FILTER = setOf(FusionPair("BRAF", "CUL1"))
        private val fusionReader = FusionReader(separators = FUSION_SEPARATORS, filterSet = FUSIONS_TO_FILTER)

        operator fun invoke(input: CivicVariantInput, evidence: Collection<CivicEvidence>): CivicRecord {
            val metadata = CivicMetadata(input.gene, input.representative_transcript)
            val additionalInfo = additionalInfo(evidence)
            val actionability = evidence.filter { it.direction == "Supports" }.flatMap { it.actionabilityItems }
            val doids = evidence.associateBy({ it.disease }, { Doid(it.doid) })
            val somaticEvents = readSomaticEvents(input)
            if (actionability.isNotEmpty() && somaticEvents.isEmpty()) {
                val aOrBLevelCount = actionability.filter { it.hmfLevel == "A" || it.hmfLevel == "B" }.size
                logger.warn("Could not extract somatic event from:\tcivic\t${input.gene}\t${input.variant}\t${input.variantTypes}\t$aOrBLevelCount")
            }
            return CivicRecord(metadata, additionalInfo, somaticEvents, actionability, doids)
        }

        private fun readSomaticEvents(input: CivicVariantInput): List<SomaticEvent> {
            return when {
                input.isFusion               -> listOfNotNull(fusionReader.read(input.gene, input.variant))
                input.isCNV                  -> listOfNotNull(readCNV(input.gene, input.variant))
                input.isVariantRecord        -> readVariants(input)
                input.isGenomicRangeMutation -> listOf(readGenomicRange(input))
                else                         -> readOtherEvents(input)
            }
        }

        private fun readVariants(input: CivicVariantInput): List<SomaticEvent> {
            val events = mutableListOf<SomaticEvent>()
            if (input.hasKnownVariant) events.add(readKnowledgebaseVariant(input))
            if (input.hasHgvs) events.add(readCDnaAnnotation(input))
            return events
        }

        private fun readKnowledgebaseVariant(input: CivicVariantInput) =
                KnowledgebaseVariant(input.gene, input.chromosome, input.start.toLong(), input.reference_bases, input.variant_bases)

        private fun readCDnaAnnotation(input: CivicVariantInput): CDnaAnnotation {
            val hgvsParts = input.hgvs!!.split(":")
            return CDnaAnnotation(hgvsParts[0], hgvsParts[1])
        }

        private fun readCNV(gene: String, variant: String): CnvEvent? = when (variant) {
            "AMPLIFICATION" -> CnvEvent.amplification(gene)
            "DELETION"      -> CnvEvent.deletion(gene)
            else            -> null
        }

        private fun readGenomicRange(input: CivicVariantInput): GenericRangeMutations {
            return GenericRangeMutations(input.gene, input.representative_transcript.substringBefore("."), input.start.toInt(),
                                         input.stop.toInt())
        }

        private fun readOtherEvents(input: CivicVariantInput): List<SomaticEvent> {
            val events = mutableListOf<SomaticEvent>()
            if (input.hasVariant) {
                events.add(OtherEvents(readVariants(input)))
            }
            return events
        }

        private fun additionalInfo(evidence: Collection<CivicEvidence>): String {
            val highestEvidenceLevel = evidence.map { it.level }.sorted().firstOrNull() ?: "N"
            return (highestEvidenceLevel == "A" || highestEvidenceLevel == "B" || highestEvidenceLevel == "C").toString()
        }
    }
}
