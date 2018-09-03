package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.civic.readers.*
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.Doid
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnownRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability

data class CivicRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                       override val events: List<SomaticEvent>, override val actionability: List<Actionability>,
                       val cancerDoids: Map<String, Doid>) : RecordMetadata by metadata, KnownRecord, ActionableRecord {
    companion object {
        private val reader = KnowledgebaseEventReader("civic", CivicVariantReader, CivicCnvReader, CivicFusionReader,
                                                      CivicRangeMutationReader, CivicMultipleEventsReader)

        operator fun invoke(input: CivicVariantInput, evidence: Collection<CivicEvidence>): CivicRecord {
            val metadata = CivicMetadata(input.gene, input.transcript)
            val additionalInfo = additionalInfo(evidence)
            val actionability = evidence.filter { it.direction == "Supports" }.flatMap { it.actionabilityItems }
            val doids = evidence.associateBy({ it.disease }, { Doid(it.doid) })
            val somaticEvents = reader.read(input)
            return CivicRecord(metadata, additionalInfo, somaticEvents, actionability, doids)
        }

        private fun additionalInfo(evidence: Collection<CivicEvidence>): String {
            val highestEvidenceLevel = evidence.map { it.level }.sorted().firstOrNull() ?: "N"
            return (highestEvidenceLevel == "A" || highestEvidenceLevel == "B" || highestEvidenceLevel == "C").toString()
        }
    }
}
