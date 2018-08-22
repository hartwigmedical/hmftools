package com.hartwig.hmftools.knowledgebaseimporter.iclusion

import com.hartwig.hmftools.apiclients.iclusion.data.IclusionStudyDetails
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers.*
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability

data class IclusionRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                          override val actionability: List<Actionability>, val iclusionEvents: List<IclusionEvent>) :
        RecordMetadata by metadata,
        ActionableRecord {
    companion object {
        private val reader = KnowledgebaseEventReader("iclusion", IclusionTransvarReader, IclusionFusionReader, IclusionCnvReader,
                                                      IclusionGeneMutationReader, IclusionExonMutationReader)

        operator fun invoke(studyDetails: IclusionStudyDetails, geneTranscript: Map<String, String?>): List<IclusionRecord> {
            val events = studyDetails.mutations.map { IclusionEvent(it, geneTranscript[it.geneName].orEmpty()) }
            // MIVO: for now, interpret each iclusion study mutation as separate record. Effectively treats the mutations as an OR predicate
            //      e.g. patient will match if ANY of the specified mutations match
            return events.map { IclusionRecord(IclusionMetadata(it.gene, it.transcript), reader.read(it), emptyList(), listOf(it)) }
        }
    }
}
