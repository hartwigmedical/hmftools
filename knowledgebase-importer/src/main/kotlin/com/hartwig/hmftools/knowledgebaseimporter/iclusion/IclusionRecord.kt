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
        private val reader = KnowledgebaseEventReader(IclusionTransvarReader, IclusionFusionReader, IclusionCnvReader,
                                                      IclusionGeneMutationReader, IclusionExonMutationReader)

        operator fun invoke(studyDetails: IclusionStudyDetails): IclusionRecord {
            val events = studyDetails.mutations.map { IclusionEvent(it) }
            val somaticEvents = events.flatMap { reader.read(it) }
            return IclusionRecord(IclusionMetadata("todo", "wip"), somaticEvents, emptyList(), events)
        }
    }
}
