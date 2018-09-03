package com.hartwig.hmftools.knowledgebaseimporter.iclusion

import com.hartwig.hmftools.apiclients.iclusion.data.IclusionStudyDetails
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.Doid
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers.*
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnownRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse

data class IclusionRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                          override val actionability: List<Actionability>, override val additionalInfo: String,
                          val iclusionEvents: List<IclusionEvent>, val doids: Map<String, Set<Doid>>) :
        RecordMetadata by metadata, KnownRecord, ActionableRecord {

    companion object {
        private val reader = KnowledgebaseEventReader("iclusion", IclusionTransvarReader, IclusionFusionReader, IclusionCnvReader,
                                                      IclusionGeneMutationReader, IclusionExonMutationReader, IclusionCodonReader)

        operator fun invoke(studyDetails: IclusionStudyDetails, geneTranscript: Map<String, String?>): List<IclusionRecord> {
            val events = studyDetails.mutations.map { IclusionEvent(it, geneTranscript[it.geneName].orEmpty()) }
            val actionability = readActionability(studyDetails)
            val doids = studyDetails.indications.map { Pair(it.indication_name_full, it.doidSet.map { Doid(it) }.toSet()) }
                    .toMap().filterValues { it.isNotEmpty() }
            // MIVO: for now, interpret each iclusion study mutation as separate record. Effectively treats the mutations as an OR predicate
            //      e.g. patient will match if ANY of the specified mutations match
            return events.map {
                IclusionRecord(IclusionMetadata(it.gene, it.transcript), reader.read(it), actionability, "", listOf(it), doids)
            }
        }

        private fun readActionability(studyDetails: IclusionStudyDetails): List<Actionability> {
            val study = studyDetails.study
            val cancerTypes = studyDetails.indications.map { it.indication_name_full }
            val drugs = listOf(HmfDrug("Study", "Study"))
            return Actionability("iclusion", study.title, cancerTypes, drugs, "study_level", "study_significance",
                                 "Predictive", HmfLevel.UNKNOWN, HmfResponse.OTHER)
        }
    }
}
