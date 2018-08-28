package com.hartwig.hmftools.knowledgebaseimporter.iclusion

import com.hartwig.hmftools.apiclients.iclusion.data.IclusionMutationDetails
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent

data class IclusionEvent(override val gene: String, override val variant: String, override val transcript: String) : KnowledgebaseEvent {
    companion object {
        operator fun invoke(mutation: IclusionMutationDetails, transcript: String): IclusionEvent {
            return IclusionEvent(mutation.geneName, mutation.variantName, transcript)
        }
    }
}
