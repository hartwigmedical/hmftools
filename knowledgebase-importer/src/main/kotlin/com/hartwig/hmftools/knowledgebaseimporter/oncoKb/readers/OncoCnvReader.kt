package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKbInput
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent

object OncoCnvReader : SomaticEventReader<OncoKbInput, CnvEvent> {
    override fun read(event: OncoKbInput): List<CnvEvent> {
        return when (event.variant) {
            "Amplification" -> listOf(CnvEvent.amplification(event.gene))
            "Deletion"      -> listOf(CnvEvent.deletion(event.gene))
            else            -> emptyList()
        }
    }
}
