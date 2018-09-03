package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent

object CivicCnvReader : SomaticEventReader<CivicVariantInput, CnvEvent> {
    override fun read(event: CivicVariantInput): List<CnvEvent> {
        return when (event.variant) {
            "AMPLIFICATION" -> listOf(CnvEvent.amplification(event.gene))
            "DELETION"      -> listOf(CnvEvent.deletion(event.gene))
            else            -> emptyList()
        }
    }
}
