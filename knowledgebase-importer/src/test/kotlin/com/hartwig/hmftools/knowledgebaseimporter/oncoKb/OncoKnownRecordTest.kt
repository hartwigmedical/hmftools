package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.HgvsVariantType
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKnownInput
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class OncoKnownRecordTest : StringSpec() {
    private val known = OncoKnownInput("ENST00000000", "KIT", "", "Loss-of-function", "")

    init {
        "filters out fusion events with loss-of-function effect" {
            OncoKnownRecord(known.copy(Alteration = "KIT Fusions")).events shouldBe emptyList<SomaticEvent>()
            OncoKnownRecord(known.copy(Alteration = "KIT-PIK3 Fusion")).events shouldBe emptyList<SomaticEvent>()
        }

        "does not filter out protein alterations with loss-of-function effect" {
            OncoKnownRecord(known.copy(Alteration = "V600E")).events shouldBe
                    listOf(ProteinAnnotation(known.transcript!!, "V600E", HgvsVariantType.SUBSTITUTION))
        }
    }
}
