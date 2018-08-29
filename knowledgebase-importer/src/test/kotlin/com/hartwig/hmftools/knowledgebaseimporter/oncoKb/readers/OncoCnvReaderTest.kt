package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class OncoCnvReaderTest : StringSpec() {
    private val actionable = OncoActionableInput("ENST00000000", "BRAF", "", "", "", "")

    init {
        "can read amplification" {
            OncoCnvReader.read(actionable.copy(Alteration = "Amplification")) shouldBe listOf(CnvEvent.amplification(actionable.gene))
        }

        "can read deletion" {
            OncoCnvReader.read(actionable.copy(Alteration = "Deletion")) shouldBe listOf(CnvEvent.deletion(actionable.gene))
        }

        "will only read exact matches" {
            OncoCnvReader.read(actionable.copy(Alteration = "amplification")) shouldBe emptyList<CnvEvent>()
            OncoCnvReader.read(actionable.copy(Alteration = "amp")) shouldBe emptyList<CnvEvent>()
            OncoCnvReader.read(actionable.copy(Alteration = "deletion")) shouldBe emptyList<CnvEvent>()
            OncoCnvReader.read(actionable.copy(Alteration = "del")) shouldBe emptyList<CnvEvent>()
        }
    }
}