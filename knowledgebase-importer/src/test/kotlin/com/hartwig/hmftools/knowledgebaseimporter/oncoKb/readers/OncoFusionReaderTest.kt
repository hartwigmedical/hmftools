package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.PromiscuousGene
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class OncoFusionReaderTest : StringSpec() {
    private val actionable = OncoActionableInput("ENST00000000", "BRAF", "", "", "", "")

    init {
        "can read actionable fusion pair" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF-PIK3 Fusion")) shouldBe
                    listOf(FusionPair(actionable.gene, "PIK3"))
        }

        "can read promiscuous gene" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF Fusions")) shouldBe
                    listOf(PromiscuousGene(actionable.gene))
        }

        "can read fusion with ? separator" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF?PIK3 Fusion")) shouldBe
                    listOf(FusionPair(actionable.gene, "PIK3"))
        }

        "can read fusion with ' - ' separator" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF - PIK3 Fusion")) shouldBe
                    listOf(FusionPair(actionable.gene, "PIK3"))
        }

        "does not match fusions without 'Fusion' keyword" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF-PIK3")) shouldBe emptyList<FusionEvent>()
        }
    }
}
