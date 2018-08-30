package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKnownInput
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.PromiscuousGene
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class OncoFusionReaderTest : StringSpec() {
    private val actionable = OncoActionableInput("ENST00000000", "BRAF", "", "", "", "")
    private val known = OncoKnownInput("ENST00000000", "KIT", "", "Loss-of-function", "")

    init {
        "can read fusion pair" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF-PIK3 Fusion")) shouldBe listOf(FusionPair(actionable.gene, "PIK3"))
            OncoFusionReader.read(known.copy(Alteration = "KIT-PIK3 Fusion")) shouldBe listOf(FusionPair(known.gene, "PIK3"))
        }

        "can read promiscuous gene" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF Fusions")) shouldBe listOf(PromiscuousGene(actionable.gene))
            OncoFusionReader.read(known.copy(Alteration = "KIT Fusions")) shouldBe listOf(PromiscuousGene(known.gene))
        }

        "can read fusion with ? separator" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF?PIK3 Fusion")) shouldBe listOf(FusionPair(actionable.gene, "PIK3"))
            OncoFusionReader.read(known.copy(Alteration = "KIT?PIK3 Fusion")) shouldBe listOf(FusionPair(known.gene, "PIK3"))
        }

        "can read fusion with ' - ' separator" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF - PIK3 Fusion")) shouldBe listOf(FusionPair(actionable.gene, "PIK3"))
            OncoFusionReader.read(known.copy(Alteration = "KIT - PIK3 Fusion")) shouldBe listOf(FusionPair(known.gene, "PIK3"))
        }

        "does not match fusions without 'Fusion' keyword" {
            OncoFusionReader.read(actionable.copy(Alteration = "BRAF-PIK3")) shouldBe emptyList<FusionEvent>()
        }

        "reads known fusion when effect is loss-of-function" {
            OncoFusionReader.read(known.copy(Alteration = "KIT-PIK3 Fusion")) shouldBe listOf(FusionPair(known.gene, "PIK3"))
        }
    }
}
