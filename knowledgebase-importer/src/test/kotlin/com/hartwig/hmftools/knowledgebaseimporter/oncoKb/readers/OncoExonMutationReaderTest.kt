package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ExonMutations
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKnownInput
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class OncoExonMutationReaderTest : StringSpec() {
    private val actionable = OncoActionableInput("ENST00000000", "BRAF", "", "", "", "")
    private val known = OncoKnownInput("ENST00000000", "KIT", "", "", "")

    init {
        "matches actionable exon mutations" {
            OncoExonMutationReader.read(actionable.copy(Alteration = "Exon 13 mutations")) shouldBe
                    listOf(ExonMutations(actionable.gene, actionable.transcript, 13))
        }

        "matches known exon mutations" {
            OncoExonMutationReader.read(known.copy(Alteration = "Exon 9 mutations")) shouldBe
                    listOf(ExonMutations(known.gene, known.transcript, 9))
        }

        "match ignores case" {
            OncoExonMutationReader.read(actionable.copy(Alteration = "exon 13 Mutations")) shouldBe
                    listOf(ExonMutations(actionable.gene, actionable.transcript, 13))
        }

        "does not match exon deletions" {
            OncoExonMutationReader.read(actionable.copy(Alteration = "Exon 13 deletions")) shouldBe emptyList<ExonMutations>()
            OncoExonMutationReader.read(known.copy(Alteration = "Exon 9 deletions")) shouldBe emptyList<ExonMutations>()
        }

        "does not match exon insertions" {
            OncoExonMutationReader.read(actionable.copy(Alteration = "Exon 13 insertions")) shouldBe emptyList<ExonMutations>()
            OncoExonMutationReader.read(known.copy(Alteration = "Exon 9 insertions")) shouldBe emptyList<ExonMutations>()
        }
    }
}
