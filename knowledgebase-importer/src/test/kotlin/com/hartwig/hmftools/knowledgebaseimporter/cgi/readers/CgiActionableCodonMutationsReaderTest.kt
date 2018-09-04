package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiActionableCodonMutationsReaderTest : StringSpec() {

    companion object {
        private val actionableInput = CgiActionableInput("BRAF", null, "BRAF:V600.", "MUT", "", "",
                                                         "", "", "", "", "", "", "V600.")
    }

    init {
        "can read codon mutation from actionable input" {
            CgiActionableCodonMutationsReader.read(actionableInput) shouldBe listOf(CodonMutations("BRAF", null, 600))
        }

        "can read codon mutation from input with generic codon" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(variant = ".537.")) shouldBe
                    listOf(CodonMutations("BRAF", null, 537))
        }

        "does not read from FUS input" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(`Alteration type` = "FUS")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from CNA input" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(`Alteration type` = "CNA")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from EXPR input" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(`Alteration type` = "EXPR")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from BIA input" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(`Alteration type` = "BIA")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with missing mutation" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(variant = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with codon range" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(variant = "500-600")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with any mutation" {
            CgiActionableCodonMutationsReader.read(actionableInput.copy(variant = ".")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
