package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonRangeMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiActionableCodonRangeMutationsReaderTest : StringSpec() {

    companion object {
        private val actionableInput = CgiActionableInput("KIT", null, "KIT:550-592", "MUT", "", "",
                                                         "", "", "", "", "", "", "550-592")
    }

    init {
        "can read codon range mutation from actionable input" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput) shouldBe
                    listOf(CodonRangeMutations("KIT", null, 550, 592))
        }

        "does not read from FUS input" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput.copy(`Alteration type` = "FUS")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from CNA input" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput.copy(`Alteration type` = "CNA")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from EXPR input" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput.copy(`Alteration type` = "EXPR")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from BIA input" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput.copy(`Alteration type` = "BIA")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with missing mutation" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput.copy(variant = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with codon mutation" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput.copy(variant = "V600.")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with any mutation" {
            CgiActionableCodonRangeMutationsReader.read(actionableInput.copy(variant = ".")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
