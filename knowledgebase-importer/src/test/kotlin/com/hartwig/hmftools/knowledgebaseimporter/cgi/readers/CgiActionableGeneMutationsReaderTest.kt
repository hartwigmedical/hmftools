package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GeneMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiActionableGeneMutationsReaderTest : StringSpec() {

    companion object {
        private val actionableInput = CgiActionableInput("ARAF", null, "ARAF:.", "MUT", "", "",
                                                         "", "", "", "", "", "", ".")
    }

    init {
        "can read gene mutations from actionable input" {
            CgiActionableGeneMutationsReader.read(actionableInput) shouldBe listOf(GeneMutations("ARAF", null))
        }

        "does not read from FUS input" {
            CgiActionableGeneMutationsReader.read(actionableInput.copy(`Alteration type` = "FUS")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from CNA input" {
            CgiActionableGeneMutationsReader.read(actionableInput.copy(`Alteration type` = "CNA")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from EXPR input" {
            CgiActionableGeneMutationsReader.read(actionableInput.copy(`Alteration type` = "EXPR")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from BIA input" {
            CgiActionableGeneMutationsReader.read(actionableInput.copy(`Alteration type` = "BIA")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with missing mutation" {
            CgiActionableGeneMutationsReader.read(actionableInput.copy(variant = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with codon mutation" {
            CgiActionableGeneMutationsReader.read(actionableInput.copy(variant = "V600.")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with codon range mutation" {
            CgiActionableGeneMutationsReader.read(actionableInput.copy(variant = "550-600")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
