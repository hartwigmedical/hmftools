package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiActionableMutationListReaderTest : StringSpec() {
    companion object {
        private val actionableInput = CgiActionableInput("ESR1", "ENST00000", "ESR1:.,E380Q,.537.,.538.,L536.,P535H,500-600",
                                                         "MUT", "", "", "", "",
                                                         "", "", "", "")

    }

    init {
        "can read multiple mutations form alteration field" {
            CgiActionableMutationListReader.read(actionableInput) shouldBe listOf(
                    GeneMutations("ESR1", "ENST00000"),
                    ProteinAnnotation("ENST00000", "E380Q", SequenceVariantType.SUBSTITUTION),
                    CodonMutations("ESR1", "ENST00000", 537),
                    CodonMutations("ESR1", "ENST00000", 538),
                    CodonMutations("ESR1", "ENST00000", 536),
                    ProteinAnnotation("ENST00000", "P535H", SequenceVariantType.SUBSTITUTION),
                    CodonRangeMutations("ESR1", "ENST00000", 500, 600))
        }

        "does not read from individual mutation field" {
            CgiActionableMutationListReader.read(actionableInput.copy(Alteration = "", individual_mutation = "ESR1:E380Q")) shouldBe
                    emptyList<SomaticEvent>()
        }

        "does not read from cdna field" {
            CgiActionableMutationListReader.read(actionableInput.copy(Alteration = "", cDNA = "c.944C>T")) shouldBe
                    emptyList<SomaticEvent>()
        }

        "does not read from gdna field" {
            CgiActionableMutationListReader.read(actionableInput.copy(Alteration = "", gDNA = "chr9:g.133748283C>T")) shouldBe
                    emptyList<SomaticEvent>()
        }

        "does not read from FUS input" {
            CgiActionableMutationListReader.read(actionableInput.copy(`Alteration type` = "FUS")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from CNA input" {
            CgiActionableMutationListReader.read(actionableInput.copy(`Alteration type` = "CNA")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from EXPR input" {
            CgiActionableMutationListReader.read(actionableInput.copy(`Alteration type` = "EXPR")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from BIA input" {
            CgiActionableMutationListReader.read(actionableInput.copy(`Alteration type` = "BIA")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
