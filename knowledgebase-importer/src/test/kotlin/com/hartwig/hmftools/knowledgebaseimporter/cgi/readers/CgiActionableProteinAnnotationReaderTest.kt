package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiActionableProteinAnnotationReaderTest : StringSpec() {
    companion object {
        private val actionableInput = CgiActionableInput("ABL1", "ENST00000318560", "ABL1:T315I", "MUT",
                                                         "chr9:g.133748283C>T", "c.944C>T", "ABL1:T315I", "", "", "", "", "")
    }

    init {
        "can read protein annotation from actionable input" {
            CgiActionableProteinAnnotationReader.read(actionableInput) shouldBe
                    listOf(ProteinAnnotation("ENST00000318560", "T315I", SequenceVariantType.SUBSTITUTION))
        }

        "can read protein annotation from actionable with empty alteration input" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(Alteration = "")) shouldBe
                    listOf(ProteinAnnotation("ENST00000318560", "T315I", SequenceVariantType.SUBSTITUTION))
        }

        "does not read from FUS input" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(`Alteration type` = "FUS")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from CNA input" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(`Alteration type` = "CNA")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from EXPR input" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(`Alteration type` = "EXPR")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from BIA input" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(`Alteration type` = "BIA")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with missing individual mutation" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(individual_mutation = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with codon mutation" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(individual_mutation = "ABL1:T315.")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with any mutation" {
            CgiActionableProteinAnnotationReader.read(actionableInput.copy(individual_mutation = "ABL1:.")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
