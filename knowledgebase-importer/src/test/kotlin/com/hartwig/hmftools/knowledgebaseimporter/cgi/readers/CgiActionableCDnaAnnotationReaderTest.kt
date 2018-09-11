package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiActionableCDnaAnnotationReaderTest : StringSpec() {
    companion object {
        private val actionableInput = CgiActionableInput("ABL1", "ENST00000318560", "ABL1:T315I", "MUT",
                                                         "chr9:g.133748283C>T", "c.944C>T", "ABL1:T315I", "", "", "", "", "")
    }

    init {
        "can read cdna annotation from actionable input" {
            CgiActionableCDnaAnnotationReader.read(actionableInput) shouldBe
                    listOf(CDnaAnnotation("ENST00000318560", "c.944C>T", SequenceVariantType.SUBSTITUTION))
        }

        "does not read from FUS input" {
            CgiActionableCDnaAnnotationReader.read(actionableInput.copy(`Alteration type` = "FUS")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from CNA input" {
            CgiActionableCDnaAnnotationReader.read(actionableInput.copy(`Alteration type` = "CNA")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from EXPR input" {
            CgiActionableCDnaAnnotationReader.read(actionableInput.copy(`Alteration type` = "EXPR")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from BIA input" {
            CgiActionableCDnaAnnotationReader.read(actionableInput.copy(`Alteration type` = "BIA")) shouldBe emptyList<SomaticEvent>()
        }

        "produces empty list for input with missing cdna" {
            CgiActionableCDnaAnnotationReader.read(actionableInput.copy(cDNA = "")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
