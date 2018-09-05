package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GDnaVariant
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiActionableVariantReaderTest : StringSpec() {
    companion object {
        private val actionableInput = CgiActionableInput("EGFR", "ENST00000275493", "EGFR:T790M", "MUT",
                                                         "chr7:g.55249071C>T", "c.2369C>T", "EGFR:T790M", "",
                                                         "", "", "", "")

        private val actionableInputMultipleAlterations = CgiActionableInput("EGFR", "ENST00000275493", "EGFR:T790M,V600E", "MUT",
                                                                            "chr7:g.55249071C>T", "c.2369C>T", "EGFR:T790M", "",
                                                                            "", "", "", "")
    }

    init {
        "can read all events from input" {
            CgiActionableVariantReader.read(actionableInput).toSet() shouldBe setOf(
                    ProteinAnnotation("ENST00000275493", "T790M", SequenceVariantType.SUBSTITUTION),
                    CDnaAnnotation("ENST00000275493", "c.2369C>T", SequenceVariantType.SUBSTITUTION),
                    GDnaVariant("chr7:g.55249071C>T"))
        }

        "can read all events from input with multiple alterations" {
            CgiActionableVariantReader.read(actionableInputMultipleAlterations).toSet() shouldBe setOf(
                    ProteinAnnotation("ENST00000275493", "T790M", SequenceVariantType.SUBSTITUTION),
                    ProteinAnnotation("ENST00000275493", "V600E", SequenceVariantType.SUBSTITUTION),
                    CDnaAnnotation("ENST00000275493", "c.2369C>T", SequenceVariantType.SUBSTITUTION),
                    GDnaVariant("chr7:g.55249071C>T"))
        }

        "does not read from FUS input" {
            CgiActionableVariantReader.read(actionableInput.copy(`Alteration type` = "FUS")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from CNA input" {
            CgiActionableVariantReader.read(actionableInput.copy(`Alteration type` = "CNA")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from EXPR input" {
            CgiActionableVariantReader.read(actionableInput.copy(`Alteration type` = "EXPR")) shouldBe emptyList<SomaticEvent>()
        }
        "does not read from BIA input" {
            CgiActionableVariantReader.read(actionableInput.copy(`Alteration type` = "BIA")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
