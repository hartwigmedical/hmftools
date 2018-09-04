package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiKnownInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GDnaVariant
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CgiGDnaReaderTest : StringSpec() {
    companion object {
        private val knownInput = CgiKnownInput("ABL1", "ENST00000318560", "", "chr9:g.133738306G>A", "p.E236K")
        private val actionableInput = CgiActionableInput("ABL1", "ENST00000318560", "ABL1:T315I", "MUT",
                                                         "chr9:g.133748283C>T", "c.944C>T", "ABL1:T315I", "", "", "", "", "")
    }

    init {
        "can read gdna from known input" {
            CgiGDnaReader.read(knownInput) shouldBe listOf(GDnaVariant("chr9:g.133738306G>A"))
        }

        "can read gdna from actionable input" {
            CgiGDnaReader.read(actionableInput) shouldBe listOf(GDnaVariant("chr9:g.133748283C>T"))
        }
    }
}
