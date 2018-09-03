package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CivicCDnaAnnotationReaderTest : StringSpec() {
    companion object {
        private const val hgvsStart = "ENST00000288602.6:c.1799T>A,NM_004333.4:c.1799T>A,NP_004324.2:p.Val600Glu,NC_000007.13:g.140453136A>T"
        private const val hgvsMiddle = "NM_004333.4:c.1799T>A,ENST00000288602.6:c.1799T>A,NP_004324.2:p.Val600Glu,NC_000007.13:g.140453136A>T"
        private const val hgvsEnd = "NM_004333.4:c.1799T>A,NP_004324.2:p.Val600Glu,NC_000007.13:g.140453136A>T,ENST00000288602.6:c.1799T>A"
        private const val noEnsemblHgvs = "NM_004333.4:c.1799T>A,NP_004324.2:p.Val600Glu,NC_000007.13:g.140453136A>T"
        val variant = CivicVariantInput("BRAF", "ENST00000288602", "", "V600E", "7",
                                        "140453136", "140453136", "A", "T", hgvsEnd, "")
        private val expectedResult = CDnaAnnotation("ENST00000288602.6", "c.1799T>A", SequenceVariantType.SUBSTITUTION)
    }

    init {
        "can read cdna annotation" {
            CivicCDnaAnnotationReader.read(variant) shouldBe listOf(expectedResult)
        }

        "reads transcript from hgvs annotation field" {
            CivicCDnaAnnotationReader.read(variant.copy(representative_transcript = "ENST00000")) shouldBe listOf(expectedResult)
        }

        "reads cdna annotation regardless of position" {
            CivicCDnaAnnotationReader.read(variant.copy(hgvs_expressions = hgvsStart)) shouldBe listOf(expectedResult)
            CivicCDnaAnnotationReader.read(variant.copy(hgvs_expressions = hgvsMiddle)) shouldBe listOf(expectedResult)
        }

        "only reads ensembl annotation" {
            CivicCDnaAnnotationReader.read(variant.copy(hgvs_expressions = noEnsemblHgvs)) shouldBe emptyList<SomaticEvent>()
        }
    }
}
