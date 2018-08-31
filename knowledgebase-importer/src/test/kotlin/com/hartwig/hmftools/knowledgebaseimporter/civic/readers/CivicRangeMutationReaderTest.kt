package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GenericRangeMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CivicRangeMutationReaderTest : StringSpec() {
    private val exonMutation = CivicVariantInput("PIK3CA", "ENST00000263967.3", "", "EXON 21 MUTATION", "3",
                                                 "178951882", "178952495", "", "", "", "exon_variant")
    private val codonMutation = CivicVariantInput("DNMT3A", "ENST00000264709.3", "", "R882", "2",
                                                  "25457241", "25457243", "", "", "", "missense_variant")

    init {
        "can read range mutation from exon_variant record" {
            CivicRangeMutationReader.read(exonMutation) shouldBe
                    listOf(GenericRangeMutations("PIK3CA", "ENST00000263967", 178951882, 178952495))
        }

        "can read range mutation from generic missense record" {
            CivicRangeMutationReader.read(codonMutation) shouldBe
                    listOf(GenericRangeMutations("DNMT3A", "ENST00000264709", 25457241, 25457243))
        }

        "does not read range mutation from record without variant types" {
            CivicRangeMutationReader.read(codonMutation.copy(variant_types = "")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
