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
        "can read range mutation from exon variant record" {
            CivicRangeMutationReader.read(exonMutation) shouldBe
                    listOf(GenericRangeMutations("PIK3CA", "ENST00000263967", 178951882, 178952495))
        }

        "can read range mutation from codon variant record" {
            CivicRangeMutationReader.read(codonMutation) shouldBe
                    listOf(GenericRangeMutations("DNMT3A", "ENST00000264709", 25457241, 25457243))
        }

        "can read range mutation from domain mutations record" {
            CivicRangeMutationReader.read(codonMutation.copy(variant = "B2 DOMAIN MUTATION")) shouldBe
                    listOf(GenericRangeMutations("DNMT3A", "ENST00000264709", 25457241, 25457243))
        }

        "can read range mutation from g12/13 record" {
            CivicRangeMutationReader.read(codonMutation.copy(variant = "G12/G13")) shouldBe
                    listOf(GenericRangeMutations("DNMT3A", "ENST00000264709", 25457241, 25457243))
        }

        "can read range mutation from mutation record" {
            CivicRangeMutationReader.read(codonMutation.copy(variant = "MUTATION")) shouldBe
                    listOf(GenericRangeMutations("DNMT3A", "ENST00000264709", 25457241, 25457243))
        }

        "does not read range mutation from record without chromosome" {
            CivicRangeMutationReader.read(codonMutation.copy(chromosome = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read range mutation from record without start" {
            CivicRangeMutationReader.read(codonMutation.copy(start = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read range mutation from record without stop" {
            CivicRangeMutationReader.read(codonMutation.copy(stop = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read range mutation from record with position and ref/alt" {
            CivicRangeMutationReader.read(codonMutation.copy(reference_bases = "A")) shouldBe emptyList<SomaticEvent>()
            CivicRangeMutationReader.read(codonMutation.copy(variant_bases = "A")) shouldBe emptyList<SomaticEvent>()

        }
    }
}
