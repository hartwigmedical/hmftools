package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.PromiscuousGene
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CivicFusionReaderTest : StringSpec() {
    private val complexVariant = CivicVariantInput("ALK", "ENST00000389048.3", "", "ALK FUSION F1245C", "2",
                                                   "29436859", "29436859", "A", "C",
                                                   "NM_004304.4:c.3734T>G,NP_004295.2:p.Phe1245Cys,NC_000002.11:g.29436859A>C,ENST00000389048.3:c.3734T>G",
                                                   "missense_variant,transcript_fusion")
    private val fusionVariant = CivicVariantInput("NTRK1", "ENST00000368300.4", "", "LMNA-NTRK1", "1",
                                                  "156084498", "156108548", "", "", "", "transcript_fusion")

    private val promiscuousGene = CivicVariantInput("FGFR1", "ENST00000425967.3", "", "FGFR1 FUSIONS", "8",
                                                    "38268656", "38325363", "", "", "", "transcript_fusion")

    init {
        "can read fusion pair" {
            CivicFusionReader.read(fusionVariant) shouldBe listOf(FusionPair("LMNA", "NTRK1"))
        }

        "can read promiscuous gene" {
            CivicFusionReader.read(promiscuousGene) shouldBe listOf(PromiscuousGene("FGFR1"))
        }

        "does not read fusions from input without variant types" {
            CivicFusionReader.read(fusionVariant.copy(variant_types = "")) shouldBe emptyList<SomaticEvent>()
            CivicFusionReader.read(promiscuousGene.copy(variant_types = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not extract fusion from complex variant" {
            CivicFusionReader.read(complexVariant) shouldBe emptyList<SomaticEvent>()
        }
    }
}
