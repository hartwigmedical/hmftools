package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseVariant
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.OtherEvents
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CivicMultipleEventsReaderTest : StringSpec() {
    private val v600eV600m = CivicVariantInput("BRAF", "ENST00000288602.6", "", "V600E+V600M", "7",
                                               "140453135", "140453137", "", "", "", "missense_variant")
    private val fusionAndVariant = CivicVariantInput("ALK", "ENST00000389048.3", "", "EML4-ALK L1196M", "2",
                                                     "29443631", "29443631", "G", "T",
                                                     "ENST00000389048.3:c.3586C>A,NC_000002.11:g.29443631G>T,NM_004304.4:c.3586C>A,NP_004295.2:p.Leu1196Met",
                                                     "missense_variant,transcript_fusion")
    private val v600eAmplification = CivicVariantInput("BRAF", "ENST00000288602.6", "", "V600E AMPLIFICATION", "7",
                                                       "140434279", "140624564", "", "", "",
                                                       "missense_variant,transcript_amplification")
    private val fusionAnd2Variants = CivicVariantInput("NTRK1", "ENST00000524377.1", "", "LMNA-NTRK1 G595R AND G667C", "1",
                                                       "156846342", "156849107", "", "", "",
                                                       "missense_variant")

    private val alkPair = FusionPair("EML4", "ALK")

    init {
        "can read complex variant with fusion pair + variant" {
            CivicMultipleEventsReader.read(fusionAndVariant) shouldBe listOf(
                    OtherEvents(listOf(alkPair, KnowledgebaseVariant("ALK", "2", 29443631, "G", "T"))),
                    OtherEvents(listOf(alkPair, CDnaAnnotation("ENST00000389048.3", "c.3586C>A", SequenceVariantType.SUBSTITUTION))))
        }

        "can't read complex variant with 2 variants" {
            CivicMultipleEventsReader.read(v600eV600m) shouldBe emptyList<SomaticEvent>()
        }

        "can't read complex variant with variant + cnv" {
            CivicMultipleEventsReader.read(v600eAmplification) shouldBe emptyList<SomaticEvent>()
        }

        "can't read complex variant with fusion and 2 variants" {
            CivicMultipleEventsReader.read(fusionAnd2Variants) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with only knowledgebase variant" {
            CivicMultipleEventsReader.read(CivicKnowledgebaseVariantReaderTest.variant) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with only cdna variant" {
            CivicMultipleEventsReader.read(CivicCDnaAnnotationReaderTest.variant) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with only fusion pair" {
            CivicMultipleEventsReader.read(CivicFusionReaderTest.fusionVariant) shouldBe emptyList<SomaticEvent>()
        }

        "does not read from input with only promiscuous gene" {
            CivicMultipleEventsReader.read(CivicFusionReaderTest.promiscuousGene) shouldBe emptyList<SomaticEvent>()
        }
    }
}
