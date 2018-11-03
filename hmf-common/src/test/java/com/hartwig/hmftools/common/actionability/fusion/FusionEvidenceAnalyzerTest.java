package com.hartwig.hmftools.common.actionability.fusion;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeReading;
import com.hartwig.hmftools.common.actionability.cancertype.ImmutableCancerTypeReading;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FusionEvidenceAnalyzerTest {

    @Test
    public void actionabilityWorksFusions() {
        ActionableFusion fusion = ImmutableActionableFusion.builder()
                .fiveGene("BCR")
                .threeGene("ABL1")
                .cancerType("Acute lymphoblastic leukemia")
                .drug("Dasatinib")
                .drugsType("BCR-ABL inhibitor")
                .level("A")
                .reference("BCR__ABL1")
                .response("Responsive")
                .source("CGI")
                .build();

        ActionablePromiscuousFive five = ImmutableActionablePromiscuousFive.builder()
                .cancerType("Acute lymphoblastic leukemia")
                .drug("FORT-1")
                .drugsType("Trial")
                .gene("FGFR1")
                .level("B")
                .reference("EXT8753 (NL62741.031.17)")
                .response("Responsive")
                .source("iclusion")
                .build();

        ActionablePromiscuousThree three = ImmutableActionablePromiscuousThree.builder()
                .cancerType("Acute lymphoblastic leukemia")
                .drug("AZD5363")
                .drugsType("AKT inhibitor")
                .gene("AKT3__.")
                .level("D")
                .reference("AKT3__.")
                .response("Responsive")
                .source("cgi")
                .build();

        FusionEvidenceAnalyzer analyzer =
                new FusionEvidenceAnalyzer((Lists.newArrayList(fusion)), Lists.newArrayList(five), Lists.newArrayList(three));

        CancerTypeReading reading = ImmutableCancerTypeReading.builder().doidSet("4159").cancerType("Skin").build();
        CancerTypeAnalyzer cancerType = new CancerTypeAnalyzer(Lists.newArrayList(reading));

        GeneFusion fusion1 =
                createFusion("TMPRSS2", "ENST00000398585", 4, 5, "PNPLA7", "ENST00000406427", 2, 3, KnownFusionsModel.CIVIC, 0.4);

        assertTrue(analyzer.evidenceForFusion(fusion1, "4159", cancerType).isEmpty());
    }

    @NotNull
    private static GeneFusion createFusion(@NotNull String startGene, @NotNull String startTranscript, int startExonUpstream,
            int startExonDownstream, @NotNull String endGene, @NotNull String endTranscript, int endExonUpstream, int endExonDownstream,
            @NotNull String source, double ploidy) {
        Transcript upstreamTranscript = createFusionLeg(true, startGene, startTranscript, startExonUpstream, startExonDownstream, ploidy);
        Transcript downstreamTranscript = createFusionLeg(false, endGene, endTranscript, endExonUpstream, endExonDownstream, ploidy);

        return ImmutableGeneFusion.builder()
                .reportable(true)
                .primarySource(source)
                .upstreamLinkedAnnotation(upstreamTranscript)
                .downstreamLinkedAnnotation(downstreamTranscript)
                .build();
    }

    @NotNull
    private static Transcript createFusionLeg(boolean isUpstream, @NotNull String gene, @NotNull String transcript, int exonUpstream,
            int exonDownstream, double ploidy) {
        EnrichedStructuralVariant variant = createEnrichedStructuralVariantBuilder().type(StructuralVariantType.BND)
                .start(createEnrichedStructuralVariantLegBuilder().orientation((byte) 1).chromosome("any").build())
                .ploidy(ploidy)
                .build();
        GeneAnnotation upstreamGene =
                new GeneAnnotation(variant, isUpstream, gene, Strings.EMPTY, 1, Lists.newArrayList(), Lists.newArrayList(), Strings.EMPTY);
        return new Transcript(upstreamGene, transcript, exonUpstream, -1, exonDownstream, -1, 10, true, null, null);
    }

    @NotNull
    private static ImmutableEnrichedStructuralVariantLeg.Builder createEnrichedStructuralVariantLegBuilder() {
        return ImmutableEnrichedStructuralVariantLeg.builder().homology(Strings.EMPTY).position(1);
    }

    @NotNull
    private static ImmutableEnrichedStructuralVariant.Builder createEnrichedStructuralVariantBuilder() {
        return ImmutableEnrichedStructuralVariant.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .qualityScore(0)
                .recovered(false);
    }
}