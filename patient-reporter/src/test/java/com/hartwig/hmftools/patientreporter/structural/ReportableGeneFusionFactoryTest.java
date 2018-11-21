package com.hartwig.hmftools.patientreporter.structural;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableGeneFusionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canConvertGeneFusions() {
        GeneFusion fusion1 =
                createFusion("TMPRSS2", "ENST00000398585", 4, 5, "PNPLA7", "ENST00000406427", 2, 3, KnownFusionsModel.CIVIC, 0.4);
        GeneFusion fusion2 = createFusion("CLCN6", "ENST00000346436", 1, 2, "BRAF", "ENST00000288602", 8, 9, KnownFusionsModel.ONCOKB, 1D);

        List<ReportableGeneFusion> reportableFusions =
                ReportableGeneFusionFactory.toReportableGeneFusions(Lists.newArrayList(fusion1, fusion2));

        assertEquals(2, reportableFusions.size());

        ReportableGeneFusion reportableFusion1 = reportableFusions.get(0);
        assertEquals("TMPRSS2", reportableFusion1.geneStart());
        assertEquals("PNPLA7", reportableFusion1.geneEnd());
        assertEquals(0.4, reportableFusion1.ploidy(), EPSILON);
    }

    @NotNull
    private static GeneFusion createFusion(@NotNull String startGene, @NotNull String startTranscript, int startExonUpstream,
            int startExonDownstream, @NotNull String endGene, @NotNull String endTranscript, int endExonUpstream, int endExonDownstream,
            @NotNull String source, double ploidy) {
        Transcript upstreamTranscript = createFusionLeg(true, startGene, startTranscript, startExonUpstream, startExonDownstream, ploidy);
        Transcript downstreamTranscript = createFusionLeg(false, endGene, endTranscript, endExonUpstream, endExonDownstream, ploidy);

        return new GeneFusion(upstreamTranscript, downstreamTranscript, source, true, true);
    }

    @NotNull
    private static Transcript createFusionLeg(boolean isUpstream, @NotNull String gene, @NotNull String transcript, int exonUpstream,
            int exonDownstream, double ploidy) {
        EnrichedStructuralVariant variant = createEnrichedStructuralVariantBuilder().type(StructuralVariantType.BND)
                .primaryKey(1)
                .start(createEnrichedStructuralVariantLegBuilder().orientation((byte) 1).chromosome("any").build())
                .end(createEnrichedStructuralVariantLegBuilder().orientation((byte) 1).chromosome("any").build())
                .ploidy(ploidy)
                .build();
        GeneAnnotation upstreamGene =
                new GeneAnnotation(variant, isUpstream, gene, Strings.EMPTY, 1, Lists.newArrayList(), Lists.newArrayList(), Strings.EMPTY);
        return new Transcript(upstreamGene, transcript, exonUpstream, -1, exonDownstream, -1, 0, 0, 10, true, 0, 0, null, null);
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