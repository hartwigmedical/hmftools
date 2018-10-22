package com.hartwig.hmftools.patientreporter.disruption;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableGeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableGeneDisruptionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canConvertPairedDisruption() {
        Pair<GeneDisruption, GeneDisruption> pairedDisruption =
                createTestDisruptionPair(StructuralVariantType.INV, "2", "q34", "ERBB4", 4, 9, 1D);

        List<ReportableGeneDisruption> reportableDisruptions =
                ReportableGeneDisruptionFactory.toReportableGeneDisruptions(Lists.newArrayList(pairedDisruption.getLeft(),
                        pairedDisruption.getRight()));

        assertEquals(1, reportableDisruptions.size());

        ReportableGeneDisruption disruption = reportableDisruptions.get(0);
        assertEquals(StructuralVariantType.INV, disruption.type());
        assertEquals("2q34", disruption.location());
        assertEquals("ERBB4", disruption.gene());
        assertEquals("Intron 4 -> Intron 9", disruption.range());
        assertEquals(4, disruption.firstAffectedExon());
        assertEquals(1D, disruption.ploidy(), EPSILON);
    }

    @Test
    private void canConvertNormalDisruptions() {
        GeneDisruption disruption1 = createDisruption(StructuralVariantType.BND, "17", "q12", "CDK12", 12, 2.3, false);
        GeneDisruption disruption2 = createDisruption(StructuralVariantType.INS, "21", "q22.12", "RUNX1", 0, 0.8, true);
        GeneDisruption disruption3 = createDisruption(StructuralVariantType.DUP, "1", "p13.1", "CD58", 2, 0.2, true);

        List<GeneDisruption> disruptions = Lists.newArrayList(disruption1, disruption2, disruption3);

        List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.toReportableGeneDisruptions(disruptions);
        assertEquals(3, reportableDisruptions.size());
    }

    @NotNull
    private static Pair<GeneDisruption, GeneDisruption> createTestDisruptionPair(@NotNull StructuralVariantType type,
            @NotNull String chromosome, @NotNull String chromosomeBand, @NotNull String gene, int startExon, int endExon, double ploidy) {
        EnrichedStructuralVariantLeg start =
                createEnrichedStructuralVariantLegBuilder().orientation((byte) 1).chromosome(chromosome).build();
        EnrichedStructuralVariantLeg end = createEnrichedStructuralVariantLegBuilder().orientation((byte) 0).chromosome(chromosome).build();
        EnrichedStructuralVariant variant =
                createEnrichedStructuralVariantBuilder().type(type).start(start).end(end).ploidy(ploidy).build();

        GeneAnnotation upstreamGeneAnnotation =
                new GeneAnnotation(variant, true, gene, "id", 1, Lists.newArrayList(), Lists.newArrayList(), chromosomeBand);

        Transcript upstreamTranscript =
                new Transcript(upstreamGeneAnnotation, "trans", startExon, -1, startExon + 1, -1, 15, true, null, null);

        GeneDisruption upstreamDisruption = ImmutableGeneDisruption.builder().reportable(true).linkedAnnotation(upstreamTranscript).build();

        GeneAnnotation downstreamGeneAnnotation =
                new GeneAnnotation(variant, false, gene, "id", 1, Lists.newArrayList(), Lists.newArrayList(), chromosomeBand);

        Transcript downstreamTranscript =
                new Transcript(downstreamGeneAnnotation, "trans", endExon, -1, endExon + 1, -1, 15, true, null, null);

        GeneDisruption downstreamDisruption =
                ImmutableGeneDisruption.builder().reportable(true).linkedAnnotation(downstreamTranscript).build();

        return Pair.of(upstreamDisruption, downstreamDisruption);
    }

    @NotNull
    private static GeneDisruption createDisruption(@NotNull StructuralVariantType type, @NotNull String chromosome,
            @NotNull String chromosomeBand, @NotNull String gene, int exonUpstream, double ploidy, boolean isUpstream) {
        byte orientation = isUpstream ? (byte) 1 : (byte) 0;
        EnrichedStructuralVariantLeg leg =
                createEnrichedStructuralVariantLegBuilder().orientation(orientation).chromosome(chromosome).build();
        ImmutableEnrichedStructuralVariant.Builder builder = createEnrichedStructuralVariantBuilder().type(type).ploidy(ploidy);
        EnrichedStructuralVariant variant = builder.start(leg).build();

        GeneAnnotation geneAnnotation =
                new GeneAnnotation(variant, true, gene, "id", 1, Lists.newArrayList(), Lists.newArrayList(), chromosomeBand);

        Transcript transcript = new Transcript(geneAnnotation, "trans", exonUpstream, -1, exonUpstream + 1, -1, 5, true, null, null);

        return ImmutableGeneDisruption.builder().reportable(true).linkedAnnotation(transcript).build();
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