package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestCopyNumberBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableGeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.Disruption;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.ImmutableDisruption;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableGeneDisruptionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    public Disruption disruptionTestData() {
        return ImmutableDisruption.builder()
                .reportable(true)
                .svId("755779")
                .chromosome("3")
                .position("125593804")
                .orientation("-1")
                .type("INV")
                .gene("ROPN1B,")
                .transcript("ENST00000514116")
                .strand("1")
                .regionType("Upstream")
                .codingType("5P_UTR")
                .biotype("protein_coding")
                .exon("0")
                .isDisruptive(false)
                .build();
    }

    @Test
    public void canConvertPairedDisruption() {
        List<GeneCopyNumber> copyNumbers =
                Lists.newArrayList(createTestCopyNumberBuilder().gene("ERBB4").minCopyNumber(1).maxCopyNumber(1).build());
        List<Disruption> pairedDisruptions = Lists.newArrayList(disruptionTestData());
//
//        List<ReportableGeneDisruption> reportableDisruptions =
//                ReportableGeneDisruptionFactory.disruptionConvertGeneDisruption(pairedDisruptions, copyNumbers);
//
//        assertEquals(1, reportableDisruptions.size());
//
//        ReportableGeneDisruption disruption = reportableDisruptions.get(0);
//        assertEquals("INV", disruption.type());
//        assertEquals("2q34", disruption.location());
//        assertEquals("ERBB4", disruption.gene());
//        assertEquals("Intron 4 -> Intron 9", disruption.range());
//        assertEquals(Integer.valueOf(1), disruption.geneMinCopies());
//        assertEquals(Integer.valueOf(1), disruption.geneMaxCopies());
//        assertEquals(4, disruption.firstAffectedExon());
//        assertEquals(1D, disruption.ploidy(), EPSILON);
    }

    @Test
    public void doesNotPairDisruptionsOnDifferentGenes() {
        List<GeneCopyNumber> copyNumbers =
                Lists.newArrayList(createTestCopyNumberBuilder().gene("ERBB4").minCopyNumber(1).maxCopyNumber(1).build(),
                        createTestCopyNumberBuilder().gene("ERBB2").minCopyNumber(1).maxCopyNumber(1).build());
        List<Disruption> pairedDisruptions = Lists.newArrayList(disruptionTestData());
//
//        List<ReportableGeneDisruption> reportableDisruptions =
//                ReportableGeneDisruptionFactory.disruptionConvertGeneDisruption(pairedDisruptions, copyNumbers);
//
//        assertEquals(2, reportableDisruptions.size());
    }

    @Test
    public void canConvertNormalDisruptionsWithoutCopyNumbers() {
        List<Disruption> disruption1 = Lists.newArrayList(disruptionTestData());
        List<Disruption> disruption2 = Lists.newArrayList(disruptionTestData());
        List<Disruption> disruption3 = Lists.newArrayList(disruptionTestData());
//
//        List<GeneDisruption> disruptions = Lists.newArrayList(disruption1, disruption2, disruption3);
//
//        List<ReportableGeneDisruption> reportableDisruptions =
//                ReportableGeneDisruptionFactory.disruptionConvertGeneDisruption(disruptions, Lists.newArrayList());
//        assertEquals(3, reportableDisruptions.size());
    }

    @NotNull
    private static List<GeneDisruption> createTestDisruptionPair(@NotNull StructuralVariantType type, @NotNull String chromosome,
            @NotNull String chromosomeBand, @NotNull String upstreamGene, @NotNull String downstreamGene, int startExon, int endExon,
            double ploidy) {
        EnrichedStructuralVariantLeg start =
                createEnrichedStructuralVariantLegBuilder().orientation((byte) 1).chromosome(chromosome).build();
        EnrichedStructuralVariantLeg end = createEnrichedStructuralVariantLegBuilder().orientation((byte) 0).chromosome(chromosome).build();
        EnrichedStructuralVariant variant =
                createEnrichedStructuralVariantBuilder().type(type).start(start).end(end).ploidy(ploidy).build();

        GeneAnnotation upstreamGeneAnnotation =
                new GeneAnnotation(variant, true, upstreamGene, "id", 1, Lists.newArrayList(), Lists.newArrayList(), chromosomeBand);

        Transcript upstreamTranscript =
                new Transcript(upstreamGeneAnnotation, 0,"trans", startExon, -1, startExon + 1, -1, 0, 0, 15, true, 0, 0, null, null);

        GeneDisruption upstreamDisruption = ImmutableGeneDisruption.builder().reportable(true).linkedAnnotation(upstreamTranscript).build();

        GeneAnnotation downstreamGeneAnnotation =
                new GeneAnnotation(variant, false, downstreamGene, "id", 1, Lists.newArrayList(), Lists.newArrayList(), chromosomeBand);

        Transcript downstreamTranscript =
                new Transcript(downstreamGeneAnnotation, 0,"trans", endExon, -1, endExon + 1, -1, 0, 0, 15, true, 0, 0, null, null);

        GeneDisruption downstreamDisruption =
                ImmutableGeneDisruption.builder().reportable(true).linkedAnnotation(downstreamTranscript).build();

        return Lists.newArrayList(upstreamDisruption, downstreamDisruption);
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

        Transcript transcript =
                new Transcript(geneAnnotation, 0,"trans", exonUpstream, -1, exonUpstream + 1, -1, 0, 0, 5, true, 0, 0, null, null);

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
                .insertSequenceAlignments(Strings.EMPTY)
                .qualityScore(0)
                .recovered(false);
    }
}