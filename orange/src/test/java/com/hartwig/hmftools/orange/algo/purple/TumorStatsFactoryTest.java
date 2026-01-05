package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.algo.purple.PurpleTestFactory.createMinimalTestPurpleDataBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.sv.StructuralVariantImpl;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.datamodel.purple.ImmutableTumorStats;
import com.hartwig.hmftools.datamodel.purple.TumorStats;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TumorStatsFactoryTest
{
    @Test
    public void canComputeStructuralVariantCount()
    {
        PurpleData purpleData = createMinimalTestPurpleDataBuilder().addAllPassingSomaticStructuralVariants(
                        withStartFragmentCount(createStructuralVariant(10, 20, StructuralVariantType.INS, false), 3),
                        withStartFragmentCount(createStructuralVariant(50, 60, StructuralVariantType.DEL, false), 5)
                )
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .structuralVariantTumorFragmentCount(8)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeHotspotStructuralVariantCount()
    {
        StructuralVariantImpl hotspotSV = createStructuralVariant(10, 20, StructuralVariantType.INS, true);
        StructuralVariantImpl nonHotspotSV1 = createStructuralVariant(30, 40, StructuralVariantType.INS, false);
        StructuralVariantImpl nonHotspotSV2 = createStructuralVariant(50, 60, StructuralVariantType.DEL, false);

        PurpleData purpleData = createMinimalTestPurpleDataBuilder()
                .addAllPassingSomaticStructuralVariants(hotspotSV, nonHotspotSV1, nonHotspotSV2)
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .hotspotStructuralVariantCount(1)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeHotspotMutationCount()
    {
        List<PurpleVariantContext> somaticVariants = List.of(
                purpleVariantContext(VariantTier.HOTSPOT, true),
                purpleVariantContext(VariantTier.HOTSPOT, false),
                purpleVariantContext(VariantTier.HIGH_CONFIDENCE, true),
                purpleVariantContext(VariantTier.LOW_CONFIDENCE, true)
        );

        List<PurpleVariantContext> germLineVariants = List.of(
                purpleVariantContext(VariantTier.HOTSPOT, true),
                purpleVariantContext(VariantTier.HOTSPOT, false)
        );

        PurpleData purpleData = createTestPurpleData(somaticVariants, germLineVariants);

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .hotspotMutationCount(2)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeSmallVariantAlleleReadCount()
    {
        List<PurpleVariantContext> somaticVariants = List.of(
                purpleVariantContext(VariantType.SNP, true, 1),
                purpleVariantContext(VariantType.SNP, false, 2),
                purpleVariantContext(VariantType.MNP, true, 4),
                purpleVariantContext(VariantType.INDEL, true, 8)
        );

        PurpleData purpleData = createTestPurpleData(somaticVariants, List.of());

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .smallVariantAlleleReadCount(3)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeBafCount()
    {
        List<Segment> segments = List.of(
                ImmutableSegment.builder().germlineStatus(
                                GermlineStatus.DIPLOID)
                        .observedTumorRatio(0.7)
                        .bafCount(1)
                        .build(),
                ImmutableSegment.builder().germlineStatus(
                                GermlineStatus.DIPLOID)
                        .observedTumorRatio(1.0)
                        .bafCount(2)
                        .build(),
                ImmutableSegment.builder().germlineStatus(
                                GermlineStatus.DIPLOID)
                        .observedTumorRatio(1.3)
                        .bafCount(4)
                        .build(),
                ImmutableSegment.builder().germlineStatus(
                                GermlineStatus.NOISE)
                        .observedTumorRatio(0.7)
                        .bafCount(8)
                        .build()
        );

        PurpleData purpleData = createMinimalTestPurpleDataBuilder()
                .segments(segments)
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .bafCount(5)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData);
        assertEquals(expectedStats, tumorStats);
    }

    @NotNull
    private StructuralVariantImpl createStructuralVariant(int start, int end, StructuralVariantType type, boolean hotspot)
    {
        return PurpleTestUtils.createStructuralVariant(
                        "chr1",
                        start,
                        "chr1",
                        end,
                        type,
                        0.0,
                        0.0)
                .hotspot(hotspot)
                .build();
    }

    private static StructuralVariantImpl withStartFragmentCount(StructuralVariantImpl variant, int count)
    {
        return ImmutableStructuralVariantImpl.builder()
                .from(variant)
                .start(ImmutableStructuralVariantLegImpl.builder().from(variant.start())
                        .tumorVariantFragmentCount(count).build())
                .build();
    }

    @NotNull
    private static PurpleVariantContext purpleVariantContext(VariantType variantType, boolean reported, int alleleReadCount)
    {
        return TestPurpleVariantFactory.contextBuilder()
                .from(TestPurpleVariantFactory.variantBuilder()
                        .type(variantType)
                        .allelicDepth(new AllelicDepth(100, alleleReadCount))
                        .reported(reported)
                        .build())
                .build();
    }

    @NotNull
    private static PurpleVariantContext purpleVariantContext(VariantTier variantTier, boolean reported)
    {
        return TestPurpleVariantFactory.contextBuilder()
                .from(TestPurpleVariantFactory.variantBuilder()
                        .tier(variantTier)
                        .reported(reported)
                        .build()
                )
                .build();
    }

    @NotNull
    public static ImmutableTumorStats.Builder createMinimalTumorStatsBuilder()
    {
        return ImmutableTumorStats.builder()
                .maxDiploidProportion(0)
                .hotspotMutationCount(0)
                .hotspotStructuralVariantCount(0)
                .smallVariantAlleleReadCount(0)
                .structuralVariantTumorFragmentCount(0)
                .bafCount(0);
    }

    @NotNull
    private PurpleData createTestPurpleData(List<PurpleVariantContext> somaticVariants, List<PurpleVariantContext> germlineVariants)
    {
        return createMinimalTestPurpleDataBuilder()
                .allSomaticVariants(somaticVariants)
                .driverSomaticVariants(somaticVariants.stream().
                        filter(PurpleVariantContext::reported)
                        .collect(Collectors.toList()))
                .allGermlineVariants(germlineVariants)
                .driverGermlineVariants(germlineVariants.stream()
                        .filter(PurpleVariantContext::reported)
                        .collect(Collectors.toList()))
                .build();
    }
}
