package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.algo.purple.PurpleTestFactory.createMinimalTestPurpleData;
import static com.hartwig.hmftools.orange.algo.purple.PurpleTestFactory.createMinimalTestPurpleDataBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

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
    public void canComputeZeroStatsFromMinimalPurpleData()
    {
        PurpleData purpleData = createMinimalTestPurpleData();
        TumorStats expectedStats = createMinimalTumorStatsBuilder().build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData, false);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeStructuralVariantCount()
    {
        StructuralVariantImpl sv1 =
                withFragmentCount(PurpleTestUtils.createStructuralVariant("chr1", 30, "chr1", 40, StructuralVariantType.INS, 0.0, 0.0)
                        .hotspot(false)
                        .build(), 3);

        StructuralVariantImpl sv2 =
                withFragmentCount(PurpleTestUtils.createStructuralVariant("chr1", 50, "chr1", 60, StructuralVariantType.DEL, 0.0, 0.0)
                        .hotspot(false)
                        .build(), 5);

        PurpleData purpleData = createMinimalTestPurpleDataBuilder()
                .addAllSomaticStructuralVariants(sv1, sv2)
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .structuralVariantTumorFragmentCount(8)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData, false);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeHotspotStructuralVariantCount()
    {
        StructuralVariantImpl hotspotSV =
                PurpleTestUtils.createStructuralVariant("chr1", 10, "chr1", 20, StructuralVariantType.INS, 0.0, 0.0)
                        .hotspot(true)
                        .build();

        StructuralVariantImpl nonHotspotSV1 =
                PurpleTestUtils.createStructuralVariant("chr1", 30, "chr1", 40, StructuralVariantType.INS, 0.0, 0.0)
                        .hotspot(false)
                        .build();

        StructuralVariantImpl nonHotspotSV2 =
                PurpleTestUtils.createStructuralVariant("chr1", 50, "chr1", 60, StructuralVariantType.DEL, 0.0, 0.0)
                        .hotspot(false)
                        .build();

        PurpleData purpleData = createMinimalTestPurpleDataBuilder()
                .addAllSomaticStructuralVariants(hotspotSV, nonHotspotSV1, nonHotspotSV2)
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .hotspotStructuralVariantCount(1)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData, false);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeHotspotMutationCount()
    {
        PurpleVariantContext hostpotVariant = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.HOTSPOT)
                .build();

        PurpleVariantContext germlineHotspotVariant = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.HOTSPOT)
                .build();

        PurpleVariantContext otherVariant1 = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.HIGH_CONFIDENCE)
                .build();

        PurpleVariantContext otherVariant2 = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.LOW_CONFIDENCE)
                .build();

        PurpleData purpleData = createMinimalTestPurpleDataBuilder()
                .addAllSomaticVariants(hostpotVariant, otherVariant1, otherVariant2)
                .addAllGermlineVariants(germlineHotspotVariant)
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .hotspotMutationCount(1)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData, false);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canConvertGermlineToSomaticForHotspotMutationCount()
    {
        PurpleVariantContext hostpotVariant = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.HOTSPOT)
                .build();

        PurpleVariantContext germlineHotspotVariant = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.HOTSPOT)
                .build();

        PurpleVariantContext otherVariant1 = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.HIGH_CONFIDENCE)
                .build();

        PurpleVariantContext otherVariant2 = TestPurpleVariantFactory.contextBuilder()
                .tier(VariantTier.LOW_CONFIDENCE)
                .build();

        PurpleData purpleData = createMinimalTestPurpleDataBuilder()
                .addAllSomaticVariants(hostpotVariant, otherVariant1, otherVariant2)
                .addAllGermlineVariants(germlineHotspotVariant)
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .hotspotMutationCount(2)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData, true);
        assertEquals(expectedStats, tumorStats);
    }

    @Test
    public void canComputeSmallVariantAlleleReadCount()
    {
        PurpleVariantContext snpVariant = TestPurpleVariantFactory.contextBuilder()
                .type(VariantType.SNP)
                .allelicDepth(new AllelicDepth(10, 2))
                .build();

        PurpleVariantContext otherVariant1 = TestPurpleVariantFactory.contextBuilder()
                .type(VariantType.MNP)
                .allelicDepth(new AllelicDepth(10, 4))
                .build();

        PurpleVariantContext otherVariant2 = TestPurpleVariantFactory.contextBuilder()
                .type(VariantType.INDEL)
                .allelicDepth(new AllelicDepth(10, 8))
                .build();

        PurpleData purpleData = createMinimalTestPurpleDataBuilder()
                .addAllSomaticVariants(snpVariant, otherVariant1, otherVariant2)
                .build();

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .smallVariantAlleleReadCount(2)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData, false);
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

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData, false);
        assertEquals(expectedStats, tumorStats);
    }

    public static StructuralVariantImpl withFragmentCount(StructuralVariantImpl variant, int count)
    {
        return ImmutableStructuralVariantImpl.builder()
                .from(variant)
                .start(ImmutableStructuralVariantLegImpl.builder().from(variant.start())
                        .tumorVariantFragmentCount(count).build())
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
}
