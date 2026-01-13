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
                .hotspotMutationCount(1)
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
                .smallVariantAlleleReadCount(1)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData);
        assertEquals(expectedStats, tumorStats);
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
                .somaticVariants(somaticVariants.stream().
                        filter(PurpleVariantContext::reported)
                        .collect(Collectors.toList()))
                .germlineVariants(germlineVariants.stream()
                        .filter(PurpleVariantContext::reported)
                        .collect(Collectors.toList()))
                .build();
    }
}
