package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.algo.purple.PurpleTestFactory.createMinimalTestPurpleDataBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.SmallVariant;
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
        List<SmallVariant> somaticVariants = List.of(
                buildPurpleVariant(VariantTier.HOTSPOT, true),
                buildPurpleVariant(VariantTier.HOTSPOT, false),
                buildPurpleVariant(VariantTier.HIGH_CONFIDENCE, true),
                buildPurpleVariant(VariantTier.LOW_CONFIDENCE, true)
        );

        List<SmallVariant> germLineVariants = List.of(
                buildPurpleVariant(VariantTier.HOTSPOT, true),
                buildPurpleVariant(VariantTier.HOTSPOT, false)
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
        List<SmallVariant> somaticVariants = List.of(
                buildPurpleVariant(VariantType.SNP, true, 1),
                buildPurpleVariant(VariantType.SNP, false, 2),
                buildPurpleVariant(VariantType.MNP, true, 4),
                buildPurpleVariant(VariantType.INDEL, true, 8)
        );

        PurpleData purpleData = createTestPurpleData(somaticVariants, List.of());

        TumorStats expectedStats = createMinimalTumorStatsBuilder()
                .smallVariantAlleleReadCount(1)
                .build();

        TumorStats tumorStats = TumorStatsFactory.compute(purpleData);
        assertEquals(expectedStats, tumorStats);
    }

    @NotNull
    private static SmallVariant buildPurpleVariant(VariantType variantType, boolean reported, int alleleReadCount)
    {
        return TestPurpleVariantFactory.variantBuilder()
                        .type(variantType)
                        .allelicDepth(new AllelicDepth(100, alleleReadCount))
                        .reported(reported)
                        .build();
    }

    @NotNull
    private static SmallVariant buildPurpleVariant(VariantTier variantTier, boolean reported)
    {
        return TestPurpleVariantFactory.variantBuilder()
                .tier(variantTier)
                .reported(reported)
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
    private PurpleData createTestPurpleData(List<SmallVariant> somaticVariants, List<SmallVariant> germlineVariants)
    {
        return createMinimalTestPurpleDataBuilder()
                .somaticVariants(somaticVariants.stream().
                        filter(SmallVariant::reported)
                        .collect(Collectors.toList()))
                .germlineVariants(germlineVariants.stream()
                        .filter(SmallVariant::reported)
                        .collect(Collectors.toList()))
                .build();
    }
}
