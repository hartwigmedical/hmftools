package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import java.util.List;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.datamodel.purple.ImmutableTumorStats;
import com.hartwig.hmftools.datamodel.purple.TumorStats;

import org.jetbrains.annotations.NotNull;

public class TumorStatsFactory
{

    @NotNull
    public static TumorStats compute(@NotNull PurpleData purpleData)
    {
        return ImmutableTumorStats.builder()
                .hotspotMutationCount(hotspotMutationCount(purpleData))
                .hotspotStructuralVariantCount(hotspotStructuralVariants(purpleData))
                .smallVariantCount(smallVariantCount(purpleData))
                .structuralVariantsCount(structuralVariantCount(purpleData))
                .sumBafCounts(sumBafCounts(purpleData))
                .build();
    }

    private static int structuralVariantCount(@NotNull PurpleData purpleData)
    {
        List<EnrichedStructuralVariant> enrichedVariants =
                new EnrichedStructuralVariantFactory().enrich(purpleData.allSomaticStructuralVariants());

        // also check the need for PASS filtering here

        return enrichedVariants.stream()
                .filter(variant -> variant.filter().equals(PASS) && variant.type() != StructuralVariantType.SGL)
                .mapToInt(variant -> variant.start().tumorVariantFragmentCount())
                .sum();
    }

    private static int smallVariantCount(@NotNull PurpleData purpleData)
    {
        // TODO verify these are already filtered for PASS?

        return purpleData.allSomaticVariants().stream()
                .filter(variant -> variant.type() == VariantType.SNP)
                .mapToInt(variant -> variant.allelicDepth().AlleleReadCount)
                .sum();
    }

    private static int hotspotMutationCount(@NotNull PurpleData purpleData)
    {
        return (int) purpleData.allSomaticVariants().stream()
                .filter(variant -> variant.tier() == VariantTier.HOTSPOT)
                .count();
    }

    private static int sumBafCounts(@NotNull PurpleData purpleData)
    {

        return purpleData.segments().stream()
                .filter(segment -> segment.germlineStatus() == GermlineStatus.DIPLOID)
                .filter(segment -> segment.observedTumorRatio() < 0.8 || segment.observedTumorRatio() > 1.2)
                .mapToInt(Segment::bafCount)
                .sum();
    }

    private static int hotspotStructuralVariants(@NotNull PurpleData purpleData)
    {
        List<EnrichedStructuralVariant> enrichedVariants =
                new EnrichedStructuralVariantFactory().enrich(purpleData.allSomaticStructuralVariants());

        // also check the need for PASS filtering here

        return (int) enrichedVariants.stream()
                .filter(variant -> variant.filter().equals(PASS) && variant.hotspot())
                .count();
    }
}
