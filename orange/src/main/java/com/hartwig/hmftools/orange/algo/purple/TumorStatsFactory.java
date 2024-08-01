package com.hartwig.hmftools.orange.algo.purple;

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
    public static TumorStats compute(@NotNull PurpleData purpleData, boolean convertGermlineToSomatic)
    {
        return ImmutableTumorStats.builder()
                .maxDiploidProportion(purpleData.purityContext().score().maxDiploidProportion())
                .hotspotMutationCount(hotspotMutationCount(purpleData, convertGermlineToSomatic))
                .hotspotStructuralVariantCount(hotspotStructuralVariants(purpleData))
                .smallVariantAlleleReadCount(smallVariantAlleleReadCount(purpleData))
                .structuralVariantTumorFragmentCount(structuralVariantTumorFragmentCount(purpleData))
                .bafCount(bafCount(purpleData))
                .build();
    }

    private static int structuralVariantTumorFragmentCount(@NotNull PurpleData purpleData)
    {
        List<EnrichedStructuralVariant> enrichedVariants =
                new EnrichedStructuralVariantFactory().enrich(purpleData.allSomaticStructuralVariants());

        return enrichedVariants.stream()
                .filter(variant -> !variant.isFiltered() && variant.type() != StructuralVariantType.SGL)
                .mapToInt(variant ->
                {
                    Integer count = variant.start().tumorVariantFragmentCount();
                    return count != null ? count : 0;
                })
                .sum();
    }

    private static int smallVariantAlleleReadCount(@NotNull PurpleData purpleData)
    {
        return purpleData.reportableSomaticVariants().stream()
                .filter(variant -> variant.type() == VariantType.SNP)
                .mapToInt(variant -> variant.allelicDepth().AlleleReadCount)
                .sum();
    }

    private static int hotspotMutationCount(@NotNull PurpleData purpleData, boolean convertGermlineToSomatic)
    {
        int count = (int) purpleData.allSomaticVariants().stream()
                .filter(variant -> variant.tier() == VariantTier.HOTSPOT)
                .count();

        if(convertGermlineToSomatic && purpleData.reportableGermlineVariants() != null)
        {
            count += (int) purpleData.reportableGermlineVariants().stream()
                    .filter(variant -> variant.tier() == VariantTier.HOTSPOT)
                    .count();
        }

        return count;
    }

    private static int bafCount(@NotNull PurpleData purpleData)
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

        return (int) enrichedVariants.stream()
                .filter(variant -> !variant.isFiltered() && variant.hotspot())
                .count();
    }
}
