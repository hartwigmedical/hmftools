package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.sv.StructuralVariant;
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
                .maxDiploidProportion(purpleData.purityContext().score().maxDiploidProportion())
                .hotspotMutationCount(hotspotMutationCount(purpleData))
                .hotspotStructuralVariantCount(hotspotStructuralVariants(purpleData))
                .smallVariantAlleleReadCount(smallVariantAlleleReadCount(purpleData))
                .structuralVariantTumorFragmentCount(structuralVariantTumorFragmentCount(purpleData))
                .bafCount(bafCount(purpleData))
                .build();
    }

    private static int structuralVariantTumorFragmentCount(@NotNull PurpleData purpleData)
    {
        int svFragmentReadCount = 0;
        for(StructuralVariant variant : purpleData.allPassingSomaticStructuralVariants())
        {
            if(variant.isFiltered() || variant.type() == StructuralVariantType.SGL)
            {
                continue;
            }

            Integer startTumorVariantFragmentCount = variant.start().tumorVariantFragmentCount();
            if(variant.end() != null && startTumorVariantFragmentCount != null)
            {
                svFragmentReadCount += startTumorVariantFragmentCount;
            }
        }
        return svFragmentReadCount;
    }

    private static int smallVariantAlleleReadCount(@NotNull PurpleData purpleData)
    {
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
        return (int) purpleData.allPassingSomaticStructuralVariants().stream()
                .filter(variant -> !variant.isFiltered() && variant.hotspot())
                .count();
    }
}
