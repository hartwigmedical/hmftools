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
                .hotspotStructuralVariantCount(0)
                .smallVariantAlleleReadCount(smallVariantAlleleReadCount(purpleData))
                .structuralVariantTumorFragmentCount(0)
                .bafCount(0)
                .build();
    }

    private static int smallVariantAlleleReadCount(@NotNull PurpleData purpleData)
    {
        return purpleData.somaticVariants().stream()
                .filter(variant -> variant.type() == VariantType.SNP)
                .mapToInt(variant -> variant.allelicDepth().AlleleReadCount)
                .sum();
    }

    private static int hotspotMutationCount(@NotNull PurpleData purpleData)
    {
        return (int) purpleData.somaticVariants().stream()
                .filter(variant -> variant.tier() == VariantTier.HOTSPOT)
                .count();
    }
}
