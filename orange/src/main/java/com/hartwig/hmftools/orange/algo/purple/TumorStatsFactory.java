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
import com.hartwig.hmftools.orange.algo.gripss.GripssData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorStatsFactory
{

    @NotNull
    public static TumorStats compute(@NotNull PurpleData purpleData, @Nullable GripssData gripsData)
    {
        return ImmutableTumorStats.builder()
                .numberHotspotMutations(hotspotMutations(purpleData))
                .numberHotspotSVs(gripss(gripsData))
                .sumSNPAlleleReadCounts(smallvariants(purpleData))
                .sumTumorVariantFragmentCountsExclSGL(structuralvariants(purpleData))
                .sumBafCount(bafCount(purpleData))
                .build();
    }

    public static int structuralvariants(@NotNull PurpleData purpleData)
    {
        List<EnrichedStructuralVariant> enrichedVariants =
                new EnrichedStructuralVariantFactory().enrich(purpleData.allSomaticStructuralVariants());

        // also check the need for PASS filtering here
        int svCount = enrichedVariants.stream()
                .filter(variant -> variant.filter().equals(PASS) && variant.type() != StructuralVariantType.SGL)
                .mapToInt(variant -> variant.start().tumorVariantFragmentCount())
                .sum();

        return svCount;
    }

    public static int smallvariants(@NotNull PurpleData purpleData)
    {
        // TODO verify these are already filtered for PASS?

        int snvCount = purpleData.allSomaticVariants().stream()
                .filter(variant -> variant.type() == VariantType.SNP)
                .mapToInt(variant -> variant.allelicDepth().AlleleReadCount)
                .sum();
        return snvCount;
    }

    public static int hotspotMutations(@NotNull PurpleData purpleData)
    {
        int hotspotCount = (int) purpleData.allSomaticVariants().stream()
                .filter(variant -> variant.tier() == VariantTier.HOTSPOT)
                .count();
        return hotspotCount;
    }

    public static int bafCount(@NotNull PurpleData purpleData)
    {

        int sumBafCount = purpleData.segments().stream()
                .filter(segment -> segment.germlineStatus() == GermlineStatus.DIPLOID)
                .filter(segment -> segment.observedTumorRatio() < 0.8 || segment.observedTumorRatio() > 1.2)
                .mapToInt(Segment::bafCount)
                .sum();

        return sumBafCount;
    }

    public static int gripss(@NotNull GripssData gripsData)
    {
        // TODO either hotspot filtering here or on load ... direct logic is;
        //        CompoundFilter filter = new CompoundFilter(true);
        //        filter.add(new PassingVariantFilter());
        //        VariantContextFilter hotspotFilter = record -> record.hasAttribute("HOTSPOT");
        //        filter.add(hotspotFilter);
        return gripsData.allSomaticStructuralVariants().size();
    }
}
