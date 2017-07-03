package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.Builder;
import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.builder;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class EnrichedSomaticVariantFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final GenomeRegionSelector<GenomeRegion> highConfidenceSelector;
    @NotNull
    private final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector;

    public EnrichedSomaticVariantFactory(double purity, double normFactor,
            @NotNull  final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull  final Multimap<String, PurpleCopyNumber> copyNumbers) {
        purityAdjuster = new PurityAdjuster(Gender.MALE, purity, normFactor);
        highConfidenceSelector = new GenomeRegionSelector<>(highConfidenceRegions);
        copyNumberSelector = new GenomeRegionSelector<>(copyNumbers);
    }

    public List<EnrichedSomaticVariant> enrich(final List<SomaticVariant> variants) {
        return variants.stream().map(this::enrich).collect(Collectors.toList());
    }

    private EnrichedSomaticVariant enrich(@NotNull final SomaticVariant variant) {
        final Builder builder = createBuilder(variant);

        highConfidenceSelector.select(variant).ifPresent(x -> inHighConfidenceRegion(builder));
        copyNumberSelector.select(variant).ifPresent(x -> addCopyNumber(builder, x, variant.alleleFrequency()));

        return builder.build();
    }

    private Builder createBuilder(@NotNull final SomaticVariant variant) {
        return builder().from(variant)
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .highConfidenceRegion(false)
                .adjustedCopyNumber(0)
                .adjustedVAF(0);
    }

    private Builder addCopyNumber(@NotNull final Builder builder, @NotNull final PurpleCopyNumber copyNumber, double alleleFrequency) {
        double adjustedVAF = purityAdjuster.purityAdjustedVAF(copyNumber.averageTumorCopyNumber(), alleleFrequency);
        return builder.adjustedCopyNumber(copyNumber.averageTumorCopyNumber()).adjustedVAF(adjustedVAF);
    }

    private Builder inHighConfidenceRegion(@NotNull final Builder builder) {
        return builder.highConfidenceRegion(true);
    }
}
