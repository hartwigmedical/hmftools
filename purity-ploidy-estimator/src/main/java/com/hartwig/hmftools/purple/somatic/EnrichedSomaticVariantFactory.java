package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.purity.PurityAdjuster.purityAdjustedVAF;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipper;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipperAllPositionsHandler;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EnrichedSomaticVariantFactory
        implements SimpleGenomeZipperAllPositionsHandler<PurpleCopyNumber, SomaticVariant> {

    @NotNull
    public static List<EnrichedSomaticVariant> create(@Nullable FittedPurity purity,
            @NotNull final List<SomaticVariant> variants, @NotNull final List<PurpleCopyNumber> copyNumbers) {
        if (purity == null || copyNumbers.isEmpty()) {
            return create(variants);
        }

        EnrichedSomaticVariantFactory handler = new EnrichedSomaticVariantFactory(purity.purity());
        SimpleGenomeZipper.zip(copyNumbers, variants, handler);
        return handler.reports();
    }

    @NotNull
    private static List<EnrichedSomaticVariant> create(@NotNull final List<SomaticVariant> variants) {
        return variants.stream().map(EnrichedSomaticVariantFactory::create).collect(Collectors.toList());
    }

    private final List<EnrichedSomaticVariant> reports = Lists.newArrayList();
    private final double purity;

    private EnrichedSomaticVariantFactory(final double purity) {
        this.purity = purity;
    }

    private List<EnrichedSomaticVariant> reports() {
        return reports;
    }

    @Override
    public void handle(@Nullable final PurpleCopyNumber copyNumber, @NotNull final SomaticVariant variant) {
        if (copyNumber != null) {
            reports.add(enrich(copyNumber, variant));
        } else {
            reports.add(create(variant));
        }
    }

    private EnrichedSomaticVariant enrich(@NotNull final PurpleCopyNumber copyNumber,
            @NotNull final SomaticVariant variant) {
        double adjustedVAF = purityAdjustedVAF(purity, copyNumber.averageTumorCopyNumber(), variant.alleleFrequency());

        return ImmutableEnrichedSomaticVariant.builder()
                .from(variant)
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .adjustedCopyNumber(copyNumber.averageTumorCopyNumber())
                .adjustedVAF(adjustedVAF)
                .build();
    }

    private static EnrichedSomaticVariant create(@NotNull final SomaticVariant variant) {
        return ImmutableEnrichedSomaticVariant.builder()
                .from(variant)
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .build();
    }
}
