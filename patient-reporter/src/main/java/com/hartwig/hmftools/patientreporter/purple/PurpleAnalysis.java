package com.hartwig.hmftools.patientreporter.purple;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.patientreporter.variants.ImmutableVariantReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleAnalysis {

    @NotNull
    public abstract FittedPurity fittedPurity();

    @NotNull
    public abstract FittedPurityScore fittedScorePurity();

    @NotNull
    public abstract List<PurpleCopyNumber> copyNumbers();

    @NotNull
    public List<CopyNumber> ploidyAdjustedCopyNumbers() {
        return copyNumbers().stream().map(x -> ploidyAdjusted(fittedPurity().ploidy(), x)).collect(Collectors.toList());
    }

    @NotNull
    public List<VariantReport> enrich(@NotNull List<VariantReport> variants) {

        final List<VariantReport> result = Lists.newArrayList();
        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, fittedPurity());
        final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers());

        for (VariantReport variant : variants) {
            final Optional<PurpleCopyNumber> optionalCopyNumber = copyNumberSelector.select(variant);
            if (optionalCopyNumber.isPresent()) {
                PurpleCopyNumber copyNumber = optionalCopyNumber.get();
                double adjustedVAF =
                        Math.min(1, purityAdjuster.purityAdjustedVAF(copyNumber.averageTumorCopyNumber(), variant.alleleFrequency()));
                result.add(ImmutableVariantReport.builder().from(variant).baf(copyNumber.descriptiveBAF()).impliedVAF(adjustedVAF).build());
            } else {
                result.add(variant);
            }
        }

        return result;
    }

    public double purityUncertainty() {
        return fittedScorePurity().maxPurity() - fittedScorePurity().minPurity();
    }

    @NotNull
    @VisibleForTesting
    static CopyNumber ploidyAdjusted(final double ploidy, final @NotNull PurpleCopyNumber copyNumber) {
        double adjustedCopyNumber =
                copyNumber.value() <= 1 ? copyNumber.averageTumorCopyNumber() : copyNumber.averageTumorCopyNumber() / ploidy * 2;

        return ImmutablePurpleCopyNumber.builder().from(copyNumber).averageTumorCopyNumber(adjustedCopyNumber).build();
    }
}
