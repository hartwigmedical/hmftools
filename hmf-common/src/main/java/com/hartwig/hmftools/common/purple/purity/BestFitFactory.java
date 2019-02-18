package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class BestFitFactory {

    private static final double PERCENT_RANGE = 0.1;
    private static final double ABS_RANGE = 0.0005;

    private final double minSomaticUnadjustedVaf;
    private final double highlyDiploidPercentage;
    private final int minVariants;
    private double minSomaticPurity;

    @NotNull
    private final BestFit bestFit;

    public BestFitFactory(double minSomaticUnadjustedVaf, int minVariants, int minPeak, double highlyDiploidPercentage,
            double minSomaticPurity, double minSomaticPuritySpread, @NotNull final List<FittedPurity> fittedPurities,
            @NotNull final List<SomaticVariant> somatics) {
        assert (!fittedPurities.isEmpty());
        this.minSomaticUnadjustedVaf = minSomaticUnadjustedVaf;
        this.minVariants = minVariants;
        this.minSomaticPurity = minSomaticPurity;
        this.highlyDiploidPercentage = highlyDiploidPercentage;

        long somaticsWithSufficientVafCount = somaticsWithSufficientVaf(somatics);

        Collections.sort(fittedPurities);
        FittedPurity lowestScore = fittedPurities.get(0);

        final FittedPurityScore score = FittedPurityScoreFactory.score(inRangeOfLowest(lowestScore.score(), fittedPurities));

        final FittedPurity fit;
        final FittedPurityStatus status;
        if (Doubles.greaterOrEqual(score.puritySpread(), minSomaticPuritySpread) && isHighlyDiploid(score)) {
            final Optional<FittedPurity> somaticFit = new SomaticFitFactory(minPeak).fromSomatics(fittedPurities, somatics);

            if (noDetectableTumor(somaticsWithSufficientVafCount)) {
                status = FittedPurityStatus.NO_TUMOR;
                fit = somaticFit.orElse(lowestScore);
            } else if (somaticsWontHelp(somatics.size(), lowestScore.purity(), somaticFit)) {
                status = FittedPurityStatus.HIGHLY_DIPLOID;
                fit = lowestScore;
            } else {
                status = FittedPurityStatus.SOMATIC;
                assert somaticFit.isPresent();
                fit = somaticFit.get();
            }

        } else {
            status = FittedPurityStatus.NORMAL;
            fit = lowestScore;
        }

        bestFit = ImmutableBestFit.builder().fit(fit).status(status).score(score).bestFitPerPurity(fittedPurities).build();
    }

    private long somaticsWithSufficientVaf(@NotNull Collection<SomaticVariant> variants) {
        return variants.stream().filter(x -> Doubles.greaterOrEqual(x.alleleFrequency(), minSomaticUnadjustedVaf)).count();
    }

    private boolean noDetectableTumor(long somaticCount) {
        return somaticCount > 0 && somaticCount < minVariants;
    }

    private boolean somaticsWontHelp(int somaticCount, double lowestScoringPurity, @NotNull final Optional<FittedPurity> somaticFit) {
        return !somaticFit.isPresent() || somaticCount == 0 || (Doubles.lessThan(lowestScoringPurity, minSomaticPurity)
                && Doubles.lessThan(somaticFit.get().purity(), minSomaticPurity));
    }

    private boolean isHighlyDiploid(@NotNull final FittedPurityScore score) {
        return Doubles.greaterOrEqual(score.maxDiploidProportion(), highlyDiploidPercentage);
    }

    @NotNull
    public BestFit bestFit() {
        return bestFit;
    }

    @NotNull
    private static List<FittedPurity> inRangeOfLowest(double lowestScore, @NotNull final List<FittedPurity> purities) {
        return purities.stream().filter(inRangeOfLowest(lowestScore)).collect(Collectors.toList());
    }

    @NotNull
    private static Predicate<FittedPurity> inRangeOfLowest(final double score) {
        return fittedPurity -> {
            double absDifference = Math.abs(fittedPurity.score() - score);
            double relDifference = Math.abs(absDifference / score);
            return lessOrEqual(absDifference, ABS_RANGE) || lessOrEqual(relDifference, PERCENT_RANGE);
        };
    }
}
