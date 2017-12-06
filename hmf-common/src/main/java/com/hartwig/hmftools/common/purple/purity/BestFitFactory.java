package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.variant.PurpleSomaticVariant;

import org.jetbrains.annotations.NotNull;

public class BestFitFactory {

    private static final double PERCENT_RANGE = 0.1;
    private static final double ABS_RANGE = 0.0005;
    private static final double HIGHLY_DIPLOID_PERCENTAGE = 0.98;
    private static final double LOWEST_SCORE_MIN_PURITY = 0.15;
    private static final double MIN_PURITY_SPREAD = 0.15;

    private final FittedPurity bestFit;
    private final FittedPurityScore score;
    private final FittedPurityStatus status;
    private final int minVariants;

    public BestFitFactory(int minVariants, int minPeak, final List<FittedPurity> fittedPurities, List<PurpleSomaticVariant> somatics) {
        assert (!fittedPurities.isEmpty());
        this.minVariants = minVariants;

        Collections.sort(fittedPurities);
        FittedPurity lowestScore = fittedPurities.get(0);

        final List<FittedPurity> candidates = candidates(lowestScore.score(), fittedPurities);
        score = FittedPurityScoreFactory.score(candidates);

        if (Doubles.greaterOrEqual(score.puritySpread(), MIN_PURITY_SPREAD) && isHighlyDiploid(score)) {
            if (noDetectableTumor(somatics.size())) {
                status = FittedPurityStatus.NO_TUMOR;
                bestFit = lowestScore;
            } else {
                final Optional<FittedPurity> somaticFit = new SomaticFitFactory(minPeak).fromSomatics(candidates, somatics);
                if (somaticsWontHelp(somatics.size(), lowestScore.purity(), somaticFit)) {
                    status = FittedPurityStatus.HIGHLY_DIPLOID;
                    bestFit = lowestScore;
                } else {
                    status = FittedPurityStatus.SOMATIC;
                    assert somaticFit.isPresent();
                    bestFit = somaticFit.get();
                }
            }

        } else {
            status = FittedPurityStatus.NORMAL;
            bestFit = lowestScore;
        }
    }

    private boolean noDetectableTumor(int somaticCount) {
        return somaticCount > 0 && somaticCount < minVariants;
    }

    private boolean somaticsWontHelp(int somaticCount, double lowestScoringPurity, @NotNull final Optional<FittedPurity> somaticFit) {
        return !somaticFit.isPresent() || somaticCount == 0 || (Doubles.lessOrEqual(lowestScoringPurity, LOWEST_SCORE_MIN_PURITY) && Doubles
                .lessOrEqual(somaticFit.get().purity(), LOWEST_SCORE_MIN_PURITY));
    }

    private boolean isHighlyDiploid(@NotNull final FittedPurityScore score) {
        return Doubles.greaterOrEqual(score.maxDiploidProportion(), HIGHLY_DIPLOID_PERCENTAGE);
    }

    public FittedPurityStatus status() {
        return status;
    }

    public FittedPurity bestFit() {
        return bestFit;
    }

    public FittedPurityScore score() {
        return score;
    }

    @NotNull
    private static List<FittedPurity> candidates(double lowestScore, @NotNull final List<FittedPurity> purities) {
        return purities.stream().filter(inRange(lowestScore)).collect(Collectors.toList());
    }

    @NotNull
    private static Predicate<FittedPurity> inRange(final double score) {
        return fittedPurity -> {
            double absDifference = Math.abs(fittedPurity.score() - score);
            double relDifference = Math.abs(absDifference / score);
            return lessOrEqual(absDifference, ABS_RANGE) || lessOrEqual(relDifference, PERCENT_RANGE);
        };
    }
}
