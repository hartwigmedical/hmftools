package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class BestFitFactory {

    private static final double PERCENT_RANGE = 0.1;
    private static final double ABS_RANGE = 0.0005;
    private static final double HIGHLY_DIPLOID_PERCENTAGE = 0.98;
    private static final double MIN_PURITY = 0.1;
    private static final int MIN_VARIANTS = 1000;

    private final FittedPurity bestFit;
    private final FittedPurityScore score;
    private final FittedPurityStatus status;

    public BestFitFactory(final List<FittedPurity> fittedPurities, Supplier<List<SomaticVariant>> somaticSupplier) {
        assert (!fittedPurities.isEmpty());

        Collections.sort(fittedPurities);
        FittedPurity lowestScore = fittedPurities.get(0);

        final List<FittedPurity> candidates = candidates(lowestScore.score(), fittedPurities);
        score = FittedPurityScoreFactory.score(candidates);

        if (Doubles.lessOrEqual(score.minPurity(), MIN_PURITY) && isHighlyDiploid(score)) {
            final List<SomaticVariant> somatics = somaticSupplier.get();
            if (somatics.size() == 0) {
                status = FittedPurityStatus.HIGHLY_DIPLOID;
                bestFit = lowestScore;
            } else if (somatics.size() < MIN_VARIANTS) {
                status = FittedPurityStatus.NO_TUMOR;
                bestFit = lowestScore;
            } else {
                status = FittedPurityStatus.SOMATIC;
                bestFit = BestFitSomatics.bestSomaticFit(somatics, candidates);
            }

        } else {
            status = FittedPurityStatus.NORMAL;
            bestFit = lowestScore;
        }
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
