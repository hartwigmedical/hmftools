package com.hartwig.hmftools.common.purple.purity;

import java.util.List;

public class BestFitFactory {

    private final FittedPurity bestFit;
    private final FittedPurityScore score;

    public BestFitFactory(final List<FittedPurity> fittedPurities) {

        final FittedPurity lowestScore = fittedPurities.get(0);
        score = FittedPurityScoreFactory.score(fittedPurities);
        bestFit = lowestScore;
    }

    public FittedPurityStatus status() {
        return FittedPurityStatus.NORMAL;
    }

    public FittedPurity bestFit() {
        return bestFit;
    }

    public FittedPurityScore score() {
        return score;
    }

}
