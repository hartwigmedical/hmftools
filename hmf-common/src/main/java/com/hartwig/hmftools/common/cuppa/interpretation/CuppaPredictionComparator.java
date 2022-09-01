package com.hartwig.hmftools.common.cuppa.interpretation;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

public class CuppaPredictionComparator implements Comparator<CuppaPrediction> {

    @Override
    public int compare(@NotNull CuppaPrediction prediction1, @NotNull CuppaPrediction prediction2) {
        int likelihoodCompare = Double.compare(prediction2.likelihood(), prediction1.likelihood());
        if (likelihoodCompare != 0) {
            return likelihoodCompare;
        }

        return prediction1.cancerType().compareTo(prediction2.cancerType());
    }
}

