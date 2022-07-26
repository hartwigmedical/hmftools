package com.hartwig.hmftools.orange.algo.cuppa;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CuppaInterpretation {

    private CuppaInterpretation() {
    }

    @NotNull
    public static CuppaPrediction best(@NotNull CuppaData cuppa) {
       CuppaPrediction best = rank(cuppa, 1);
       if (best == null) {
           throw new IllegalStateException("Could not determine best prediction from cuppa: " + cuppa);
       }
       return best;
    }

    @Nullable
    public static CuppaPrediction rank(@NotNull CuppaData cuppa, int rank) {
        List<CuppaPrediction> predictions = Lists.newArrayList();
        predictions.addAll(cuppa.predictions());
        predictions.sort(new CuppaPredictionComparator());

        return predictions.size() >= rank ? predictions.get(rank - 1) : null;
    }
}
