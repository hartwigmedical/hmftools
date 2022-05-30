package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaPredictionComparatorTest {

    @Test
    public void canSortCuppaPredictions() {
        CuppaPrediction prediction1 = create("cancer C", 0.4);
        CuppaPrediction prediction2 = create("cancer B", 0.3);
        CuppaPrediction prediction3 = create("cancer A", 0.3);

        List<CuppaPrediction> predictions = Lists.newArrayList(prediction1, prediction2, prediction3);
        predictions.sort(new CuppaPredictionComparator());

        assertEquals(prediction1, predictions.get(0));
        assertEquals(prediction3, predictions.get(1));
        assertEquals(prediction2, predictions.get(2));
    }

    @NotNull
    private static CuppaPrediction create(@NotNull String cancerType, double likelihood) {
        return ImmutableCuppaPrediction.builder().cancerType(cancerType).likelihood(likelihood).build();
    }
}