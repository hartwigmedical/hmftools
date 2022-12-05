package com.hartwig.hmftools.orange.algo.cuppa;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CuppaInterpretationTest {

    @Test
    public void canFindPredictionsByRank() {
        CuppaPrediction cancerA = create("Cancer A", 0.3);
        CuppaPrediction cancerB = create("Cancer B", 0.4);
        CuppaPrediction cancerC = create("Cancer C", 0.1);
        CuppaPrediction cancerD = create("Cancer D", 0.2);

        CuppaData cuppa = withPredictions(Lists.newArrayList(cancerA, cancerB, cancerC, cancerD));

        assertEquals(cancerB, CuppaInterpretation.best(cuppa));
        assertEquals(cancerB, CuppaInterpretation.rank(cuppa, 1));
        assertEquals(cancerA, CuppaInterpretation.rank(cuppa, 2));
        assertEquals(cancerD, CuppaInterpretation.rank(cuppa, 3));
        assertEquals(cancerC, CuppaInterpretation.rank(cuppa, 4));
        assertNull(CuppaInterpretation.rank(cuppa, 5));
    }

    @NotNull
    private static CuppaData withPredictions(@NotNull List<CuppaPrediction> predictions) {
        return ImmutableCuppaData.builder().from(TestCuppaFactory.createMinimalCuppaData()).predictions(predictions).build();
    }

    @NotNull
    private static CuppaPrediction create(@NotNull String cancerType, double likelihood) {
        return ImmutableCuppaPrediction.builder().cancerType(cancerType).likelihood(likelihood).build();
    }
}