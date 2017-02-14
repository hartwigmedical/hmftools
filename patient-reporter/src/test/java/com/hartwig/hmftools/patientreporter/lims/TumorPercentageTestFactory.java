package com.hartwig.hmftools.patientreporter.lims;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class TumorPercentageTestFactory {

    private TumorPercentageTestFactory() {
    }

    @NotNull
    public static TumorPercentages buildTestTumorPercentages() {
        return new TumorPercentages(Maps.newHashMap());
    }
}
