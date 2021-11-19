package com.hartwig.hmftools.common.chord;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ChordTestFactory {

    private ChordTestFactory() {
    }

    @NotNull
    public static ChordAnalysis createMinimalTestChordAnalysis() {
        return ImmutableChordAnalysis.builder()
                .BRCA1Value(0D)
                .BRCA2Value(0D)
                .hrdValue(0D)
                .hrStatus(ChordStatus.HR_PROFICIENT)
                .hrdType(Strings.EMPTY)
                .remarksHrStatus(Strings.EMPTY)
                .remarksHrdType(Strings.EMPTY)
                .build();
    }
}
