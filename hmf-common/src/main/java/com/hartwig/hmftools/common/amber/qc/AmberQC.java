package com.hartwig.hmftools.common.amber.qc;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class AmberQC {

    private static final double MAX_FAIL = 0.51;
    private static final double MAX_WARN = 0.50;

    private static final double MIN_FAIL = 0.48;
    private static final double MIN_WARN = 0.487;

    public AmberQCStatus status() {

        if (Doubles.greaterThan(meanBAF(), MAX_FAIL) || Doubles.lessThan(meanBAF(), MIN_FAIL)) {
            return AmberQCStatus.FAIL;
        }

        if (Doubles.greaterThan(meanBAF(), MAX_WARN) || Doubles.lessThan(meanBAF(), MIN_WARN)) {
            return AmberQCStatus.WARN;
        }

        return AmberQCStatus.PASS;
    }

    public abstract double meanBAF();

}
