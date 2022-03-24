package com.hartwig.hmftools.common.amber.qc;

import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class AmberQC {

    private static final double CONTAMINATION_MAX_FAIL = 0.1;

    public AmberQCStatus status() {

        if (Doubles.greaterThan(contamination(), CONTAMINATION_MAX_FAIL)) {
            return AmberQCStatus.FAIL;
        }

        if (Doubles.greaterThan(contamination(), 0)) {
            return AmberQCStatus.WARN;
        }

        return AmberQCStatus.PASS;
    }

    public abstract double contamination();

    public abstract double consanguinityProportion();

    @Nullable
    public abstract String uniparentalDisomy();
}
