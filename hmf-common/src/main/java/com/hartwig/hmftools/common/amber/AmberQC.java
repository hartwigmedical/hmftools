package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public record AmberQC(double contamination, double consanguinityProportion, String uniparentalDisomy)
{
    private static final double CONTAMINATION_MAX_FAIL = 0.1;

    public AmberQCStatus status()
    {
        if(Doubles.greaterThan(contamination(), CONTAMINATION_MAX_FAIL))
        {
            return AmberQCStatus.FAIL;
        }

        if(Doubles.greaterThan(contamination(), 0))
        {
            return AmberQCStatus.WARN;
        }

        return AmberQCStatus.PASS;
    }
}
