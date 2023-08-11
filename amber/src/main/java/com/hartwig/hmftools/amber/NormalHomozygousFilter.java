package com.hartwig.hmftools.amber;

import java.util.function.Predicate;

public class NormalHomozygousFilter implements Predicate<PositionEvidence>
{
    @Override
    public boolean test(final PositionEvidence bafEvidence)
    {
        return bafEvidence.isValid() && bafEvidence.AltSupport == 0;
    }
}
