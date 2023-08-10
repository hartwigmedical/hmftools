package com.hartwig.hmftools.amber;

import java.util.function.Predicate;

public class NormalHomozygousFilter implements Predicate<BaseDepth>
{
    @Override
    public boolean test(final BaseDepth bafEvidence)
    {
        return bafEvidence.isValid() && bafEvidence.AltSupport == 0;
    }
}
