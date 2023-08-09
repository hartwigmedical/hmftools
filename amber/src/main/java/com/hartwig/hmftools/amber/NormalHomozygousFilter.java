package com.hartwig.hmftools.amber;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.amber.BaseDepth;

public class NormalHomozygousFilter implements Predicate<BaseDepth>
{
    @Override
    public boolean test(final BaseDepth bafEvidence)
    {
        return bafEvidence.isValid() && bafEvidence.AltSupport == 0;
    }
}
