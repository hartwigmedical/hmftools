package com.hartwig.hmftools.amber;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.NormalHeterozygousCheck;

public class NormalHeterozygousFilter implements Predicate<BaseDepth>
{
    private final NormalHeterozygousCheck mCheck;

    public NormalHeterozygousFilter(final double minHetAFPercentage, final double maxHetAFPercentage)
    {
        mCheck = new NormalHeterozygousCheck(minHetAFPercentage, maxHetAFPercentage);
    }

    @Override
    public boolean test(final BaseDepth bafEvidence)
    {
        return mCheck.test(bafEvidence.ReadDepth, bafEvidence.RefSupport, bafEvidence.AltSupport, bafEvidence.IndelCount);
    }
}
