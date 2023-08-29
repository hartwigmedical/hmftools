package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.Collection;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.Integers;

public class BaseDepthFilter implements Predicate<PositionEvidence>
{
    private final int mMinDepth;
    private final int mMaxDepth;

    public BaseDepthFilter(
            final double minDepthPercentage, final double maxDepthPercentage, final Multimap<Chromosome, PositionEvidence> evidence)
    {
        this(minDepthPercentage, maxDepthPercentage, evidence.values());
    }

    public BaseDepthFilter(double minDepthPercentage, double maxDepthPercentage, final Collection<PositionEvidence> evidence)
    {
        if(minDepthPercentage == 0 && maxDepthPercentage == 0)
        {
            mMaxDepth = 0;
            mMinDepth = 0;
            return;
        }

        int medianDepth = medianDepth(evidence);
        mMinDepth = (int) Math.round(medianDepth * minDepthPercentage);
        mMaxDepth = (int) Math.round(medianDepth * maxDepthPercentage);

        AMB_LOGGER.info("median normal depth({}) reads, filtered(min={} max={})", medianDepth, mMinDepth, mMaxDepth);
    }

    @Override
    public boolean test(final PositionEvidence bafEvidence)
    {
        if(mMinDepth == 0 && mMaxDepth == 0)
            return true;

        return bafEvidence.ReadDepth > 0 && bafEvidence.ReadDepth >= mMinDepth && bafEvidence.ReadDepth <= mMaxDepth;
    }

    private int medianDepth(final Collection<PositionEvidence> evidence)
    {
        return Integers.medianPositiveValue(evidence.stream().map(x -> x.ReadDepth).collect(Collectors.toList()));
    }
}
