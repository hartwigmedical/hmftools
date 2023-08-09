package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.Collection;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.Integers;

import org.jetbrains.annotations.NotNull;

public class BaseDepthFilter implements Predicate<BaseDepth>
{
    private final int mMinDepth;
    private final int mMaxDepth;

    public BaseDepthFilter(
            final double minDepthPercentage, final double maxDepthPercentage, final Multimap<Chromosome, BaseDepth> evidence)
    {
        this(minDepthPercentage, maxDepthPercentage, evidence.values());
    }

    public BaseDepthFilter(final double minDepthPercentage, final double maxDepthPercentage, final Collection<BaseDepth> evidence)
    {
        int medianDepth = medianDepth(evidence);
        mMinDepth = (int) Math.round(medianDepth * minDepthPercentage);
        mMaxDepth = (int) Math.round(medianDepth * maxDepthPercentage);
        AMB_LOGGER.info("Median normal depth is {} reads: filtering reads outside of {} and {}", medianDepth, mMinDepth, mMaxDepth);
    }

    @Override
    public boolean test(final BaseDepth bafEvidence)
    {
        return bafEvidence.ReadDepth > 0 && bafEvidence.ReadDepth >= mMinDepth && bafEvidence.ReadDepth <= mMaxDepth;
    }

    private int medianDepth(final Collection<BaseDepth> evidence)
    {
        return Integers.medianPositiveValue(evidence.stream().map(x -> x.ReadDepth).collect(Collectors.toList()));
    }
}
