package com.hartwig.hmftools.common.amber;

import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.Integers;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BaseDepthFilter implements Predicate<BaseDepth> {

    private static final Logger LOGGER = LogManager.getLogger(BaseDepthFilter.class);

    private final int minDepth;
    private final int maxDepth;

    public BaseDepthFilter(final double minDepthPercentage, final double maxDepthPercentage,
            @NotNull final Multimap<Chromosome, BaseDepth> evidence) {
        int medianDepth = medianDepth(evidence);
        minDepth = (int) Math.round(medianDepth * minDepthPercentage);
        maxDepth = (int) Math.round(medianDepth * maxDepthPercentage);
        LOGGER.info("Median normal depth is {} reads: filtering reads outside of {} and {}", medianDepth, minDepth, maxDepth);
    }

    @Override
    public boolean test(final BaseDepth bafEvidence) {
        return bafEvidence.readDepth() > 0 && bafEvidence.readDepth() >= minDepth && bafEvidence.readDepth() <= maxDepth;
    }

    private int medianDepth(@NotNull final Multimap<Chromosome, BaseDepth> evidence) {
        return Integers.medianPositiveValue(evidence.values().stream().map(BaseDepth::readDepth).collect(Collectors.toList()));
    }
}
