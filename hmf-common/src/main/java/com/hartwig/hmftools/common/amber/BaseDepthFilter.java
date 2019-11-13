package com.hartwig.hmftools.common.amber;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

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
        final List<Integer> reads =
                evidence.values().stream().map(BaseDepth::readDepth).filter(x -> x > 0).sorted().collect(Collectors.toList());
        int count = reads.size();
        return count % 2 == 0 ? (reads.get(count / 2) + reads.get(count / 2 - 1)) / 2 : reads.get(count / 2);
    }
}
