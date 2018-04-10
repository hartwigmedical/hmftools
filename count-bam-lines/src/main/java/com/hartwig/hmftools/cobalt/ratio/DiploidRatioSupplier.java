package com.hartwig.hmftools.cobalt.ratio;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.cobalt.ReferenceRatioStatistics;
import com.hartwig.hmftools.common.cobalt.ReferenceRatioStatisticsFactory;

import org.jetbrains.annotations.NotNull;

class DiploidRatioSupplier {

    private static final long ROLLING_MEDIAN_MAX_DISTANCE = 5_000;
    private static final long ROLLING_MEDIAN_MIN_COVERAGE = 1_000;

    private final ListMultimap<String, ReadRatio> result = ArrayListMultimap.create();

    DiploidRatioSupplier(@NotNull final ListMultimap<String, ReadRatio> normalRatios) {

        final ReferenceRatioStatistics stats = ReferenceRatioStatisticsFactory.fromReferenceRatio(normalRatios);

        for (String chromosome : normalRatios.keySet()) {
            final List<ReadRatio> ratios = normalRatios.get(chromosome);
            final List<ReadRatio> adjustedRatios;
            if (chromosome.equals("Y")) {
                adjustedRatios = ratios;
            } else {
                double expectedRatio = chromosome.equals("X") && !stats.containsTwoXChromosomes() ? 0.5 : 1;
                adjustedRatios = new DiploidRatioNormalization(expectedRatio,
                        ROLLING_MEDIAN_MAX_DISTANCE,
                        ROLLING_MEDIAN_MIN_COVERAGE,
                        ratios).get();
            }
            result.replaceValues(chromosome, adjustedRatios);
        }
    }

    @NotNull
    ListMultimap<String, ReadRatio> result() {
        return result;
    }
}
