package com.hartwig.hmftools.cobalt.ratio;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

class DiploidRatioSupplier {

    private static final long ROLLING_MEDIAN_MAX_DISTANCE = 5_000;
    private static final long ROLLING_MEDIAN_MIN_COVERAGE = 1_000;

    private final ListMultimap<String, ReadRatio> result = ArrayListMultimap.create();

    DiploidRatioSupplier(@NotNull final Gender gender, @NotNull final ListMultimap<String, ReadRatio> normalRatios) {
        for (String chromosome : normalRatios.keySet()) {
            double expectedRatio = HumanChromosome.fromString(chromosome).isHomologous(gender) ? 1 : 0.5;
            final List<ReadRatio> ratios = normalRatios.get(chromosome);
            final List<ReadRatio> adjustedRatios =
                    new DiploidRatioNormalization(expectedRatio, ROLLING_MEDIAN_MAX_DISTANCE, ROLLING_MEDIAN_MIN_COVERAGE, ratios).get();
            result.replaceValues(chromosome, adjustedRatios);
        }
    }

    @NotNull
    ListMultimap<String, ReadRatio> result() {
        return result;
    }
}
