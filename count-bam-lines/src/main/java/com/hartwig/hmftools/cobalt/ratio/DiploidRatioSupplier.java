package com.hartwig.hmftools.cobalt.ratio;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.ratio.ReadRatio;

import org.jetbrains.annotations.NotNull;

class DiploidRatioSupplier {

    private static final long ROLLING_MEDIAN_MAX_DISTANCE = 5_000;
    private static final long ROLLING_MEDIAN_MIN_COVERAGE = 1_000;

    private final ListMultimap<String, ReadRatio> result = ArrayListMultimap.create();

    DiploidRatioSupplier(@NotNull final ListMultimap<String, ReadRatio> normalRatios) {

        final Gender gender = determineGender(normalRatios.get("X"));

        for (String chromosome : normalRatios.keySet()) {

            double expectedRatio = HumanChromosome.fromString(chromosome).isHomologous(gender) ? 1 : 0.5;
            final List<ReadRatio> ratios = normalRatios.get(chromosome);
            final List<ReadRatio> adjustedRatios =
                    new DiploidRatioNormalization(expectedRatio, ROLLING_MEDIAN_MAX_DISTANCE, ROLLING_MEDIAN_MIN_COVERAGE, ratios).get();
            result.replaceValues(chromosome, adjustedRatios);
        }
    }

    ListMultimap<String, ReadRatio> result() {
        return result;
    }

    private static Gender determineGender(@NotNull final Collection<ReadRatio> ratios) {
        return Doubles.greaterThan(median(ratios), 0.75) ? Gender.FEMALE : Gender.MALE;
    }

    private static double median(Collection<ReadRatio> readRatios) {
        final List<Double> ratios = readRatios.stream().map(ReadRatio::ratio).collect(Collectors.toList());
        Collections.sort(ratios);
        int count = ratios.size();
        return ratios.size() % 2 == 0 ? (ratios.get(count / 2) + ratios.get(count / 2 - 1)) / 2 : ratios.get(count / 2);
    }

}
