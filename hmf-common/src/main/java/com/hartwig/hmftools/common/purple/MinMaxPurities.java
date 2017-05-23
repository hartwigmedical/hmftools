package com.hartwig.hmftools.common.purple;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.numeric.Doubles;

class MinMaxPurities {

    static PurityBounds bounds(List<FittedPurity> purities) {
        ImmutablePurityBounds.Builder builder = ImmutablePurityBounds.builder().minPloidy(0).minPurity(0).maxPloidy(0).maxPurity(0);

        if (!purities.isEmpty()) {
            Collections.sort(purities);

            FittedPurity best = purities.get(0);

            List<FittedPurity> within10Percent = purities.stream().filter(inPercentRange(0.10, best.score())).collect(
                    Collectors.toList());
            within10Percent.stream().max(MinMaxPurities::comparePloidy).ifPresent(x -> builder.maxPloidy(x.ploidy()));
            within10Percent.stream().min(MinMaxPurities::comparePloidy).ifPresent(x -> builder.minPloidy(x.ploidy()));
            within10Percent.stream().max(MinMaxPurities::comparePurity).ifPresent(x -> builder.maxPurity(x.purity()));
            within10Percent.stream().min(MinMaxPurities::comparePurity).ifPresent(x -> builder.minPurity(x.purity()));

        }

        return builder.build();
    }

    static Map<String, FittedPurity> minMax(List<FittedPurity> purities) {
        Map<String, FittedPurity> results = new LinkedHashMap<>();

        if (!purities.isEmpty()) {
            Collections.sort(purities);

            FittedPurity best = purities.get(0);
            results.put("MinScore", best);

            List<FittedPurity> within10Percent = purities.stream().filter(inPercentRange(0.10, best.score())).collect(
                    Collectors.toList());
            within10Percent.stream().max(MinMaxPurities::comparePloidy).ifPresent(x -> results.put("MaxPloidy10%", x));
            within10Percent.stream().min(MinMaxPurities::comparePloidy).ifPresent(x -> results.put("MinPloidy10%", x));
            within10Percent.stream().max(MinMaxPurities::comparePurity).ifPresent(x -> results.put("MaxPurity10%", x));
            within10Percent.stream().min(MinMaxPurities::comparePurity).ifPresent(x -> results.put("MinPurity10%", x));


            List<FittedPurity> within20Percent = purities.stream().filter(inPercentRange(0.20, best.score())).collect(
                    Collectors.toList());

            within20Percent.stream().max(MinMaxPurities::comparePloidy).ifPresent(x -> results.put("MaxPloidy20%", x));
            within20Percent.stream().min(MinMaxPurities::comparePloidy).ifPresent(x -> results.put("MinPloidy20%", x));
            within20Percent.stream().max(MinMaxPurities::comparePurity).ifPresent(x -> results.put("MaxPurity20%", x));
            within20Percent.stream().min(MinMaxPurities::comparePurity).ifPresent(x -> results.put("MinPurity20%", x));
        }

        return results;
    }

    private static int comparePurity(FittedPurity o1, FittedPurity o2) {
        return Double.compare(o1.purity(), o2.purity());
    }

    private static int comparePloidy(FittedPurity o1, FittedPurity o2) {
        return Double.compare(o1.ploidy(), o2.ploidy());
    }

    private static Predicate<FittedPurity> inPercentRange(double percent, double score) {
        return fittedPurity -> Doubles.lessOrEqual(Math.abs((fittedPurity.score() - score) / score), percent);
    }

}
