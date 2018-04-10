package com.hartwig.hmftools.common.cobalt;

import static java.util.stream.Collectors.toList;

import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public class ReferenceRatioStatisticsFactory {

    public static ReferenceRatioStatistics fromReferenceRatio(@NotNull final Multimap<String, ReadRatio> readRatios) {
        return fromRatio(readRatios, ReadRatio::ratio);
    }

    public static ReferenceRatioStatistics fromCobalt(@NotNull final Multimap<String, CobaltRatio> readRatios) {
        return fromRatio(readRatios, CobaltRatio::referenceGCRatio);
    }

    @NotNull
    @VisibleForTesting
    static <T> ReferenceRatioStatistics fromRatio(@NotNull final Multimap<String, T> readRatios, @NotNull Function<T, Double> transform) {

        final List<Double> xRatios = readRatios.get("X").stream().map(transform).filter(x -> !Doubles.equal(x, -1)).collect(toList());
        int xCount = xRatios.size();
        double xMedian = xCount > 0 ? median(xRatios) : 0;

        final List<Double> yRatios = readRatios.get("Y").stream().map(transform).filter(x -> !Doubles.equal(x, -1)).collect(toList());
        int yCount = yRatios.size();
        double yMedian = yCount > 0 ? median(yRatios) : 0;

        return ImmutableReferenceRatioStatistics.builder()
                .xCount(xCount)
                .xMedian(xMedian)
                .yCount(yCount)
                .yMedian(yMedian)
                .build();
    }

    private static double median(@NotNull List<Double> ratios) {
        Collections.sort(ratios);
        int count = ratios.size();
        return ratios.size() % 2 == 0 ? (ratios.get(count / 2) + ratios.get(count / 2 - 1)) / 2 : ratios.get(count / 2);
    }
}
