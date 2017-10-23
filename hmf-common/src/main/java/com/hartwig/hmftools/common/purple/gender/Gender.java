package com.hartwig.hmftools.common.purple.gender;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.baf.TumorBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public enum Gender {
    MALE,
    FEMALE;

    private static final int MIN_BAF_COUNT = 1000;

    @NotNull
    public static Gender fromAmber(@NotNull final Multimap<String, TumorBAF> bafs) {
        return bafs.get("X").stream().filter(x -> x.position() > 2_699_520 && x.position() < 155_260_560).count() > MIN_BAF_COUNT ? FEMALE : MALE;
    }

    @NotNull
    public static Gender fromReferenceReadRatios(@NotNull final Multimap<String, ReadRatio> readRatios) {
        return fromRatio(readRatios, ReadRatio::ratio);
    }

    @NotNull
    public static Gender fromCobalt(@NotNull final Multimap<String, CobaltRatio> readRatios) {
        return fromRatio(readRatios, CobaltRatio::referenceGCRatio);
    }

    @NotNull
    private static <T> Gender fromRatio(@NotNull final Multimap<String, T> readRatios, @NotNull Function<T, Double> transform) {
        return Doubles.greaterThan(median(readRatios.get("X"), transform), 0.75) ? Gender.FEMALE : Gender.MALE;
    }

    private static <T> double median(@NotNull Collection<T> readRatios, @NotNull Function<T, Double> transform) {
        return median(readRatios.stream().map(transform).filter(x -> !Doubles.equal(x, -1)).collect(Collectors.toList()));
    }

    private static double median(@NotNull List<Double> ratios) {
        Collections.sort(ratios);
        int count = ratios.size();
        return ratios.size() % 2 == 0 ? (ratios.get(count / 2) + ratios.get(count / 2 - 1)) / 2 : ratios.get(count / 2);
    }
}
