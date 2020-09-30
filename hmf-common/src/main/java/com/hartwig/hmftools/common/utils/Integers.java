package com.hartwig.hmftools.common.utils;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public final class Integers {

    private Integers() {
    }

    public static int medianPositiveValue(@NotNull final List<Integer> values) {
        return median(values, x -> x > 0);
    }

    public static int median(@NotNull final List<Integer> values, Predicate<Integer> filter) {
        final List<Integer> sortedFiltered = values.stream().filter(filter).sorted().collect(Collectors.toList());
        int count = sortedFiltered.size();
        if (count == 0) {
            return 0;
        }

        return count % 2 == 0 ? (sortedFiltered.get(count / 2) + sortedFiltered.get(count / 2 - 1)) / 2 : sortedFiltered.get(count / 2);
    }
}



