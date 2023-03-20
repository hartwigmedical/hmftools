package com.hartwig.hmftools.orange.conversion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collectors;

public class ConversionUtil {

    private ConversionUtil() {
    }

    @NotNull
    public static <T, R> Iterable<R> mapToIterable(@Nullable Collection<T> collection, Function<T, R> mapper) {
        var nonNull = Objects.requireNonNullElseGet(collection, List::<T>of);
        return () -> nonNull.stream().map(mapper).iterator();
    }

    @NotNull
    public static <T, R> List<R> mapToList(@Nullable Collection<T> collection, Function<T, R> mapper) {
        var nonNull = Objects.requireNonNullElseGet(collection, List::<T>of);
        return nonNull.stream().map(mapper).collect(Collectors.toList());
    }
}
