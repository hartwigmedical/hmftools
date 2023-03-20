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
    public static <T, R> Iterable<R> convertCollection(@Nullable Collection<T> collection, Function<T, R> conversion) {
        var nonNull = Objects.requireNonNullElseGet(collection, List::<T>of);
        return () -> nonNull.stream().map(conversion).iterator();
    }

    @NotNull
    public static <T, R> List<R> mapCollection(@Nullable Collection<T> collection, Function<T, R> conversion) {
        var nonNull = Objects.requireNonNullElseGet(collection, List::<T>of);
        return nonNull.stream().map(conversion).collect(Collectors.toList());
    }

}
