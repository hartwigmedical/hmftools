package com.hartwig.hmftools.orange.conversion;

import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ConversionUtil
{
    @NotNull
    public static <T, R> Iterable<R> mapToIterable(@Nullable Collection<T> collection, @NotNull Function<T, R> mapper)
    {
        Collection<T> nonNull = Objects.requireNonNullElseGet(collection, List::of);
        return () -> nonNull.stream().map(mapper).iterator();
    }

    @NotNull
    public static <T, R> List<R> mapToList(@Nullable Collection<T> collection, @NotNull Function<T, R> mapper)
    {
        Collection<T> nonNull = Objects.requireNonNullElseGet(collection, List::of);
        return nonNull.stream().map(mapper).collect(Collectors.toList());
    }
}