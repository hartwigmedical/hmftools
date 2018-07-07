package com.hartwig.hmftools.breakpointinspector;

import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public final class Util {

    private Util() {
    }

    @NotNull
    public static List<String> prefixList(final List<String> list, final String prefix) {
        return list.stream().map(s -> prefix + s).collect(Collectors.toList());
    }

    @NotNull
    public static <U> List<String> toStrings(final List<U> list) {
        return list.stream().map(Objects::toString).collect(Collectors.toList());
    }
}