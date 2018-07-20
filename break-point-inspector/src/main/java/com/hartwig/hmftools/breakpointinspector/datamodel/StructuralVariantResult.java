package com.hartwig.hmftools.breakpointinspector.datamodel;

import java.util.Collection;
import java.util.Collections;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.apache.commons.lang3.tuple.Pair;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantResult {

    @NotNull
    public abstract SampleStats refStats();

    @NotNull
    public abstract SampleStats tumorStats();

    @NotNull
    public abstract Pair<Location, Location> breakpoints();

    @NotNull
    public abstract Pair<Double, Double> alleleFrequency();

    @NotNull
    public abstract Collection<String> filters();

    @NotNull
    public abstract QueryInterval[] queryIntervals();

    @NotNull
    public static StructuralVariantResult buildForBreakpointError(@NotNull Pair<Location, Location> breakPoints,
            @NotNull final QueryInterval[] queryIntervals) {
        // KODU: This is fairly ugly, but needed to prevent NPE.
        return ImmutableStructuralVariantResult.builder()
                .refStats(SampleStats.emptyStats())
                .tumorStats(SampleStats.emptyStats())
                .breakpoints(breakPoints)
                .alleleFrequency(Pair.of(0D, 0D))
                .filters(Collections.singletonList(FilterType.BREAKPOINT_ERROR.toString()))
                .queryIntervals(queryIntervals)
                .build();
    }

    @NotNull
    @Value.Derived
    public String filterString() {
        return filters().isEmpty() ? "PASS" : String.join(";", filters());
    }
}
