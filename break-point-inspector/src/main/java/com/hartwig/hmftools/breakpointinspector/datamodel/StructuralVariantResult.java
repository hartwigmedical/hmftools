package com.hartwig.hmftools.breakpointinspector.datamodel;

import java.util.Collection;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;

public class StructuralVariantResult {

    @NotNull
    private final SampleStats refStats;
    @NotNull
    private final SampleStats tumorStats;
    @NotNull
    private final Pair<Location, Location> breakpoints;
    @NotNull
    private final Pair<Double, Double> alleleFrequency;
    @NotNull
    private final Collection<String> filters;
    @NotNull
    private final QueryInterval[] queryIntervals;

    @NotNull
    public static StructuralVariantResult buildFromPartialData(@NotNull Pair<Location, Location> breakPoints,
            @NotNull final Collection<String> filters, @NotNull final QueryInterval[] queryIntervals) {
        // KODU: This is fairly ugly, but needed to prevent NPE.
        return new StructuralVariantResult(SampleStats.emptyStats(),
                SampleStats.emptyStats(),
                breakPoints,
                Pair.of(0D, 0D),
                filters,
                queryIntervals);
    }

    public StructuralVariantResult(@NotNull final SampleStats refStats, @NotNull final SampleStats tumorStats,
            @NotNull final Pair<Location, Location> breakpoints, @NotNull final Pair<Double, Double> alleleFrequency,
            @NotNull final Collection<String> filters, @NotNull final QueryInterval[] queryIntervals) {
        this.refStats = refStats;
        this.tumorStats = tumorStats;
        this.breakpoints = breakpoints;
        this.alleleFrequency = alleleFrequency;
        this.filters = filters;
        this.queryIntervals = queryIntervals;
    }

    @NotNull
    public SampleStats refStats() {
        return refStats;
    }

    @NotNull
    public SampleStats tumorStats() {
        return tumorStats;
    }

    @NotNull
    public Pair<Location, Location> breakpoints() {
        return breakpoints;
    }

    @NotNull
    public Pair<Double, Double> alleleFrequency() {
        return alleleFrequency;
    }

    @NotNull
    public Collection<String> filters() {
        return filters;
    }

    @NotNull
    public QueryInterval[] queryIntervals() {
        return queryIntervals;
    }

    @NotNull
    public String filterString() {
        return filters.isEmpty() ? "PASS" : String.join(";", filters);
    }

}
