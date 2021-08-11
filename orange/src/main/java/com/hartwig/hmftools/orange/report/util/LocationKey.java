package com.hartwig.hmftools.orange.report.util;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public class LocationKey {

    @NotNull
    private final String chromosome;
    @NotNull
    private final String chromosomeBand;

    public LocationKey(@NotNull final String chromosome, @NotNull final String chromosomeBand) {
        this.chromosome = chromosome;
        this.chromosomeBand = chromosomeBand;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final LocationKey that = (LocationKey) o;
        return chromosome.equals(that.chromosome) && chromosomeBand.equals(that.chromosomeBand);
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromosome, chromosomeBand);
    }
}
