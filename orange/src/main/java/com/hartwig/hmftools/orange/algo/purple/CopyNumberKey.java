package com.hartwig.hmftools.orange.algo.purple;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

class CopyNumberKey {

    @NotNull
    private final String chromosome;
    @NotNull
    private final String chromosomeBand;

    public CopyNumberKey(@NotNull final String chromosome, @NotNull final String chromosomeBand) {
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
        final CopyNumberKey that = (CopyNumberKey) o;
        return chromosome.equals(that.chromosome) && chromosomeBand.equals(that.chromosomeBand);
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromosome, chromosomeBand);
    }
}
