package com.hartwig.hmftools.common.drivercatalog;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DriverCatalogKey {

    @Nullable
    private final String gene;
    @NotNull
    private final String transcript;

    @NotNull
    public static DriverCatalogKey create(@NotNull String gene, @NotNull String transcript) {
        return new DriverCatalogKey(gene, transcript);
    }

    private DriverCatalogKey(@Nullable final String gene, @NotNull final String transcript) {
        this.gene = gene;
        this.transcript = transcript;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final DriverCatalogKey that = (DriverCatalogKey) o;
        return Objects.equals(gene, that.gene) && transcript.equals(that.transcript);
    }

    @Override
    public int hashCode() {
        return Objects.hash(gene, transcript);
    }

    @Override
    public String toString() {
        return String.format("gene(%s) transcript(%s)", gene, transcript);
    }
}