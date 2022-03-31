package com.hartwig.hmftools.common.drivercatalog;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DriverCatalogKey {

    @Nullable
    private final String gene;
    @NotNull
    private final String transript;

    private DriverCatalogKey(@Nullable final String gene, @NotNull final String transript) {
        this.gene = gene;
        this.transript = transript;
    }

    @NotNull
    public static DriverCatalogKey create(@NotNull String gene, @NotNull String transcript) {
        return new DriverCatalogKey(gene, transcript);
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
        return Objects.equals(gene, that.gene) && transript.equals(that.transript);
    }

    @Override
    public int hashCode() {
        return Objects.hash(gene, transript);
    }

    @Override
    public String toString() {
        return "DriverCatalogKey{" + "gene='" + gene + '\'' + ", transript='" + transript + '\'' + '}';
    }
}