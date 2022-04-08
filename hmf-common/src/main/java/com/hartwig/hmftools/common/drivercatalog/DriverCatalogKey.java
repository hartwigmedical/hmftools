package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DriverCatalogKey {

    @Nullable
    private final String gene;
    @NotNull
    private final String transcript;

    private DriverCatalogKey(@Nullable final String gene, @NotNull final String transcript) {
        this.gene = gene;
        this.transcript = transcript;
    }

    @NotNull
    public static Set<DriverCatalogKey> buildUniqueKeysSet(@NotNull Map<DriverCatalogKey, DriverCatalog> geneDriverMap) {
        Set<DriverCatalogKey> keys = Sets.newHashSet();
        for (Map.Entry<DriverCatalogKey, DriverCatalog> entry : geneDriverMap.entrySet()) {
            keys.add(create(entry.getKey().gene, entry.getKey().transcript));
        }
        return keys;
    }

    @NotNull
    public static Set<DriverCatalogKey> buildUniqueKeysSet(@NotNull List<DriverCatalog> drivers) {
        Set<DriverCatalogKey> keys = Sets.newHashSet();
        for (DriverCatalog driver : drivers) {
            keys.add(create(driver.gene(), driver.transcript()));
        }
        return keys;
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