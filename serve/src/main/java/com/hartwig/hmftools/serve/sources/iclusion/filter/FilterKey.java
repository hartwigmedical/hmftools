package com.hartwig.hmftools.serve.sources.iclusion.filter;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

class FilterKey {

    @NotNull
    private final String gene;
    @NotNull
    private final String name;

    public FilterKey(@NotNull final String gene, @NotNull final String name) {
        this.gene = gene;
        this.name = name;
    }

    @NotNull
    @VisibleForTesting
    String gene() {
        return gene;
    }

    @NotNull
    @VisibleForTesting
    String name() {
        return name;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final FilterKey filterKey = (FilterKey) o;
        return gene.equals(filterKey.gene) && name.equals(filterKey.name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(gene, name);
    }

    @Override
    public String toString() {
        return "FilterKey{" + "gene='" + gene + '\'' + ", name='" + name + '\'' + '}';
    }
}
