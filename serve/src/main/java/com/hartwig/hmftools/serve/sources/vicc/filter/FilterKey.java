package com.hartwig.hmftools.serve.sources.vicc.filter;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

class FilterKey {

    @NotNull
    private final ViccSource source;
    @NotNull
    private final String gene;
    @NotNull
    private final String name;

    public FilterKey(@NotNull final ViccSource source, @NotNull final String gene, @NotNull final String name) {
        this.source = source;
        this.gene = gene;
        this.name = name;
    }

    @NotNull
    @VisibleForTesting
    ViccSource source() {
        return source;
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
        return source == filterKey.source && gene.equals(filterKey.gene) && name.equals(filterKey.name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(source, gene, name);
    }

    @Override
    public String toString() {
        return "FilterKey{" + "source=" + source + ", gene='" + gene + '\'' + ", name='" + name + '\'' + '}';
    }
}
