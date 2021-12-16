package com.hartwig.hmftools.serve.sources.actin.curation;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public class CurationEntry {

    @NotNull
    private final String gene;
    @NotNull
    private final String name;

    public CurationEntry(@NotNull final String gene, @NotNull final String name) {
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
        final CurationEntry
                that = (CurationEntry) o;
        return gene.equals(that.gene) && name.equals(that.name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(gene, name);
    }

    @Override
    public String toString() {
        return "CurationEntry{" + "gene='" + gene + '\'' + ", name='" + name + '\'' + '}';
    }
}
