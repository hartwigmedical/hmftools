package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

class CurationEntry {

    @NotNull
    private final String geneSymbol;
    @NotNull
    private final String variant;

    public CurationEntry(@NotNull final String geneSymbol, @NotNull final String variant) {
        this.geneSymbol = geneSymbol;
        this.variant = variant;
    }

    @NotNull
    @VisibleForTesting
    String geneSymbol() {
        return geneSymbol;
    }

    @NotNull
    @VisibleForTesting
    String variant() {
        return variant;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final CurationEntry that = (CurationEntry) o;
        return geneSymbol.equals(that.geneSymbol) && variant.equals(that.variant);
    }

    @Override
    public int hashCode() {
        return Objects.hash(geneSymbol, variant);
    }

    @Override
    public String toString() {
        return "CurationEntry{" + "geneSymbol='" + geneSymbol + '\'' + ", variant='" + variant + '\'' + '}';
    }
}
