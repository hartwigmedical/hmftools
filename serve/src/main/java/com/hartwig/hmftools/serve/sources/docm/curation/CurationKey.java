package com.hartwig.hmftools.serve.sources.docm.curation;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

class CurationKey {

    @NotNull
    private final String gene;
    @NotNull
    private final String transcript;
    @NotNull
    private final String proteinAnnotation;

    public CurationKey(@NotNull final String gene, @NotNull final String transcript, @NotNull final String proteinAnnotation) {
        this.gene = gene;
        this.transcript = transcript;
        this.proteinAnnotation = proteinAnnotation;
    }

    @NotNull
    @VisibleForTesting
    String gene() {
        return gene;
    }

    @NotNull
    @VisibleForTesting
    public String transcript() {
        return transcript;
    }

    @NotNull
    @VisibleForTesting
    public String proteinAnnotation() {
        return proteinAnnotation;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final CurationKey that = (CurationKey) o;
        return gene.equals(that.gene) && transcript.equals(that.transcript) && proteinAnnotation.equals(that.proteinAnnotation);
    }

    @Override
    public int hashCode() {
        return Objects.hash(gene, transcript, proteinAnnotation);
    }

    @Override
    public String toString() {
        return "CurationKey{" + "gene='" + gene + '\'' + ", transcript='" + transcript + '\'' + ", proteinAnnotation='" + proteinAnnotation
                + '\'' + '}';
    }
}
