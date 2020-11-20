package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class CurationKey {

    @NotNull
    private final ViccSource source;
    @NotNull
    private final String gene;
    @Nullable
    private final String transcript;
    @NotNull
    private final String featureName;

    public CurationKey(@NotNull final ViccSource source, @NotNull final String gene, @Nullable final String transcript,
            @NotNull final String featureName) {
        this.source = source;
        this.gene = gene;
        this.transcript = transcript;
        this.featureName = featureName;
    }

    @VisibleForTesting
    @NotNull
    ViccSource source() {
        return source;
    }

    @VisibleForTesting
    @NotNull
    String gene() {
        return gene;
    }

    @VisibleForTesting
    @Nullable
    String transcript() {
        return transcript;
    }

    @VisibleForTesting
    @NotNull
    String featureName() {
        return featureName;
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
        return source == that.source && gene.equals(that.gene) && Objects.equals(transcript, that.transcript)
                && featureName.equals(that.featureName);
    }

    @Override
    public int hashCode() {
        return Objects.hash(source, gene, transcript, featureName);
    }

    @Override
    public String toString() {
        return "CurationKey{" + "source=" + source + ", gene='" + gene + '\'' + ", transcript='" + transcript + '\'' + ", featureName='"
                + featureName + '\'' + '}';
    }
}
