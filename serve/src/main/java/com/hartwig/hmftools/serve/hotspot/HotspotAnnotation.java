package com.hartwig.hmftools.serve.hotspot;

import java.util.Set;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HotspotAnnotation {

    @NotNull
    private final Set<String> sources;
    @NotNull
    private final String gene;
    @Nullable
    private final String transcript;
    @NotNull
    private final String proteinAnnotation;

    public HotspotAnnotation(@NotNull final Set<String> sources, @NotNull final String gene, @Nullable final String transcript,
            @NotNull final String proteinAnnotation) {
        this.sources = sources;
        this.gene = gene;
        this.transcript = transcript;
        this.proteinAnnotation = proteinAnnotation;
    }

    @NotNull
    public Set<String> sources() {
        return sources;
    }

    @NotNull
    public String gene() {
        return gene;
    }

    @Nullable
    public String transcript() {
        return transcript;
    }

    @NotNull
    public String proteinAnnotation() {
        return proteinAnnotation;
    }

    @Override
    public String toString() {
        return "HotspotAnnotation{" + "sources=" + sources + ", gene='" + gene + '\'' + ", transcript='" + transcript + '\''
                + ", proteinAnnotation='" + proteinAnnotation + '\'' + '}';
    }
}
