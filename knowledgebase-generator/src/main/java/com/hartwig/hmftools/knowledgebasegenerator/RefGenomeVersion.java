package com.hartwig.hmftools.knowledgebasegenerator;

import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion {
    HG19("hg19");

    @NotNull
    private final String refVersionString;

    RefGenomeVersion(@NotNull final String refVersionString) {
        this.refVersionString = refVersionString;
    }

    @NotNull
    public String refVersionString() {
        return refVersionString;
    }
}
