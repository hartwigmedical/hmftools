package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import org.jetbrains.annotations.NotNull;

public enum RefVersion {
    HG19("hg19");

    @NotNull
    private final String refVersionString;

    RefVersion(@NotNull final String refVersionString) {
        this.refVersionString = refVersionString;
    }

    @NotNull
    public String refVersionString() {
        return refVersionString;
    }
}
