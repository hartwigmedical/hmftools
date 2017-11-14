package com.hartwig.hmftools.patientdb.curators;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SearchToken {
    @NotNull
    abstract String term();

    abstract int startOffset();

    abstract int endOffset();

    @Value.Derived
    int length() {
        return term().length();
    }
}
