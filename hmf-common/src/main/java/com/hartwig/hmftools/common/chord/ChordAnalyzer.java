package com.hartwig.hmftools.common.chord;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ChordAnalyzer {

    @NotNull
    public abstract ChordAnalysis chordAnalysis();
    @NotNull
    public abstract ChordStatus hrdStatus();
}
