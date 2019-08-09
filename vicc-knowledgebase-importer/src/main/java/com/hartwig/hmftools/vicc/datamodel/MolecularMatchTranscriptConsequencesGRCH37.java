package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTranscriptConsequencesGRCH37 {

    @Nullable
    public abstract String aminoAcidChange();

    @NotNull
    public abstract List<String> txSites();

    @Nullable
    public abstract List<String> exonNumber();

    @Nullable
    public abstract String intronNumber();

    @NotNull
    public abstract String transcript();

    @Nullable
    public abstract String cdna();
}
