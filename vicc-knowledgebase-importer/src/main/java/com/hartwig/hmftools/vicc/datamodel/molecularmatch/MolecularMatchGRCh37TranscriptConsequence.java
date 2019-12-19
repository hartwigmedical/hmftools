package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchGRCh37TranscriptConsequence {

    @NotNull
    public abstract String transcript();

    @Nullable
    public abstract String cdna();

    @Nullable
    public abstract String aminoAcidChange();

    @NotNull
    public abstract List<String> txSites();

    @Nullable
    public abstract String intronNumber();

    @NotNull
    public abstract List<String> exonNumbers();
}
