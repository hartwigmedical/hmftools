package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicCoordinates {

    @Nullable
    public abstract String chromosome2();

    @Nullable
    public abstract String referenceBases();

    @Nullable
    public abstract String start2();

    @Nullable
    public abstract String variantBases();

    @NotNull
    public abstract String stop();

    @Nullable
    public abstract String stop2();

    @Nullable
    public abstract String representativeTranscript2();

    @NotNull
    public abstract String start();

    @NotNull
    public abstract String representativeTranscript();

    @NotNull
    public abstract String ensemblVersion();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String referenceBuild();

}
