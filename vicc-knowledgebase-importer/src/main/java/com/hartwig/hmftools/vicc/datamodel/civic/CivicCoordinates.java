package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicCoordinates {

    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String start();

    @Nullable
    public abstract String stop();

    @Nullable
    public abstract String referenceBases();

    @Nullable
    public abstract String variantBases();

    @Nullable
    public abstract String representativeTranscript();

    @Nullable
    public abstract String ensemblVersion();

    @Nullable
    public abstract String referenceBuild();

    @Nullable
    public abstract String chromosome2();

    @Nullable
    public abstract String start2();

    @Nullable
    public abstract String stop2();

    @Nullable
    public abstract String representativeTranscript2();

}
