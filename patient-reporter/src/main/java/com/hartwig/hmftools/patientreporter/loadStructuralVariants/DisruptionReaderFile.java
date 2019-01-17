package com.hartwig.hmftools.patientreporter.loadStructuralVariants;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DisruptionReaderFile {

    @NotNull
    public abstract Boolean reportable();

    @NotNull
    public abstract String svId();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String position();

    @NotNull
    public abstract String orientation();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String strand();

    @NotNull
    public abstract String regionType();

    @NotNull
    public abstract String codingType();

    @NotNull
    public abstract String biotype();

    @NotNull
    public abstract String exon();

    @NotNull
    public abstract Boolean isDisruptive();
}
