package com.hartwig.hmftools.patientreporter.algo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LohGenesReporting {

    @NotNull
    public abstract String location();

    @NotNull
    public abstract String gene();

    @Nullable
    public abstract Long minorAlleleCopies();

    @Nullable
    public abstract Long tumorCopies();

}
