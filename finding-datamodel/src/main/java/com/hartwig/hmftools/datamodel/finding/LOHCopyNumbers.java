package com.hartwig.hmftools.datamodel.finding;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LOHCopyNumbers extends Finding {

    @NotNull
    String location();

    @NotNull
    String gene();

    @Nullable
    Integer tumorCopies();

    @Nullable
    Integer tumorMinorAlleleCopies();
}
