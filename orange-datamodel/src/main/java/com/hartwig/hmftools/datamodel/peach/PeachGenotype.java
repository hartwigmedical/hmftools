package com.hartwig.hmftools.datamodel.peach;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PeachGenotype
{
    @NotNull
    String gene();

    @NotNull
    String haplotype();

    @NotNull
    String function();

    @NotNull
    String linkedDrugs();

    @NotNull
    String urlPrescriptionInfo();

    @Deprecated
    @Nullable
    String panelVersion();

    @Deprecated
    @Nullable
    String repoVersion();
}