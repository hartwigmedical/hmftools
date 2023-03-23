package com.hartwig.hmftools.datamodel.isofox;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(allParameters = true, passAnnotations = { NotNull.class, Nullable.class })
public interface IsofoxRnaStatistics {

    long totalFragments();

    long duplicateFragments();

    @NotNull
    String qcStatus();
}
