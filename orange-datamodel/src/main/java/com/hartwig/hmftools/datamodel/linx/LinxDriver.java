package com.hartwig.hmftools.datamodel.linx;

import com.hartwig.hmftools.datamodel.driver.ReportedStatus;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxDriver
{
    @NotNull
    String gene();

    @NotNull
    LinxDriverType type();
}
