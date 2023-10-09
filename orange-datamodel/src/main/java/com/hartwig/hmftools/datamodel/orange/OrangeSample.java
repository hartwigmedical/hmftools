package com.hartwig.hmftools.datamodel.orange;

import com.hartwig.hmftools.datamodel.flagstat.Flagstat;
import com.hartwig.hmftools.datamodel.metrics.WGSMetrics;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeSample
{
    @NotNull
    WGSMetrics metrics();

    @NotNull
    Flagstat flagstat();
}
