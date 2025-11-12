package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.linx.LinxFusion;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Fusion extends Driver, LinxFusion {
}
