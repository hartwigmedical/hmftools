package com.hartwig.hmftools.datamodel.orange;

import com.hartwig.hmftools.datamodel.flagstat.Flagstat;
import com.hartwig.hmftools.datamodel.metrics.WGSMetrics;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OrangeSample {

    @NotNull
    public abstract WGSMetrics metrics();

    @NotNull
    public abstract Flagstat flagstat();
}
