package com.hartwig.hmftools.common.hla;

import java.util.List;
import java.util.Map;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class HlaAllelesReportingData {

    @NotNull
    public abstract Map<String, List<HlaReporting>> hlaAllelesReporting();

    @NotNull
    public abstract String hlaQC();
}