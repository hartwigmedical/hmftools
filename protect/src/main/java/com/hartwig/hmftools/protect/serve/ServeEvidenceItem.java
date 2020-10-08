package com.hartwig.hmftools.protect.serve;

import com.hartwig.hmftools.serve.sources.Source;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ServeEvidenceItem {

    // This should display the event which led to the evidence, eg "EGFR Amplification".
    @NotNull
    public abstract String event();

    @NotNull
    public abstract Source source();

    @NotNull
    public abstract String treatment();

    @NotNull
    public abstract String level();

    @NotNull
    public abstract String response();

    public abstract boolean isOnLabel();

    @NotNull
    public abstract String cancerType();

}
