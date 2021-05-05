package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.protect.purple.ReportableVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableVariantNotify {

    @NotNull
    public abstract ReportableVariant reportableVariant();

    public abstract boolean notifyVariant();
}
