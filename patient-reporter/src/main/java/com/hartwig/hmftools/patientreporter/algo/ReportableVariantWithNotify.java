package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.protect.purple.ReportableVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ReportableVariantWithNotify {

    @NotNull
    public abstract ReportableVariant variant();

    public abstract boolean notifyVariant();
}
