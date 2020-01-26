package com.hartwig.hmftools.protect.common;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableVariantAnalysis {
    @NotNull
    public abstract List<ReportableVariant> variantsToReport();
}
