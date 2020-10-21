package com.hartwig.hmftools.protect.reportable;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.protect.variants.ReportableVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ReportableData {

    @NotNull
    List<ReportableVariant> variants();

    @NotNull
    List<ReportableGainLoss> copyNumberAlterations();
}
