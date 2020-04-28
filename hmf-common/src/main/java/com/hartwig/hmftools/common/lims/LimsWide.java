package com.hartwig.hmftools.common.lims;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsWide {

    @NotNull
    public abstract String studyName();

    @NotNull
    public abstract String ReportReceiverName();

    @NotNull
    public abstract String ReportReceiverEmail();

}
