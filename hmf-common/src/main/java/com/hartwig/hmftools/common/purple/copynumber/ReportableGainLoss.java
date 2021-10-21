package com.hartwig.hmftools.common.purple.copynumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGainLoss {

    @NotNull
    @Value.Derived
    public String genomicEvent() {
        return this.gene() + " " + this.interpretation().display();
    }

    @NotNull
    public abstract CopyNumberInterpretation interpretation();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String chromosomeBand();

    @NotNull
    public abstract String gene();

    public abstract long minCopies();

    public abstract long maxCopies();
}
