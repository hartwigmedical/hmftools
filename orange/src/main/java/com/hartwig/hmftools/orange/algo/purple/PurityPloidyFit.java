package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.PurpleQC;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurityPloidyFit {

    @NotNull
    public abstract PurpleQC qc();

    public abstract boolean hasSufficientQuality();

    @NotNull
    public abstract FittedPurityMethod fittedPurityMethod();

    public abstract boolean containsTumorCells();

    public abstract double purity();

    public abstract double minPurity();

    public abstract double maxPurity();

    public abstract double ploidy();

    public abstract double minPloidy();

    public abstract double maxPloidy();
}
