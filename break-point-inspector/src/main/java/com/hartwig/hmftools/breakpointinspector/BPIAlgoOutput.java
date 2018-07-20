package com.hartwig.hmftools.breakpointinspector;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;
import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BPIAlgoOutput {

    @NotNull
    public abstract List<VariantContext> variants();

    @NotNull
    public abstract List<String> tsvOutput();

    @NotNull
    public abstract QueryInterval[] optimizedIntervals();

}
