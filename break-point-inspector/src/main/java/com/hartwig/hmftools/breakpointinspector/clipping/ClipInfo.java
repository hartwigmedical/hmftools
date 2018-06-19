package com.hartwig.hmftools.breakpointinspector.clipping;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClipInfo {

    @NotNull
    public abstract SAMRecord record();

    @NotNull
    public abstract Location alignment();

    public abstract int length();

    @NotNull
    public abstract String sequence();

    public abstract boolean hardClipped();
    public abstract boolean left();
    public abstract boolean right();
}
