package com.hartwig.hmftools.common.teal;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface TelomereLength
{
    String type();
    double rawTelomereLength();
    double finalTelomereLength();
    int fullFragments();
    int cRichPartialFragments();
    int gRichPartialFragments();
    int totalTelomericReads();
    double purity();
    double ploidy();
    double duplicateProportion();
    double meanReadDepth();
    double gc50ReadDepth();
}
