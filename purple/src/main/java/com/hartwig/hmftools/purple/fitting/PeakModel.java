package com.hartwig.hmftools.purple.fitting;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PeakModel {

    double peak();

    double bucket();

    double bucketWeight();

    double peakAvgWeight();

    boolean isValid();

    boolean isSubclonal();
}
