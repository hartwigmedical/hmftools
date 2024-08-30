package com.hartwig.hmftools.common.vchord;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VChordPrediction
{
    double breastCancerHrdScore();
    double ovarianCancerHrdScore();
    double pancreaticCancerScore();
    double prostateCancerScore();
    double otherCancerScore();
}
