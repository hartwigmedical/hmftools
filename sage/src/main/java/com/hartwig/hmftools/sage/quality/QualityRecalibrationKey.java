package com.hartwig.hmftools.sage.quality;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface QualityRecalibrationKey
{

    byte ref();

    byte alt();

    byte[] trinucleotideContext();

    byte qual();
}
