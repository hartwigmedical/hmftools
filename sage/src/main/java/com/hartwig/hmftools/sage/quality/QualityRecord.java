package com.hartwig.hmftools.sage.quality;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface QualityRecord {
    byte ref();

    byte alt();

    byte qual();

    int position();

    boolean firstOfPair();

    byte[] trinucleotideContext();
}
