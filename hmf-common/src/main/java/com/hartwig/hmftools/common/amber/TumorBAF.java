package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface TumorBAF extends GenomePosition {
    @NotNull
    String ref();

    @NotNull
    String alt();

    int normalReadDepth();

    int normalRefSupport();

    int normalAltSupport();

    int tumorReadDepth();

    int tumorRefSupport();

    int tumorAltSupport();

    int tumorAltQuality();
}
