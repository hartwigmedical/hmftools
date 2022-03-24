package com.hartwig.hmftools.common.cobalt;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CobaltRatio extends CobaltCount
{
    double referenceGCRatio();

    double referenceGCDiploidRatio();

    double tumorGCRatio();
}
