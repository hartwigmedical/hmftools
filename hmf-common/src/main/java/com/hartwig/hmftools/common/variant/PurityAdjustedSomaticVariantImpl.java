package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurityAdjustedSomaticVariantImpl implements PurityAdjustedSomaticVariant
{
    public abstract static class Builder implements PurityAdjustedSomaticVariantBuilder
    {
    }
}
