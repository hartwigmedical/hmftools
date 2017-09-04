package com.hartwig.hmftools.common.purple.variant;

import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantBuilder;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class PurityAdjustedPurpleSomaticVariant implements PurityAdjustedSomaticVariant {
    public abstract static class Builder implements PurityAdjustedSomaticVariantBuilder {
    }

}
