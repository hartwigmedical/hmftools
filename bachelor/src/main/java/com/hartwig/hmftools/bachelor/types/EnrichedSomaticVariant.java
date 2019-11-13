package com.hartwig.hmftools.bachelor.types;

import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedSomaticVariant implements PurityAdjustedSomaticVariant {

}
