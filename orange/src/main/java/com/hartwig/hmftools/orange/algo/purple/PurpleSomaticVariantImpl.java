package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             visibility = Value.Style.ImplementationVisibility.PUBLIC)
public abstract class PurpleSomaticVariantImpl implements PurpleSomaticVariant {
}
