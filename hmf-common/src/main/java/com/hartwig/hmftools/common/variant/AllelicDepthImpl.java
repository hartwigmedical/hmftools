package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             visibility = Value.Style.ImplementationVisibility.PUBLIC)
abstract class AllelicDepthImpl implements AllelicDepth
{

}
