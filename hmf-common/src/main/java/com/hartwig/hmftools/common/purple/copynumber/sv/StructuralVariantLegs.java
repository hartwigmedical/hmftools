package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class StructuralVariantLegs {

    public abstract Optional<StructuralVariantLeg> start();

    public abstract Optional<StructuralVariantLeg> end();

}
