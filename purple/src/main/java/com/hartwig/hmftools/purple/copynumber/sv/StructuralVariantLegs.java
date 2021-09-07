package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Optional;

import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantLegs {

    public abstract Optional<StructuralVariantLeg> start();

    public abstract Optional<StructuralVariantLeg> end();

}
