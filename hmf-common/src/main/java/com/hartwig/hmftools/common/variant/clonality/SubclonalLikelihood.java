package com.hartwig.hmftools.common.variant.clonality;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SubclonalLikelihood {

    double bucket();

    double likelihood();

    double clonalWeight();

    double subclonalWeight();

}
