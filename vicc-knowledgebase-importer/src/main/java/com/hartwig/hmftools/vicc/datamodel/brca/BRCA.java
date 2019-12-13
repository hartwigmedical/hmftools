package com.hartwig.hmftools.vicc.datamodel.brca;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BRCA implements KbSpecificObject {

    @NotNull
    public abstract BRCApart1 brcaPart1();

    @NotNull
    public abstract BRCApart2 brcaPart2();

}
