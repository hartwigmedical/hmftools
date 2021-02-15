package com.hartwig.hmftools.ckb.datamodel.common.druginterpretation;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugInterpretation {

    @NotNull
    public abstract List<Drug> drugs();
}
