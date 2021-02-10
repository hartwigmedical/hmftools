package com.hartwig.hmftools.ckb.interpretation;

import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.Drug;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drugclass.DrugClass;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugsInterpretation {

    @NotNull
    public abstract Drug drug();

    @NotNull
    public abstract DrugClass drugClass();
}
