package com.hartwig.hmftools.ckb.interpretation.common.druginterpretation;

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

    @Nullable
    public abstract DrugClass drugClass(); //dome drugs has none drugsclass eg. drugs 9758
}
