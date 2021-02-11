package com.hartwig.hmftools.ckb.interpretation.common.druginterpretation;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.Drug;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugsInterpretation {

    @NotNull
    public abstract List<Drug> drug();
}
