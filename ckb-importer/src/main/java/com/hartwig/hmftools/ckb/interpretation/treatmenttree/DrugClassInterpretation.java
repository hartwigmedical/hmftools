package com.hartwig.hmftools.ckb.interpretation.treatmenttree;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.drugclass.DrugClass;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClassInterpretation {

    @NotNull
    public abstract List<DrugClass> drugClasses();

    @Nullable
    public abstract Drug drug();
}
