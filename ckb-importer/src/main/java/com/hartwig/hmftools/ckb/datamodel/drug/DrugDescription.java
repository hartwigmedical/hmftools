package com.hartwig.hmftools.ckb.datamodel.drug;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.Reference;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugDescription {

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<Reference> references();
}
