package com.hartwig.hmftools.ckb.datamodelinterpretation.treatmentApproch;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.common.ReferenceExtend;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentApproch {

    public abstract int id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract List<ReferenceExtend> references();

    @NotNull
    public abstract Date createDate();

    @NotNull
    public abstract Date updateDate();
}