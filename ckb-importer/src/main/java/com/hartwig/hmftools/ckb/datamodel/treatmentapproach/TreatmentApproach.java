package com.hartwig.hmftools.ckb.datamodel.treatmentapproach;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatmentApproach {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String profileName();

    @Nullable
    public abstract DrugClassInfo drugClass();

    @Nullable
    public abstract TherapyInfo therapy();

    @NotNull
    public abstract List<ReferenceInfo> reference();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract String updateDate();
}
